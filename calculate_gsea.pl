#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-09-29 08:04:58 -0700 (Tue, 29 Sep 2009) $';

our ($verbose, $help, $man);
our ($rnkfile, $gmtfile);
our ($seed, $cycle, $setmin, $setmax, $weight, $distance, $mapfile, $logfile, $permfile, $setstatfile,
	$traditional, $skipfdr, $large_es, $pvalue_flag, $leout);


GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'seed=i'=>\$seed, 'cycle|c=i'=>\$cycle, 'setmin=i'=>\$setmin, 'setmax=i'=>\$setmax, 'weight=f'=>\$weight,
	'distance=s'=>\$distance, 'mapfile=s'=>\$mapfile, 'logfile=s'=>\$logfile, 'permfile=s'=>\$permfile, 'traditional'=>\$traditional, 
	'skipfdr'=>\$skipfdr, 'large_es'=>\$large_es, 'pvalue_flag'=>\$pvalue_flag, 'setstatfile=s'=>\$setstatfile, 'leout=s'=>\$leout) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($rnkfile, $gmtfile) = @ARGV;
$permfile or $cycle ||= 1000;
if (defined $cycle) {
	$cycle >= 5 or pod2usage ("Error: the --cyle argument should be at least 5 to produce reasonable results");
}

$setmin ||= 20;
$setmax ||= 200;
defined $weight or $weight = 1;
$distance and $distance =~ s/k$/000/;
$distance and $distance =~ s/m$/000000/;
$distance and $distance =~ m/^\d+$/ || pod2usage ("Error in argument: --distance must be specified as numbers (suffix of k and m are allowable)");
defined $seed or $seed = 1;
srand ($seed);
$traditional and $large_es and pod2usage ("Error in argument: --large_es cannot be used together with --traditional");

main ();

sub main {
	#gene_stat is a hash, where key is gene id, and value is the statistic for the gene from an experiment
	#gene_rank is a hash, where key is gene id, and value is the rank (from large to small) of the gene from all genes
	my ($gene_stat, $gene_rank, $gene_index, $allgenestat, $allgenestatsort, $snp_stat, $snp_pos, $allsnpstatsort);
	my ($geneset, $gene_snp_map, $gene_snp_pos);
	my ($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi);
	my (@cycle_max, @cycle_min, @result);
	
	if ($mapfile) {									#if --mapfile is specified, the items in rnkfile is SNPs, rather than genes
		($snp_stat) = readRnkFile ($rnkfile);
		($gene_snp_map, $gene_snp_pos, $snp_pos, $allsnpstatsort) = readMapFile ($mapfile, $snp_stat);
		($geneset) = readGmtFile ($gmtfile, $gene_snp_pos);
		if (not %$geneset) {
			print STDERR "ERROR: No gene set is found from $gmtfile that match the inclusion criteria. Program exits.\n";
			exit (0);
		}
		
		if ($setstatfile) {							#write the top SNP and its stat, for each gene, for each gene set
			print STDERR "NOTICE: Writting top SNPs for each gene in each gene set to $setstatfile (one gene set per line)\n";
			outputGeneStat ($setstatfile, $geneset, $gene_snp_map, $snp_stat);
		}
		
		($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi) = analyzeAllSnpSet ($geneset, $snp_stat, $snp_pos, $allsnpstatsort, $weight, $gene_snp_pos, $leout, $gene_snp_map);
	} else {
		($gene_stat, $gene_index, $allgenestat) = readRnkFile ($rnkfile);	#gene_index is the index of gene in rankfile (not sorted, but in the same order as in the rnkfile)
		($geneset) = readGmtFile ($gmtfile, $gene_stat);
		if (not %$geneset) {
			print STDERR "ERROR: No gene set is found from $gmtfile that match the inclusion criteria. Program exits.\n";
			exit (0);
		}
		
		my (@gene_sort, $gene_pos, $allgenestatsort);				#gene_pos is the index of gene in SORTED stat values (sort from big to small)
		for my $geneid (keys %$gene_index) {
			push @gene_sort, [$geneid, $gene_stat->{$geneid}];
		}
		@gene_sort = sort {$b->[1] <=> $a->[1]} @gene_sort;
		for my $i (0 .. @gene_sort-1) {
			$gene_pos->{$gene_sort[$i]->[0]} = $i;
		}
		@$allgenestatsort = map {$_->[1]} @gene_sort;				#allgenestatsort is the sorted list of all genes
	
		($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi) = analyzeAllGeneSet ($geneset, $gene_stat, $gene_pos, $allgenestatsort, $weight);
	}

	#calculate max and min for each cycle, for the FWER calculation
	for my $i (0 .. $cycle-1) {
		my ($cycle_max, $cycle_min);						#maximum and minimum (over all gene sets) value of this cycle
		for my $setid (keys %$geneset) {
			defined $cycle_max or $cycle_max = $set_nespi->{$setid}->[$i];
			defined $cycle_min or $cycle_min = $set_nespi->{$setid}->[$i];
			$cycle_max < $set_nespi->{$setid}->[$i] and $cycle_max = $set_nespi->{$setid}->[$i];
			$cycle_min > $set_nespi->{$setid}->[$i] and $cycle_min = $set_nespi->{$setid}->[$i];
		}
		push @cycle_max, $cycle_max;
		push @cycle_min, $cycle_min;
	}
	
	#calculate statistical significance (FDR and FWER) for this gene set
	for my $setid (keys %$geneset) {
		my ($fdr, $fwer);
		if ($skipfdr) {
			$fdr = 'NA'; $fwer = 'NA'; next;
		}
		if ($traditional) {
			$fdr = calculateFDR_TRADITIONAL ($set_nes->{$setid}, $geneset, $set_nes, $set_nespi);
			$fwer = calculateFWER_TRADITIONAL ($set_nes->{$setid}, \@cycle_max, \@cycle_min);
		} else {
			$fdr = calculateFDR1 ($set_nes->{$setid}, $geneset, $set_nes, $set_nespi);
			$fwer = calculateFWER1 ($set_nes->{$setid}, \@cycle_max, \@cycle_min);
		}
		push @result, [$setid, scalar (@{$geneset->{$setid}}), $set_es->{$setid}, $set_nes->{$setid}, $set_nominalp->{$setid}, $fdr, $fwer];
	}
	#print out the formatted results
	outputResult (@result);
	$logfile and outputLog ($logfile, $geneset, $set_es, $set_nes, $set_nominalp, $set_espi);
}

#the ES (rather than NES!) can be recorded to log file so that multiple log files can be combined to give more accurate statistics
sub outputLog {
	my ($logfile, $geneset, $set_es, $set_nes, $set_nominalp, $set_espi) = @_;
	open (LOG, ">$logfile") or confess "Error: cannot write to logfile $logfile: $!";
	print STDERR "NOTICE: The NES and NESpi values are written to logfile $logfile\n";
	for my $setid (keys %$set_nes) {
		print LOG $setid, "\tsize=", scalar (@{$geneset->{$setid}}), "\tES=", sprintf ("%.3f", $set_es->{$setid}), "\tNES=", sprintf ("%.3f", $set_nes->{$setid}), "\tnominalP=", sprintf ("%.3f", $set_nominalp->{$setid}), "\t", join ("\t", map {sprintf("%.3f", $_)} @{$set_espi->{$setid}}), "\n";
	}
	close (LOG);
}	

#given a list of snpstat values, and the gene-snp mapping, generate an array of genestat (maximum value of associated SNPs), and a hash containing the index of each gene in the array
sub convertSnp2Gene {
	my ($allsnpstat, $gene_snp_pos) = @_;
	my ($gene_pos, $allgenestatsort, @temp);
	for my $geneid (keys %$gene_snp_pos) {
		my @snpstat = map {$allsnpstat->[$_]} @{$gene_snp_pos->{$geneid}};
		@snpstat = sort {$b<=>$a} @snpstat;
		push @temp, [$geneid, $snpstat[0]];
	}
	@temp = sort {$b->[1]<=>$a->[1]} @temp;
	@$allgenestatsort = map {$_->[1]} @temp;
	for my $i (0 .. @temp-1) {
		$gene_pos->{$temp[$i]->[0]} = $i;
	}
	return ($gene_pos, $allgenestatsort);
}

sub outputLeadingEdge {
	my ($leout, $geneset, $gene_snp_map, $snp_stat, $le) = @_;
	
	open (LEOUT, ">$leout") or confess "Error: cannot write to Leading Edge output file $leout: $!";
	for my $setid (sort keys %$geneset) {
		print LEOUT $setid, "\t", $setid;
		my @geneid = @{$geneset->{$setid}};
		my @genesort;			#array element1=output element2=stat
		for my $geneid (@geneid) {
			my @snpsort;		#array element1=output element2=stat
			for my $snpid (@{$gene_snp_map->{$geneid}}) {
				push @snpsort, [$snpid, $snp_stat->{$snpid}];
			}
			@snpsort = sort {$a->[1] <=> $b->[1]} @snpsort;		#sort all SNPs for each gene
			push @genesort, [$geneid, $snpsort[$#snpsort]->[1]];
		}
		@genesort = sort {$b->[1] <=> $a->[1]} @genesort;
		
		$le->{$setid} < @genesort or confess "Error: the gene set $setid contains only " . scalar (@genesort). " genes but the leading edge has $le->{$setid} genes";
		for my $i (0 .. $le->{$setid}) {
			print LEOUT "\t$genesort[$i]->[0]";
		}
		
		print LEOUT "\n";
	}
	close (LEOUT);
}

sub outputGeneStat {
	my ($setstatfile, $geneset, $gene_snp_map, $snp_stat) = @_;
	
	open (SETSTAT, ">$setstatfile") or confess "Error: cannot write to gene set statistics file $setstatfile: $!";
	for my $setid (sort keys %$geneset) {
		print SETSTAT $setid;
		my @geneid = @{$geneset->{$setid}};
		my @genesort;			#array element1=output element2=stat
		for my $geneid (@geneid) {
			my @snpsort;		#array element1=output element2=stat
			for my $snpid (@{$gene_snp_map->{$geneid}}) {
				push @snpsort, [$snpid, $snp_stat->{$snpid}];
			}
			@snpsort = sort {$a->[1] <=> $b->[1]} @snpsort;		#sort all SNPs for each gene
			if ($pvalue_flag) {					#make sure that the output P-value is correct by using 1-P formula
				push @genesort, ["$geneid,$snpsort[$#snpsort]->[0],".(2 ** (-$snpsort[$#snpsort]->[1])), $snpsort[$#snpsort]->[1]];
			} else {
				push @genesort, ["$geneid,$snpsort[$#snpsort]->[0],$snpsort[$#snpsort]->[1]", $snpsort[$#snpsort]->[1]];
			}
		}
		@genesort = sort {$b->[1] <=> $a->[1]} @genesort;
		for my $i (0 .. @genesort-1) {
			print SETSTAT "\t$genesort[$i]->[0]";
		}
		
		print SETSTAT "\n";
	}
	close (SETSTAT);
}

sub analyzeAllSnpSet {
	my ($geneset, $snp_stat, $snp_pos, $allsnpstatsort, $weight, $gene_snp_pos, $leout, $gene_snp_map) = @_;
	my ($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi);

	my ($gene_pos, $allgenestatsort) = convertSnp2Gene ($allsnpstatsort, $gene_snp_pos);	#gene_pos: hash containing index of genes in genestatsort
	my %pos_gene = reverse %$gene_pos;
	my (%le);										#leading edge index for each gene set
	for my $setid (sort keys %$geneset) {
		my (@set_posindex, @set_stat);
		@set_posindex = map {$gene_pos->{$_}} @{$geneset->{$setid}};
		@set_stat = @$allgenestatsort[@set_posindex];
		my $es;
		if ($large_es) {
			$es = calculateESFast_LARGE (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
		} else {
			$es = calculateESFast (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
		}
		$set_es->{$setid} = $es;
		$verbose and print STDERR "NOTICE: Set=$setid ES=$es\n";
		if ($setid eq 'DEBUG_GO0043085') {							#for debugging purposes only
			print "setposindex=@set_posindex\nsetstat=@set_stat\n";
			print "set=$setid gene=", join(",", map {$pos_gene{$_}} sort {$a<=>$b} @set_posindex), "\n";
			print "set=$setid genepos=", join(",", sort {$a<=>$b} @set_posindex), "\n";
			print "set=$setid stat=", join(",", sort {$b<=>$a} @set_stat), "\n";
			print "set=$setid topsnp=", join (",", @{$gene_snp_pos->{'M3KL4_HUMAN'}}), "\n";
			print STDERR "NOTICE: set=$setid ES=$es\n";
		}
		
		$leout and $le{$setid} = calculateLeadingEdge (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
	}
	
	if ($leout) {
		print STDERR "NOTICE: Writting Leading Edge output file to $leout as a new gene set file\n";
		outputLeadingEdge ($leout, $geneset, $gene_snp_map, $snp_stat, \%le);
	}
	
	my @allsnpstat = @$allsnpstatsort;
	my $perm_snpstat;									#contains permutated snpstat values for each cycle
	$permfile and $perm_snpstat = readPermFile ($permfile, $snp_pos);
	for my $current_cycle (1 .. $cycle) {
		if ($permfile) {
			@allsnpstat = ();
			for my $i (0 .. @$perm_snpstat-1) {
				if ($pvalue_flag) {
					#push @allsnpstat, 1-(split (/,/, $perm_snpstat->[$i]))[$current_cycle-1];		#2008sep24: this is a bad treatment of the test statistic
					push @allsnpstat, -log ((split (/,/, $perm_snpstat->[$i]))[$current_cycle-1]) / log(2);
				} else {
					push @allsnpstat, (split (/,/, $perm_snpstat->[$i]))[$current_cycle-1];
				}
			}
		} else {
			fisher_yates_shuffle (\@allsnpstat);					#this is NOT recommended, because it disrupts SNP-SNP correlation structure
		}
		my ($gene_pos, $allgenestatsort) = convertSnp2Gene (\@allsnpstat, $gene_snp_pos);
		for my $setid (sort keys %$geneset) {
			my (@set_posindex, @set_stat, $es);
			@set_posindex = sort {$a<=>$b} map {$gene_pos->{$_}} @{$geneset->{$setid}};
			@set_stat = @$allgenestatsort[@set_posindex];
			if ($large_es) {
				$es = calculateESFast_LARGE (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
			} else {
				$es = calculateESFast (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
			}
			push @{$set_espi->{$setid}}, $es;
			if ($setid eq 'DEBUG_GO0043085') {						#for debugging purposes only
				#print STDERR "NOTICE: setid=$setid cycle=$current_cycle ES=$es posindex=@set_posindex stat=@set_stat\n";
				my @genesize = sort {$b<=>$a} map {scalar (@{$gene_snp_pos->{$_}})} @{$geneset->{$setid}};
				#print STDERR "genesize=@genesize\n";
			}
		}
		$current_cycle =~ /0$/ and print STDERR "NOTICE: permutation cycle $current_cycle done!\n";
	}
	
	for my $setid (sort keys %$geneset) {
		my ($nominalp, $nes, $nespi);
		if ($traditional) {
			($nominalp, $nes, $nespi) = calculateNES_TRADITIONAL ($set_es->{$setid}, $set_espi->{$setid});
		} else {
			($nominalp, $nes, $nespi) = calculateNES1 ($set_es->{$setid}, $set_espi->{$setid});
		}
		$set_nominalp->{$setid} = $nominalp;
		$set_nes->{$setid} = $nes;
		$set_nespi->{$setid} = $nespi;
	}
	return ($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi);
}

sub readPermFile {
	my ($permfile, $snp_pos) = @_;
	my ($count_found_snp, @perm_snpstat, @column_match) = (0);
	open (PERM, $permfile) or confess "Error: cannot read from permutation file $permfile: $!";
	print STDERR "NOTICE: Reading permutation file $permfile ... ";
	
	$_ = <PERM>;
	s/\s*[\r\n]+$//;
	my @header = split (/\s+/, $_);
	my ($marker_index, $column);
	for my $i (0 .. @header-1) {
		if ($pvalue_flag) {
			$header[$i] =~ m/_P_PERM$/ and $column = $i;
		} else {
			defined $column and $header[$i] =~ m/(CHI2|F)_PERM$/ and confess "Error: the permutation file $permfile cannot contain multiple types of permutation results ($header[$column] and $header[$i])";
			$header[$i] =~ m/(CHI2|F)_PERM$/ and $column = $i;		#look for the column "CHI2_PERM" or "F_PERM" in the permutation file to retrieve permutated test chi2 statistic values
		}
		$header[$i] eq 'Name' and $marker_index = $i;
		$header[$i] eq 'Marker' and $marker_index = $i;
		$header[$i] eq 'SNP' and $marker_index = $i;
	}
	defined $column or confess "Error: unable to find the column (either CHI2_PERM or P_PERM) containing permutation statistics in the header line of permfile $permfile";
	defined $marker_index or confess "Error: unable to find the column (Name) containing SNP names";

	while (<PERM>) {
		my @record = split (/\s+/, $_);
		exists $snp_pos->{$record[$marker_index]} or next;
		$record[$column] =~ s/^,//;					#sometimes there is a leading comma in the group of test statistic values (comma-separated string)
		push @perm_snpstat, $record[$column];
		push @column_match, [$count_found_snp, $snp_pos->{$record[$marker_index]}];
		$count_found_snp++;
		
		if (not defined $cycle) {
			$cycle = ($record[$column] =~ tr/,/,/) + 1;
			print STDERR "Automatically setting --cycle argument to $cycle ... ";
		} else {
			$cycle == ($record[$column] =~ tr/,/,/) + 1 or confess "Error: discordant number of cycles ($cycle) with the number of permutated statistic in current line ($record[$column])";
		}
	}
	print STDERR "Done with permuted statistic values for $count_found_snp SNPs!\n";
	$count_found_snp == scalar (keys %$snp_pos) or confess "Error: discordant number of SNPs $count_found_snp (in $permfile) and ", scalar (keys %$snp_pos), " (in $mapfile with --dist argument)\n";

	print STDERR "Rearranging SNP stat values ... ";
	@column_match = sort {$a->[1] <=> $b->[1]} @column_match;
	@column_match = map {$_->[0]} @column_match;
	@perm_snpstat = @perm_snpstat[@column_match];
	print STDERR "Done!\n";

	return \@perm_snpstat;
}

sub analyzeAllGeneSet {
	my ($geneset, $gene_stat, $gene_pos, $allgenestatsort, $weight) = @_;		#gene_pos is the index of the gene in allgenestatsort array
	my ($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi);
	
	my %pos_gene = reverse %$gene_pos;
	for my $setid (sort keys %$geneset) {
		my @set_posindex = map {$gene_pos->{$_}} @{$geneset->{$setid}};
		my @set_stat = @$allgenestatsort[@set_posindex];
		my $es;
		if ($large_es) {
			$es = calculateESFast_LARGE (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
		} else {
			$es = calculateESFast (\@set_posindex, \@set_stat, $weight, @$allgenestatsort-@set_stat);
		}
		$set_es->{$setid} = $es;
		$verbose and print STDERR "NOTICE: Set=$setid ES=$es\n";
		if ($setid eq 'DEBUG_GO0043085') {						#for debugging purposes only
			print "set=$setid gene=", join(",", @{$geneset->{$setid}}), "\n";
			print "set=$setid gene=", join(",", map {$pos_gene{$_}} sort {$a<=>$b} @set_posindex), "\n";
			print "set=$setid stat=", join(",", sort {$b<=>$a} @set_stat), "\n";
		}
	}
	print STDERR "NOTICE: calculating all ESpi ... ";
	$set_espi = calculateAllESpi ($geneset, $gene_stat, $gene_pos, $allgenestatsort, $weight);
	print STDERR "Done!\n";
	
	for my $setid (keys %$geneset) {
		my ($nominalp, $nes, $nespi);
		if ($traditional) {
			($nominalp, $nes, $nespi) = calculateNES_TRADITIONAL ($set_es->{$setid}, $set_espi->{$setid});
		} else {
			($nominalp, $nes, $nespi) = calculateNES1 ($set_es->{$setid}, $set_espi->{$setid});
		}
		$set_nominalp->{$setid} = $nominalp;
		$set_nes->{$setid} = $nes;
		$set_nespi->{$setid} = $nespi;
	}
	return ($set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi);
}

sub calculateNES_TRADITIONAL {
	my ($es, $espi) = @_;
	my ($nominalp, $nes, @nespi) = (0);
	my (@pos_espi, @neg_espi, $mean_posespi, $mean_negespi);
	if ($es >= 0) {
		@pos_espi = grep {$_>=0} @$espi;
		@neg_espi = grep {$_<0} @$espi;
		@pos_espi or push (@pos_espi, 0.1) and print STDERR "FATAL ERROR: unable to find positive mean (increase --cycle argument might help) but 0.1 is used here to prevent program exit\n";
		$mean_posespi = mean (\@pos_espi);
		@neg_espi and $mean_negespi = mean (\@neg_espi);
		$nes = $es / $mean_posespi;

		for (@$espi) {
			if ($_ >= 0) {
				push @nespi, $_/$mean_posespi;
			} else {
				push @nespi, -$_/$mean_negespi;
			}
		}
		for (@pos_espi) {
			$_ >= $es and $nominalp++;
		}
		$nominalp /= @pos_espi;
	} else {
		@pos_espi = grep {$_>=0} @$espi;
		@neg_espi = grep {$_<0} @$espi;
		@neg_espi or push (@neg_espi, -0.1) and print STDERR "FATAL ERROR: unable to find negative mean (increase --cycle argument might help) but -0.1 is used here to prevent program exit\n";
		$mean_negespi = mean (\@neg_espi);
		@pos_espi and $mean_posespi = mean (\@pos_espi);
		$nes = -$es / $mean_negespi;

		for (@$espi) {
			if ($_ >= 0) {
				push @nespi, $_/$mean_posespi;
			} else {
				push @nespi, -$_/$mean_negespi;
			}
		}

		for (@neg_espi) {
			$_ <= $es and $nominalp++;
		}
		$nominalp /= @neg_espi;
	}
	return ($nominalp, $nes, \@nespi);
}

sub calculateNES1 {
	my ($es, $espi) = @_;
	my ($nominalp, $nes, @nespi) = (0);
	my $mean_espi = mean ($espi);
	my $sd_espi = sd ($espi);
	$nes = ($es-$mean_espi)/$sd_espi;
	for (@$espi) {
		push @nespi, ($_-$mean_espi)/$sd_espi;
		$_>=$es and $nominalp++;
	}
	$nominalp /= @$espi;
	return ($nominalp, $nes, \@nespi);
}

sub calculateFDR1 {
	my ($current_nes, $geneset, $set_nes, $set_nespi) = @_;
	my ($count_nespi, $count_nespi_flag, $count_nes_flag, $count_num_geneset, $fdr) = (0, 0, 0, 0);
	for my $nextsetid (keys %$geneset) {						#retrieve information for all these gene sets
		for (@{$set_nespi->{$nextsetid}}) {
			$_ >= $current_nes and $count_nespi_flag++;
			$count_nespi++;
		}
	}

	$count_nes_flag = scalar (grep {$_ >= $current_nes} values %$set_nes);
	$count_num_geneset = scalar (values %$set_nes);

	$fdr = $count_nespi_flag / $count_nespi / ($count_nes_flag / $count_num_geneset);
	if ($fdr >= 1) {
		#printf ("WARNING: HUGE FDR: %i %i %i %i $fdr\n", $count_nespi_flag, $count_nespi, $count_nes_flag, $count_num_geneset);
		$fdr = 1;
	}
	return $fdr;
}

#traditional ways to calculate FDR that dichotimize the positive and negative values
sub calculateFDR_TRADITIONAL {
	my ($current_nes, $geneset, $set_nes, $set_nespi) = @_;
	my ($count_nespi, $count_nespi_flag, $count_nes_flag, $count_num_geneset, $fdr) = (0, 0, 0, 0);
	for my $nextsetid (keys %$geneset) {						#retrieve information for all these gene sets
		for (@{$set_nespi->{$nextsetid}}) {
			if ($current_nes >= 0) {					#process positive side
				$_ >= $current_nes and $count_nespi_flag++;
				#$_ >= 0 and $count_nespi++;				#there is no need to dichotomize the count_nespi (2008Oct)
				$count_nespi++;
			} else {
				$_ <= $current_nes and $count_nespi_flag++;
				#$_ < 0 and $count_nespi++;				#there is no need to dichotomize the count_nespi (2008Oct)
				$count_nespi++;
			}
		}
	}
	if ($current_nes >= 0) {
		$count_nes_flag = scalar (grep {$_ >= $current_nes} values %$set_nes);
		$count_num_geneset = scalar (grep {$_>=0} values %$set_nes);
	} else {
		$count_nes_flag = scalar (grep {$_ <= $current_nes} values %$set_nes);
		$count_num_geneset = scalar (grep {$_<0} values %$set_nes);
	}
	if (!$count_nespi or !$count_nes_flag or !$count_num_geneset) {
		print STDERR "WARNING: FDR cannot be calculated so 1 is assumed\n";
		$fdr = 1;
	} else {
		$fdr = $count_nespi_flag / $count_nespi / ($count_nes_flag / $count_num_geneset);
		if ($fdr >= 1) {
			print STDERR "WARNING: HUGE FDR: $count_nespi_flag, $count_nespi, $count_nes_flag, $count_num_geneset\n";
			$fdr = 1;
		}
	}
	return $fdr;
}

sub calculateFWER1 {
	my ($current_nes, $cycle_max, $cycle_min) = @_;
	my (@more_extreme, $fwer);
	@more_extreme = grep {$_ >= $current_nes} @$cycle_max;
	$fwer = @more_extreme/@$cycle_max;
	return $fwer;
}

#traditional FWER calculate dichotomize the NES values so positive and negative values have different method
sub calculateFWER_TRADITIONAL {
	my ($current_nes, $cycle_max, $cycle_min) = @_;
	my (@more_extreme, $fwer);
	my $cycle_pos_max = scalar (grep {$_ >= 0} @$cycle_max);
	my $cycle_neg_min = scalar (grep {$_ <0 } @$cycle_min);
	if ($current_nes >= 0) {
		@more_extreme = grep {$_ >= $current_nes} @$cycle_max;
		$fwer = @more_extreme/$cycle_pos_max;
	} else {
		@more_extreme = grep {$_ <= $current_nes} @$cycle_min;
		$fwer = @more_extreme/$cycle_neg_min;
	}
	return $fwer;
}

sub outputResult {
	my (@result) = @_;
	@result = sort {$b->[3] <=> $a->[3]} @result;
	if ($traditional) {
		print "<-----------------Over-represented in head of ranked list-------------------------->\n";
		for my $result (@result) {
			$result->[3] >= 0 or next;
			print "Geneset=$result->[0]\tSize=", sprintf ("%-4d", $result->[1]), "\tES=", sprintf ("%.3f", $result->[2]), "\tNES=", sprintf ("%.3f", $result->[3]), "\tNominalP=", sprintf ("%.5f", $result->[4]), "\tFDR=", sprintf ("%.3f", $result->[5]), "\tFWER=", sprintf ("%.3f", $result->[6]), "\n";
		}
		print "<-----------------Over-represented in tail of ranked list-------------------------->\n";
		@result = sort {$a->[3] <=> $b->[3]} @result;
		for my $result (@result) {
			$result->[3] < 0 or next;
			print "Geneset=$result->[0]\tSize=", sprintf ("%-4d", $result->[1]), "\tES=", sprintf ("%.3f", $result->[2]), "\tNES=", sprintf ("%.3f", $result->[3]), "\tNominalP=", sprintf ("%.5f", $result->[4]), "\tFDR=", sprintf ("%.3f", $result->[5]), "\tFWER=", sprintf ("%.3f", $result->[6]), "\n";
		}
	} else {
		print "<-----------------Ranked list of over-represented gene sets/pathways-------------------------->\n";
		for my $result (@result) {
			print "Geneset=$result->[0]\tSize=", sprintf ("%-4d", $result->[1]), "\tES=", sprintf ("%.3f", $result->[2]), "\tNES=", sprintf ("%.3f", $result->[3]), "\tNominalP=", sprintf ("%.5f", $result->[4]), "\tFDR=", sprintf ("%.3f", $result->[5]), "\tFWER=", sprintf ("%.3f", $result->[6]), "\n";
		}
	}
}

sub calculateLeadingEdge {
	my ($posindex, $stat, $weight, $num_miss) = @_;
	my @posindex = sort {$a <=> $b} @$posindex;				#double check to make sure that posindex is sorted
	my @stat = sort {$b <=> $a} @$stat;					#double check to make sure that stat is sorted

	my ($current_sum_rj, $n_r, $es, $current_es, $p_hit, $p_miss) = (0);
	my ($leading_edge) = (0);

	$n_r += abs ($_) ** $weight for @stat;		#use absolute value here (sometimes the correlation might be negative)
	for my $i (0 .. @posindex-1) {
		$p_miss = ($posindex[$i] - $i) / $num_miss;
		
		#"maximum deviation from zero" could be achieved by previous gene (especially when ES is negative)
		$p_hit = $current_sum_rj / $n_r;
		$current_es = $p_hit - $p_miss;
		defined $es or $es = $current_es;
		if (abs ($es) < abs ($current_es)) {
			$es = $current_es;
			$leading_edge = $i;
		}
		
		$current_sum_rj += abs ($stat[$i]) ** $weight;
		$p_hit = $current_sum_rj / $n_r;
		$current_es = $p_hit - $p_miss;
		if (abs ($es) < abs ($current_es)) {
			$es = $current_es;
			$leading_edge = $i;
		}
	}
	return $leading_edge;
}
	

#fast calculation of ES. Originally used in the espi calculation, but I decided to merge the subroutine to the espi calculation
sub calculateESFast {
	my ($posindex, $stat, $weight, $num_miss) = @_;
	my @posindex = sort {$a <=> $b} @$posindex;				#double check to make sure that posindex is sorted
	my @stat = sort {$b <=> $a} @$stat;					#double check to make sure that stat is sorted

	my ($current_sum_rj, $n_r, $es, $current_es, $p_hit, $p_miss) = (0);

	$n_r += abs ($_) ** $weight for @stat;		#use absolute value here (sometimes the correlation might be negative)
	for my $i (0 .. @posindex-1) {
		$p_miss = ($posindex[$i] - $i) / $num_miss;
		
		#"maximum deviation from zero" could be achieved by previous gene (especially when ES is negative)
		$p_hit = $current_sum_rj / $n_r;
		$current_es = $p_hit - $p_miss;
		defined $es or $es = $current_es;
		abs ($es) < abs ($current_es) and $es = $current_es;
		
		$current_sum_rj += abs ($stat[$i]) ** $weight;
		$p_hit = $current_sum_rj / $n_r;
		$current_es = $p_hit - $p_miss;
		abs ($es) < abs ($current_es) and $es = $current_es;
	}
	return $es;
}

#fast calculation of largest ES only (i.e., only consider large ES values, as opposed of "maximum deviation from zero")
sub calculateESFast_LARGE {
	my ($posindex, $stat, $weight, $num_miss) = @_;
	my @posindex = sort {$a <=> $b} @$posindex;				#double check to make sure that posindex is sorted
	my @stat = sort {$b <=> $a} @$stat;					#double check to make sure that stat is sorted

	my ($n_r, $es, $current_es, $current_sum_rj, $p_hit, $p_miss);

	$n_r += abs ($_) ** $weight for @stat;					#use absolute value here (sometimes the correlation might be negative)
	for my $i (0 .. @posindex-1) {
		$p_miss = ($posindex[$i] - $i) / $num_miss;
		
		$current_sum_rj += abs ($stat[$i]) ** $weight;
		$p_hit = $current_sum_rj / $n_r;
		$current_es = $p_hit - $p_miss;
		defined $es or $es = $current_es;
		$es < $current_es and $es = $current_es;
	}
	return $es;
}

sub calculateAllESpi {
	my ($geneset, $gene_stat, $gene_pos, $allgenestatsort, $weight) = @_;
	my (@random_allstat, %set_posindex, $set_espi);
	my (@set_posindex, @set_stat, $es);
	
	for my $i (0 .. @$allgenestatsort-1) {
		push @random_allstat, [$allgenestatsort->[$i], $i];		#this array will be shuffled later
	}
	
	for my $setid (keys %$geneset) {
		@{$set_posindex{$setid}} = map {$gene_pos->{$_}} @{$geneset->{$setid}};
	}
	
	for my $current_cycle (1 .. $cycle) {
		fisher_yates_shuffle (\@random_allstat);
		for my $setid (keys %$geneset) {
			@set_posindex = @{$set_posindex{$setid}};
			@set_stat = map {$_->[0]} @random_allstat[@set_posindex];
			@set_posindex = map {$_->[1]} @random_allstat[@set_posindex];
			if ($large_es) {
				$es = calculateESFast_LARGE (\@set_posindex, \@set_stat, $weight, @random_allstat-@set_posindex);
			} else {
				$es = calculateESFast (\@set_posindex, \@set_stat, $weight, @random_allstat-@set_posindex);
			}
			push @{$set_espi->{$setid}}, $es;
		}
	}
	return $set_espi;
}

#read the stat values for genes or SNPs sequentially.
sub readRnkFile {
	my ($rnkfile) = @_;
	my (@gene_stat, %gene_stat, %gene_index, @allstat);
	my ($count_invalid_record, $index) = (0, 0);
	open (RNK, $rnkfile) or confess "Error: cannot read from rnkfile $rnkfile: $!";
	print STDERR "NOITCE: Reading gene/snp-stat-file $rnkfile ... ";
	while (<RNK>) {
		s/[\r\n]+$//;
		m/^(\S+)\t([\d\.\-\+eE]+)/ or ++$count_invalid_record and next;
		my ($gene, $stat) = ($1, $2);
		exists $gene_stat{$gene} and print STDERR "WARNING: The stat values for $gene occur more than once in $rnkfile: replacing old=$gene_stat{$gene} with new=$stat\n";
		if ($pvalue_flag) {
			$stat >=0 and $stat <= 1 or ++$count_invalid_record and next;
			#$stat = 1 - $stat;			#2008sep24: this is not a good treatment of stat, we should try use -log2(P) formula
			$stat = -log($stat)/log(2);
		}
		$gene_stat{$gene} = $stat;
		$gene_index{$gene} = $index++;
		push @allstat, $stat;
	}
	print STDERR "Done with ${\(scalar @allstat)} records ($count_invalid_record records skipped due to unrecognizable format)\n";
	return (\%gene_stat, \%gene_index, \@allstat);
}

sub readGmtFile {
	my ($gmtfile, $genehash) = @_;
	my %geneset;
	open (GMT, $gmtfile) or confess "Error: cannot read from gmtfile $gmtfile: $!";
	print STDERR "NOTICE: Reading GeneSet file $gmtfile ... ";
	while (<GMT>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);
		my $setid = shift @record;
		my $setname = shift @record;
		my (@newrecord, %newrecord);
		if (exists $geneset{$setid}) {
			confess "Error: the geneset $setid occur more than once in gmtfile $gmtfile. Please use unique gene set identifier";
		}
		for my $nextrecord (@record) {
			my $ucnextrecord;
			
			#$ucnextrecord = uc $nextrecord;
			$ucnextrecord = $nextrecord;
			$newrecord{$ucnextrecord} and next;			#in case the same gene occur multiple times in a gene set
			if (exists $genehash->{$ucnextrecord}) {
				push @newrecord, $ucnextrecord;
				$newrecord{$ucnextrecord} = 1;
			}
		}
		if (@newrecord >= $setmin and @newrecord <= $setmax) {
			$geneset{$setid} = \@newrecord;
			$verbose and print STDERR "NOTICE: Processing gene set $setid with ${\(scalar @newrecord)} genes\n";
		} else {
			$verbose and print STDERR "WARNING: Skipping gene set $setid due to having ${\(scalar @newrecord)} genes\n";
		}
	}
	print STDERR "Done with ", scalar (keys %geneset), " gene sets that meet the size criteria (min=$setmin max=$setmax)\n";
	return (\%geneset);
}

sub fisher_yates_shuffle {
	my $array = shift;
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i, $j] = @$array[$j, $i];
	}
}

sub restart_shuffle {
	my $array = shift;
	my $restart = int rand (@$array);
	$restart ||= 1;
	@$array = @$array[$restart .. (@$array-1), 0 .. ($restart-1)];
}	

sub sd {
	my ($score) = @_;
	@$score >= 2 or confess "Error: cannot perform SD calculation due to lack of data";
	my $mean = mean ($score);
	my $sum;
	for my $i (0 .. @$score-1) {
		$sum += ($score->[$i]-$mean)*($score->[$i]-$mean);
	}
	$sum /= (@$score-1);
	return sqrt ($sum);
}

#this subroutine calculate the arithmatic mean of a list of numbers.
sub mean {
	my @score = @{$_[0]};
	my $sum = 0;
	$sum += $_ for (@score);
	return $sum / @score;
}

#read the snp-gene map file, and delete any SNP that does not have gene associated
sub readMapFile {
	my ($mapfile, $snp_stat) = @_;
	my (%found_snp, %gene_snp_map, %gene_snp_pos);
	open (MAP, $mapfile) or confess "Error: unable to read from mapfile $mapfile: $!";
	print STDERR "NOTICE: Reading snp-gene-map file $mapfile ... ";
	while (<MAP>) {
		s/[\r\n]+$//;
		my @record = split (/\t/, $_);			#the tab-delimited fields are: snp, gene, gene annotation, distance, etc.
		$record[2] =~ m/^\d+$/ or next;			#the distance is not a number (possbily the header line, or possibly NOT_FOUND, or other mis-annotation)
		my ($snpid, $geneid, $snpgenedist) = @record[0, 1, 2];
		defined $distance and defined $snpgenedist and $snpgenedist > $distance and next;
		$geneid = uc $geneid;				#use capital letters for gene identifier (since most analysis was done for human association study)
		exists $snp_stat->{$snpid} or next;		#this SNP is not in the rnkfile so should be ignored
		
		my @subgene = split (/,/, $geneid);		#in case geneid contains several overlapping genes
		for my $subgene (@subgene) {
			push @{$gene_snp_map{$subgene}}, $snpid;
		}
		$found_snp{$snpid}++;
	}
	close (MAP);
	
	#update the snp_stat, snp_rank and allsnpstatsort variables, since now we are only focusing on these particular SNPs (that can be assigned to genes)
	my (@snp_stat, %snp_pos, @allsnpstatsort);
	my ($count_delete_snp) = (0);
	for my $snpid (keys %$snp_stat) {
		if ($found_snp{$snpid}) {
			push @snp_stat, [$snpid, $snp_stat->{$snpid}];
		} else {
			delete $snp_stat->{$snpid};		#SNP from rankfile not found in map file, so do not use this snp
			$count_delete_snp++;
		}
	}
	@snp_stat = sort {$b->[1] <=> $a->[1]} @snp_stat;	#sort from larger values to smaller values
	@allsnpstatsort = map {$_->[1]} @snp_stat;
	for my $i (0 .. @snp_stat-1) {
		$snp_pos{$snp_stat[$i]->[0]} = $i;
	}
	for my $geneid (keys %gene_snp_map) {
		@{$gene_snp_pos{$geneid}} = map {$snp_pos{$_}} @{$gene_snp_map{$geneid}};
	}
	print STDERR "Done with ", scalar (keys %gene_snp_pos), " genes and ${\(scalar @allsnpstatsort)} SNPs ($count_delete_snp SNPs droppped due to lack of gene mapping)\n";
	return (\%gene_snp_map, \%gene_snp_pos, \%snp_pos, \@allsnpstatsort);
}

=head1 SYNOPSIS

 calculate_gsea.pl [arguments] <snp/gene-stat-file> <gene-set-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --seed <int>		randomization seed
            --cycle <int>		cycle of permutation (default=auto_detect)
            --setmin <int>		minimum number of genes in a gene set to be considered (default=20)
            --setmax <int>		maximum number of genes in a gene set to be considered (default=200)
            --weight <float>		the weighting parameter p in calculating ES score (default=1)
            --distance <int>		maximum distance between SNP and gene for their association in --mapfile
            --mapfile <file>		a file that contains SNP-gene mapping
            --permfile <file>		a file containing test statistic values for permutations
            --logfile <file>		write NES values to this file
            --setstatfile <file>	write top SNP for each gene for each gene set to this file
            --leout <file>		write new gene-set-file with leading edge genes in the original gene-set-file
            
            --large_es			calculate largest ES (may improve GWAS analysis)
            --traditional		traditional method that dichotimize ES (useful in microarray analysis)
            --pvalue_flag		flag that snp/gene-stat-file contains P-values (-log2(P) transformation is used)
            --skip_fdr			skip FDR and FWER calculation (program finishes faster)

 Function: calculate gene set enrichment statistics from the statistics (such as 
 P-value or chi2 statistic) for all genes and a given gene-set file.

 Example: calculate_gsea.pl temp.assoc temp.gmt -map temp.snpgenemap -c 10000 -log temp.log
          calculate_gsea.pl gsea.chi2 gsea.gmt -map gsea.snpgenemap -perm gsea.cc10

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--seed>

specify the randomization seed (default=1). Two instances of the program should 
produce the same results if applied on same data set with same arguments 
(including -seed).

=item B<--cycle>

specify the number of permutations used to calculate ES, NES, FDR and FWER 
statistics.

=item B<--setmin>

a threshold (default=20) specifying the minimum number of genes in a gene set 
that are also observed in the gene-stat-file. Smaller gene sets are not 
considered by this program.

=item B<--setmax>

a threshold (default=200) specifying the maximum number of genes in a gene set 
that are also observed in the gene-stat-file. Larger gene sets are not 
considered by this program.

=item B<--weight>

a parameter used in enrichment score calculation (default=1). When this argument 
is zero, the GSEA method reduces to a Kolmogorov-Smirnov statistic. The original 
authors have found that weight=1 works well for gene expression data sets. In a 
sense, this parameter can be considered as an inflation factor that inflates the 
contribution of genes at the extreme (head or tail) of distributions.

=item B<--mapfile>

a tab-delimited text file containing SNP to gene mapping. The first 4 columns 
are SNP id, gene id, gene description and snp-gene distance. Only the first 2 
columns are manditory.

=item B<--distance>

a threshold specifying the maximum distance allowable to establish association 
between a SNP and a gene.

=item B<--mapfile>

a file that contains SNP-gene mapping and, optionally, gene description and SNP-
gene distance. When SNP-gene distance is given, the --distance argument can be 
used to further filter out SNP and genes that are far apart.

=item B<--logfile>

write NES and NESpi information to this file. This is important in cases where 
multiple runs of the program are used (for example, 1000 cycles are used in each 
of 10 CPUs in a computer cluster), and their results can be combined together to 
generate more accurate FDR and FWER values. The combine_gsea_output.pl program 
to combine these log files together. This way the program can be parallelized, 
which is essential for large data sets such as genome-wide association studies.

=item B<--setstatfile>

write the SNP identifier and its statistic values for each gene within each gene 
set to this file. In the output file, each gene set occupies one line, with tab-
delimited records for each gene, where each record contains gene ID, top SNP ID 
and the statistic value separated by comma. This file is useful for examining 
the top SNPs for each pathway/gene set.

=item B<--large_es>

calculate largest ES values for the gene set. By default this program calculate 
the ES with maximum absolute value ("maximum deviation from zero"), which is 
used in Subramanian et al. However, in many cases when one care only about 
enrichment in one direction, the calculation of largest ES makes more sense.

=item B<--traditional>

apply the same dichotimization techniques used in the original GSEA publication.

=item B<--pvalue_flag>

flag that the snp/gene-stat-file contains P-values, rather than raw test 
statistic values (such as t test statistic, or chi2 test static). When reading 
the stat-file, a formula (-log2(Pvalue)) will be used to convert the P-values so that 
higher value indicate higher significance. However, it is always recommended to 
use raw test statistic, rather than P-values for the GSEA calculation.

=back

=head1 DESCRIPTION

This program is used to calculate Gene Set Enrichment Analysis (GSEA) statistic 
and extended GSEA statistic. The input file <snp/gene-stat-file> contains the 
statistics values (higher means more significant, unless --pvalue_flag argument is 
set) for all genes or all SNPs (when --mapfile is specified) in an experiment. 
The input file <gene-set-file> contains pre-defined sets of genes. The --mapfile 
contains the mapping from SNP to its corresponding genes.

There are two ways the GSEA can be normalized. The first way, or the traditional 
way in Subramanian et al, is by dividing ES score for each gene set by the ES 
scores calculated from all permutations for this gene set (ES/mean(ES_perm)). 
The second way, or the default way in Wang et al, is the Z-score (ES-
mean(ES_perm))/sd(ES_perm). An earlier version of this program was used in the 
paper: B<Wang K, Li M, Bucan M.> I<Pathway-based approaches for analysis of 
genome- wide association studies.> American Journal of Human Genetics, 81:1278-
1283, 2007

Below is a brief description of the application of this program on several 
different scenarios:

=over 8

=item B<Gene expression analysis>

When applied on statistics on genes (such as t-test statistic or P-values for 
differential expression), this program does the same thing as the GSEA program 
(implemented in Java and R), and they should produce highly similar results.

A group of genes can all show higher or lower expression levels in the 
comparison of two physiological conditions. However, due to the use of a 
weighting parameter (as an exponent), the significance levels are calculated 
separately for the positively and negatively scoring gene sets via a permutation 
approach.

Example: 

	calculate_gsea.pl gsea.rnk gsea.gmt -c 100 -log ~/temp

This command calculates gene set enrichment using the pre-ranked module. The 
temp.rnk file contains 1 minus P-values for gene differential expression. The --
log argument specify that the log information be written to ~/temp file. This 
file is important, because you can use combine_gsea.pl program to re-calculate 
enrichment statistics using this log file, or using multiple log files together 
(these multiple files can be generated by specifying differnet randomization 
seeds).

The --traditional argument is generally recommended to be used in gene 
expression analysis. It applies the same dichotimization techniques used in 
Subramanian et al, so that the enrichment statistics is calculated for 
positively and negatively enriched gene sets separately.

=item B<Genome-wide association analysis>

To use this software on genome-wide association (GWA) analysis, several input 
files are necessary: an inputfile containing statistic values for all SNPs, a 
map file containing SNP-gene mapping and the distances between SNPs and genes, a 
gene-set file containing sets of genes to be tested, a permutation file 
containing the permuted test statistic values for all SNPs using phenotype 
permutation.

It is always recommended to use the permutation file that contains real test 
statistic values generated by phenotype permutation. The SNPs in genome is 
highly correlated with each other, so any permutation that disrupt this 
correlation structure will be biased. For this reason, I recommend always using 
the permutation file, rather than permutating statistics between SNPs. The 
permutation file can be generated from the B<calculate_association.pl> program that 
I wrote for GWA analysis. A typical permutation file is shown below (10 lines 
only):

	Marker       Chr   Position          A:B   case_AF  control_AF    CHI2    CHI2_P CHI2_PERM      CHI2_P_PERM
	rs2977670      1     763754          C:G   0.02988   0.01245     3.386     0.184        ,2.495,1.518,1.499,2.567,1.053,0.9854,5.399,1.518,0.9991,6.984  ,0.2872,0.4682,0.4725,0.2771,0.5907,0.611,0.06725,0.4682,0.6068,0.03044
	rs3934834      1    1045729          C:T    0.8669    0.8523     2.628    0.2687        ,1.237,2.286,3.006,2.212,1.365,2.322,6.413,1.143,4.763,0.9403   ,0.5388,0.3189,0.2225,0.3309,0.5053,0.3132,0.04049,0.5647,0.09239,0.6249
	rs3766193      1    1057093          C:G    0.5685    0.5775     1.573    0.4554        ,0.1177,1.826,0.1623,1.696,1.341,1.842,2.355,0.2308,4.827,3.819 ,0.9428,0.4014,0.9221,0.4283,0.5114,0.3982,0.3081,0.891,0.08949,0.1482
	rs3737728      1    1061338          C:T    0.7491    0.7362    0.5062    0.7764        ,0.8074,2.252,0.6411,0.6818,2.305,0.5168,0.9075,3.028,2.964,0.3785      ,0.6678,0.3244,0.7257,0.7111,0.3158,0.7723,0.6352,0.22,0.2272,0.8276
	rs6687776      1    1070488          C:T    0.8352    0.8376    0.8587    0.6509        ,1.643,2.162,1.336,5.999,1.039,1.994,0.6144,6.577,0.1056,0.307  ,0.4399,0.3393,0.5128,0.0498,0.5949,0.369,0.7355,0.0373,0.9486,0.8577
	rs6678318      1    1070556          A:G    0.1648    0.1611    0.8016    0.6698        ,1.358,1.973,1.246,6.602,0.9562,2.009,0.7184,6.66,0.09339,0.4395        ,0.5072,0.3729,0.5364,0.03684,0.62,0.3663,0.6982,0.03579,0.9544,0.8027
	rs9651273      1    1071463          A:G    0.2519    0.2889     2.178    0.3365        ,0.271,0.184,13.01,1.524,0.04583,0.3561,1.945,2.916,1.081,0.2136        ,0.8733,0.9121,0.001493,0.4667,0.9773,0.8369,0.3782,0.2327,0.5825,0.8987
	rs4970405      1    1088878          A:G    0.8848    0.8875    0.5791    0.7486        ,2.743,2.609,0.5305,5.471,2.086,2.068,2.609,2.117,0.8489,0.9476 ,0.2538,0.2713,0.767,0.06486,0.3523,0.3556,0.2713,0.3469,0.6541,0.6226
	rs12726255     1    1089873          A:G    0.8519    0.8487     5.813   0.05467        ,0.503,2.784,2.892,9.033,0.7313,1.991,0.4316,2.114,0.2846,0.2335        ,0.7776,0.2486,0.2355,0.01093,0.6938,0.3695,0.8059,0.3475,0.8674,0.8898

The first line of the permutation file is referred to as the "header" line, 
which contains information on the actual meaning of each column. This program 
will automatically search for a column called CHI2_PERM, and then read in data 
for this column in subsequent lines. For example, for the SNP rs2977670, 10 chi2 
values will be read into memory, starting from 2.495 and ending at 6.984.

Although the above permutation file is generated by the GWA program that I wrote 
(calculate_association.pl), you can also generate one by yourself. For example, 
you can run the PLINK program 1000 times, each time using a differnet phenotype 
definition file, and then collect the chi2 values for each marker at each 
permutation, and then generate a permutation file by yourself with similar 
format as above (basically making sure the Marker column and the CHI2_PERM 
column are present in the permutaion file, since all other columns will be 
ignored by this program)

To run the program, you need to supply a SNP-statistics file (which contains 
chi2 values in the absence of any permutation). In the example above, this file 
can simply be generated by taking the first and the seventh columns in the 
permutation file. If you use your own GWA analysis program, basically you can 
generate a two-column file with one marker per line, and each line 
containsmarker ID and chi2 statistics separated by tab character.

To run the program, you also need to provide a gene set file. This can be 
downloaded from the GSEA website using a series of criteria, but generally the 
file contains highly heterogeneous data. Alternatively, I supply a few GMT files 
in the software distribution, since I believe it has better quality and 
annotation.

A SNP-gene mapping file should be also supplied. This file associate each SNP to 
one gene (occasionally several genes), and gives information on the SNP-gene 
distances. The file can be generated using various gene annotations, such as 
those from ENSEMBL database, those from UCSC Known Gene, or those from RefSeq 
genes.

There are several example files included in the software distribution. You can run:

	calculate_gsea.pl gsea.chi2 gsea.gmt -c 10 -map gsea.snpgenemap -perm gsea.cc10

This command calculates gene set enrichment for SNP statistics contained with 
the gsea.assoc file. The gsea.map file contains gene-SNP mapping. The gsea.cc10 
file contains test statistics for 10 permutation cycles that permutes disease 
labels.

In some rare cases the permutation file is not available. Therefore, this 
software support permutating the statistic values among all SNPs:

	calculate_gsea.pl gsea.chi2 gsea.gmt -c 10 -map gsea.snpgenemap

In the absence of --permfile argument, the program automatically switch the test 
statistic values for all SNPs during each permutation cycle. Again, this is not 
recommended in any way. It is implemented in the program for completeness.

The --large_es argument is generally recommended to be used in GWA analysis. 
It calculate only largest ES values, rather than ES with largest absolute values 
("maximum deviation from zero"). This is because in GWA analysis, we only care 
about gene sets that are enriched at the top of the enrichment scores. Gene sets 
that rank low are simply due to random chance, not due to biological reasons (in 
comparison, in gene expression analysis, gene sets that rank low might be due to 
negative regulation of all genes in the set).

=back

=cut                                                                                                                                                                                                                                                                                
