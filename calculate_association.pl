#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

BEGIN {($_=$0)=~s{[^\\\/]+$}{};$_||="./"}
use lib $_, $_."kext";
use GWA;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-07-10 12:30:41 -0700 (Fri, 10 Jul 2009) $';

our ($verbose, $help, $man);
our ($pedheaderfile, $gtfile);
our ($qt, $tdt, $tsp, $pdt, $remove, $keep, $exclude, $extract, $output, $mind_threshold, $geno_threshold, $maf_threshold, $mme_threshold, $fme_threshold, $ime_threshold, 
	$hwe_threshold, $perfam, $force_tdt, $alleleAB, $cc, $cycle, $seed, $cellsize, $canonical_gtfile_format, $illumina_gtfile_format, $affymetrix_gtfile_format, 
	$plink_tpedfile_format, $simple_format,
	$marker_start, $marker_end, $snppropfile, $allmarker, $allind, $flush_output, $perm_method);
our ($orig_command_line);

processArgument ();
main ();


sub processArgument {
	$orig_command_line = $0 . " @ARGV";
	GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'qt'=>\$qt, 'tdt'=>\$tdt, 'tsp'=>\$tsp, 'pdt'=>\$pdt, 'cc|assoc'=>\$cc, 'force_tdt'=>\$force_tdt, 
		'remove=s@'=>\$remove, 'keep=s@'=>\$keep, 'exclude=s@'=>\$exclude, 'extract=s@'=>\$extract, 'output=s'=>\$output,
		'allmarker'=>\$allmarker, 'allind'=>\$allind, 'mind_threshold=f'=>\$mind_threshold,
		'geno_threshold=f'=>\$geno_threshold, 'maf_threshold=f'=>\$maf_threshold, 'mme_threshold=f'=>\$mme_threshold,
		'fme_threshold=f'=>\$fme_threshold, 'ime_threshold=f'=>\$ime_threshold, 'hwe_threshold=f'=>\$hwe_threshold,
		'perfam'=>\$perfam, 'alleleAB|ab'=>\$alleleAB, 'cycle=i'=>\$cycle, 'seed=i'=>\$seed,
		'cellsize=i'=>\$cellsize, 'canonical_gtfile_format'=>\$canonical_gtfile_format, 'illumina_gtfile_format'=>\$illumina_gtfile_format,
		'affymetrix_gtfile_format'=>\$affymetrix_gtfile_format, 'plink_tpedfile_format'=>\$plink_tpedfile_format, 'simple_format'=>\$simple_format,
		'mstart=i'=>\$marker_start, 'mend=i'=>\$marker_end,
		'snppropfile=s'=>\$snppropfile, 'flush!'=>\$flush_output, 'perm_method=i'=>\$perm_method) or pod2usage ();

	$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
	$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
	@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
	@ARGV == 2 or pod2usage ("Syntax error");
	
	if (not defined $flush_output or $flush_output == 1) {
		$| = 1;							#by default force a flush after every read or write to minotor program progress in real-time (it however increase disk overhead when running the program in parallel in multiple machines so --noflush can be used)
	}
	
	($pedheaderfile, $gtfile) = @ARGV;

	if ($allmarker) {
		defined $geno_threshold || defined $maf_threshold || defined $mme_threshold || defined $hwe_threshold and pod2usage ("Error in argument: please do not specify any marker exclusion criteria when --allmarker is specified");
		$geno_threshold = 1;
		$maf_threshold = 0;
		$mme_threshold = 1;
		$hwe_threshold = 0;
		print STDERR "NOTICE: the --allmarker argument set all marker exclusion criteria as: --geno 1 --maf 0 --mme 0 --hwe 0\n";
	}
	if ($allind) {
		defined $mind_threshold || defined $fme_threshold || defined $ime_threshold and pod2usage ("Error in argument: please do not specify any sample exclusion criteria when the --allsample argument is specified");
		$mind_threshold = 1;
		$fme_threshold = 1;
		$ime_threshold = 1;
		print STDERR "NOTICE: the --allind argument set all individual/sample exclusion criteria as: --mind 1 --fme 1 --ime 1\n";
	}

	defined $mind_threshold or $mind_threshold = 0.1;		#missingness (non-call rate) for individual (same as PLINK)
	defined $geno_threshold or $geno_threshold = 0.1;		#genotyping call rate threshold for SNP (same as PLINK)
	defined $maf_threshold or $maf_threshold = 0.01;		#minor allele frequency threshold for SNP (same as PLINK)
	defined $mme_threshold or $mme_threshold = 0.1;			#marker mendelian error threshold for SNP (same as PLINK)
	defined $fme_threshold or $fme_threshold = 1;			#(OPTION NOT WELL-DEFINED YET) #family mendelian error threshold, where one error in one individual is counted as one error for the family (different as PLINK)
	defined $ime_threshold or $ime_threshold = 0.02;		#individual (offspring) mendelian error threshold (for 500K marker, this equals to less than 10K mendel errors for an offspring)
	defined $hwe_threshold or $hwe_threshold = 0.001;		#hardy-weinburg equilibrium threshold for SNP with exact test (same as PLINK)

	#prepare various output file handles
	if (defined $output) {
		open (STDOUT, ">$output") or confess "Error: cannot write to outputfile $output: $!";
		eval {
			open(STDERR, " | tee $output.log 1>&2");
		};
		if ($@) {
			print STDERR "WARNING: Failed to redirect standard error to log file (Error message will be printed in STDERR only but not recorded in log file)\n";
		} else {
			print STDERR "NOTICE: program notification/warning messages that appear in STDERR will be also recorded in log file $output.log\n";
		}
		open (MINFO, ">$output.minfo") or confess "Error: cannot write to marker information file $output.minfo: $!";
		print STDERR "NOTICE: detailed genotyping summary information for all markers (even if not passing quality control threshold) will be written to $output.minfo\n";
		outputMarkerInfo ('Name', 'Chr', 'Position', 'genotype_missing', 'minor_allele_frequency', 'hwe_p_value', 'marker_mendel_error', 'comments');
		open (SINFO, ">$output.sinfo") or confess "Error: cannot write to sample information file $output.sinfo: $!";
		print SINFO "family_id\tindividual_id\tnum_analyzed_marker\tnocall_count\tnocall_rate\tmendelian_error_count\tmendelian_error_rate\n";
		print STDERR "NOTICE: detailed genotyping summary information for all samples will be written to $output.sinfo\n";
		open (FINFO, ">$output.finfo") or confess "Error: cannot write to family information file $output.finfo: $!";
		print FINFO "family_id\tnum_analyzed_marker\tmendelian_error_count\tmendelian_error_rate\n";
		print STDERR "NOTICE: detailed genotyping summary information for all families with Mendelian errors will be written to $output.finfo\n";

		if ($perfam) {
			$tdt or $tsp or pod2usage ("Error: the --evidence argument is supported only for --tdt and --tsp analysis");
			open (EVI, ">$output.evidence") or confess "Error: cannot write to evidence file $output.evidence";
			print STDERR "NOTICE: detailed statistical evidence on association for each family will be written to $output.evidence\n";
			print EVI '';		#prevent Perl from complaining "file handle used only once"
		}
	} else {
		$perfam and pod2usage ("Error in argument: please specify --output argument when using --evidence argument to record detailed association evidence");
	}



	$maf_threshold >= 0 and $maf_threshold <= 1 or pod2usage ("Error in argument: the --maf_threshold argument must be between 0 and 1 inclusive");
	$hwe_threshold >= 0 and $hwe_threshold <= 1 or pod2usage ("Error in argument: the --hwe_threshold argument must be between 0 and 1 inclusive");
	$mme_threshold >= 0 and $mme_threshold <= 1 or pod2usage ("Error in argument: the --mme_threshold argument must be between 0 and 1 inclusive");
	$geno_threshold >= 0 and $geno_threshold <= 1 or pod2usage ("Error in argument: the --geno_threshold argument must be between 0 and 1 inclusive");
	$fme_threshold >= 0 and $fme_threshold <= 1 or pod2usage ("Error in argument: the --fme_threshold argument must be between 0 and 1 inclusive");
	$ime_threshold >= 0 and $ime_threshold <= 1 or pod2usage ("Error in argument: the --ime_threshold argument must be between 0 and 1 inclusive");
	$mind_threshold >= 0 and $mind_threshold <= 1 or pod2usage ("Error in argument: the --mind_threshold argument must be between 0 and 1 inclusive");
	
	defined $marker_start and $marker_start >= 1 || pod2usage ("Error in argument: the --mstart argument must be greater or equal to 1 (firts marker)");
	$marker_start ||= 1;
	defined $marker_end and $marker_end >= $marker_start || pod2usage ("Error in argument: the --mend argument must be greater than or equal to --mstart argument");

	$allmarker or print STDERR "NOTICE: the following marker exclusion criteria is in effect: --geno $geno_threshold --maf $maf_threshold --mme $mme_threshold --hwe $hwe_threshold\n";
	$allind or print STDERR "NOTICE: the following individual exclusion criteria will *NOT* be used in association test\n        However, individuals failing these criteria will be recorded for re-running program (with --remove argument): --mind $mind_threshold --fme $fme_threshold --ime $ime_threshold\n";
	
	defined $cycle or $cycle = 0;					#by default no permutation cycle is used
	$cycle >= 0 or pod2usage ("Error in argument: --cycle must be a positive integer or zero");
	defined $seed or $seed = 1;					#default randomization seed is 1
	defined $cellsize or $cellsize = 5;				#default minimum cellsize (in contingency table) is five (same as PLINK)
	$cellsize >= 0 or pod2usage ("Error in argument: --cellsize must be a positive integer or zero");
	
	my $count_format = 0;
	$illumina_gtfile_format and $count_format++;
	$affymetrix_gtfile_format and $count_format++;
	$plink_tpedfile_format and $count_format++;
	$canonical_gtfile_format and $count_format++;
	$simple_format and $count_format++;
	
	if ($count_format > 1) {
		pod2usage ("Error in argument: please specify only one argument from --illumina_gtfile_format, --affymetrix_gtfile_format, --plink_tpedfile_format, --simple_format and --canonical_gtfile_format");
	} elsif ($count_format == 0) {
		print STDERR "NOTICE: Assuming the input data is in canonical gtfile format, with columns as SNP, Chr, Position and genotypes\n";
		$canonical_gtfile_format = 1;
	}
	
	my $analysis = 0;
	$tdt and $analysis++;
	$tsp and $analysis++;
	$pdt and $analysis++;
	$cc and $analysis++;
	$qt and $analysis++;
	$analysis > 1 and pod2usage ("Error in argument: please specify only one of the --tdt, --tsp, --pdt or --cc analysis options");
	if (not $analysis) {
		$cc = 1;						#default is to perform case-control association tests
		print STDERR "NOTICE: The default analysis type is case-control comparisons\n";
	}
	$force_tdt and $tdt || pod2usage ("Error in argument: the --force_tdt argument can only be used in conjunction with --tdt argument");
	
	if ($perfam) {
		$tdt or $tsp or pod2usage ("Error in argument: the --evidence argument is only applicable for --tdt or --tsp operation to record per-family transmission ratio");
	}
	
	$perm_method ||= 1;
	$perm_method >= 1 and $perm_method <= 5 or pod2usage ("Error in argument: --perm_method should be an integer between 1-5");
}

sub main {
	#reading ped_header file to make sure that this file is specified correctly, before doing any analysis
	my ($ped, $ped_index) = GWA::readPedHeaderFile ($pedheaderfile, $qt);
	my $ped_stat;
	my (%extract_snp, %exclude_snp);

	#read snppropfile
	my ($snpprop, $proplen);
	$snppropfile and ($snpprop, $proplen) = readSNPPropFile ($snppropfile);
	


	#handle the --keep, --remove, --extract and --exclude argument (they specify the initial inclusion/exclusion criteria for samples and markers)
	if ($keep) {
		my @keepfile = split (/,/, join (',', @$keep));
		my ($newped, %keep_not_found);
		my $count_keep = 0;
		for my $nextfile (@keepfile) {
			open (IND, $nextfile) or confess "Error: cannot read from file $nextfile: $!";
			while (<IND>) {
				s/^\s*|\s*[\r\n]+$//g;			#delete heading and trailing spaces and return characters
				m/\S/ or next;				#skip empty lines
				m/^(\S+)\s+(\S+)$/ or confess "Error: invalid record in $nextfile (two non-blank fields expected): <$_>";
				if ($ped->{$1}{$2}) {
					$newped->{$1}{$2} and next;	#duplicate sample identifiers found in the keep files
					$newped->{$1}{$2} = $ped->{$1}{$2};
					$count_keep++;
				} else {
					$keep_not_found{"$1:$2"}++;
				}
			}
			close (IND);
		}
		$newped or print STDERR "FATAL ERROR: no individual found in ped-header-file $pedheaderfile based on the --keep argument (@keepfile)\n" and exit (1000);
		if (%keep_not_found) {
			my @example = keys %keep_not_found;
			@example = splice (@example, 0, 5);
			print STDERR "WARNING: ${\(scalar keys %keep_not_found)} individuals (for example: @example) specified by --keep argument are not found in pedheader file $pedheaderfile\n";
		}
		$ped = $newped;
		print STDERR "NOTICE: Finished keeping $count_keep individuals in analysis based on --keep argument (@keepfile)\n";
	}
	if ($remove) {
		my @removefile = split (/,/, join (',', @$remove));
		my (%removed_ind, %remove_not_found);
		for my $nextfile (@removefile) {
			open (IND, $nextfile) or confess "Error: cannot read from file $nextfile: $!";
			while (<IND>) {
				s/^\s*|\s*[\r\n]+$//g;			#delete heading and trailing spaces and return characters
				m/\S/ or next;				#skip empty lines
				m/^(\S+)\s+(\S+)$/ or confess "Error: invalid record in $remove (two fields expected): <$_>";
				if ($ped->{$1}{$2}) {
					delete $ped->{$1}{$2};
					$removed_ind{$1, $2}++;
				} elsif ($removed_ind{$1, $2}) {
					1;
				} else {
					$remove_not_found{"$1:$2"}++;
				}
			}
			close (IND);
		}
		if (%remove_not_found) {
			my @example = keys %remove_not_found;
			@example = splice (@example, 0, 5);
			print STDERR "WARNING: ${\(scalar keys %remove_not_found)} individuals (for example: @example) specified by --remove argument are not found in pedheader file $pedheaderfile\n";
		}
		print STDERR "NOTICE: Finished removing ${\(scalar keys %removed_ind)} individuals from analysis based on --remove argument (@removefile)\n";
	}
	if ($extract) {
		my @extractfile = split (/,/, join (',', @$extract));
		for my $nextfile (@extractfile) {
			open (EXTRACT, $nextfile) or confess "Error: cannot read from file $nextfile: $!";
			while (<EXTRACT>) {
				s/^\s*|\s*[\r\n]+$//g;
				m/\S/ or next;			#skip empty lines
				m/^(\S+)$/ or confess "Error: invalid record found in $nextfile (one marker name expected per line): <$_>";
				$extract_snp{$1}++;
			}
			close (EXTRACT);
		}
		print STDERR "NOTICE: Finished reading ${\(scalar keys %extract_snp)} markers to use in analysis based on --extract argument (@extractfile)\n";
	}
	if ($exclude) {
		my @excludefile = split (/,/, join (',', @$exclude));
		for my $nextfile (@excludefile) {
			open (EXCLUDE, $nextfile) or confess "Error: cannot read from file $nextfile: $!";
			while (<EXCLUDE>) {
				s/^\s*|\s*[\r\n]+$//;
				m/\S/ or next;			#skip empty lines
				m/^(\S+)$/ or confess "Error: invalid record found in $nextfile (one marker name expected per line): <$_>";
				$exclude_snp{$1}++;
			}
			close (EXCLUDE);
		}
		print STDERR "NOTICE: Finished reading ${\(scalar keys %exclude_snp)} markers to exclude from analysis based on --exclude argument (@excludefile)\n";
	}
	
	#inspect the pedigree, analyze family structure and calculate certain statistics
	$ped_stat = GWA::inspectPedigree ($ped, $qt);
	
	#perform actual association analysis
	my ($fam_mendel_error, $ind_mendel_error, $ind_missing, $num_analyzed_marker, $nf_index) = calculateAssociation ($gtfile, $ped, $ped_index, $ped_stat, \%extract_snp, \%exclude_snp, $snpprop);
	
	#output sample statistics (missing genotypes and mendelian errors for offsprings)
	outputSampleStat ($ped_index, $fam_mendel_error, $ind_mendel_error, $ind_missing, $num_analyzed_marker);

	#identifying individuals failing the inclusion criteria and print them out (to be included in --remove argument to re-run the program)
	outputRemoveIndividual ($ped, $ped_index, $fam_mendel_error, $ind_mendel_error, $ind_missing, $num_analyzed_marker);
	
	print STDERR "NOTICE: Program finished!\n";
	close (STDOUT); close (STDERR);	
}

#write individual remove files to help re-execute the program
sub outputRemoveIndividual {
	my ($ped, $ped_index, $fam_mendel_error, $ind_mendel_error, $ind_missing, $num_analyzed_marker) = @_;
	my %bad_ind;
	if ($output) {
		open (IND, ">$output.remove") or confess "Error: cannot write individuals that should be removed from analysis to the file $output.remove: $!";
		print STDERR "NOTICE: Writting individuals who fail to meet the --mind, --fme_threshold or --ime_threshold criteria to $output.remove file\n";
	} else {
		*IND = *STDERR;
		print STDERR "NOTICE: Identifying individuals who fail to meet the --mind --fme_threshold or --ime_threshold criteria\n";
	}
	for my $i (0 .. @$ped_index-1) {
		my ($fam, $ind) = ($ped_index->[$i][0], $ped_index->[$i][1]);
		$ped->{$fam}{$ind} or next;					#this individual is excluded from analysis (not in %ped), possibly due to --remove argument
		if ($ind_missing->[$i] and $ind_missing->[$i]/$num_analyzed_marker > $mind_threshold or $ind_mendel_error->[$i] and $ind_mendel_error->[$i]/$num_analyzed_marker > $ime_threshold) {
			$bad_ind{$i} and next;			#this individual is already written to file
			$bad_ind{$i}++;
			print IND $ped_index->[$i][0], "\t", $ped_index->[$i][1], "\n";
		}
		if ($fam_mendel_error->{$fam} and $fam_mendel_error->{$fam} / $num_analyzed_marker > $fme_threshold) {
			$bad_ind{$i} and next;			#this individual is already written to file
			$bad_ind{$i}++;
			print IND $ped_index->[$i][0], "\t", $ped_index->[$i][1], "\n";
		}
	}
	scalar (keys %bad_ind) and print STDERR "NOTICE: A total of ", scalar (keys %bad_ind), " individuals that fail to meet the --mind, --fme_threshold or --ime_threshold criteria\n";
	scalar (keys %bad_ind) and print STDERR "NOTICE: Please consider to re-run the program using the --remove argument to remove these individuals\n";

	if ($output) {
		scalar (keys %bad_ind) and print STDERR "\tCommand to execute: $orig_command_line --remove $output.remove\n";
		close (IND);
	}
}


sub analyzeHeaderLine {
	my ($header) = @_;
	$header =~ s/\s*[\r\n]+$//;
	my @header = split (/\t/, $header);
	my ($name_index, $chr_index, $pos_index, @ac_index, @gt_index, @allindex);
	
	if ($canonical_gtfile_format) {
		($name_index, $chr_index, $pos_index) = (0, 1, 2);
		if ($alleleAB) {
			@gt_index = (3..@header-1);
		} else {
			@ac_index = (3..@header-1);
		}
	} elsif ($illumina_gtfile_format) {
		for my $i (0 .. @header-1) {
			lc $header[$i] eq 'name' and $name_index = $i;
			lc $header[$i] eq 'chr' and $chr_index = $i;
			lc $header[$i] eq 'position' || lc $header[$i] eq 'pos' and $pos_index = $i;
			$header[$i] =~ m/\.Allele Calls$/i and push @ac_index, $i;
			$header[$i] =~ m/\.GType/i and push @gt_index, $i;
		}
	} elsif ($affymetrix_gtfile_format) {		#this format is created by apt-probeset-genotype, with header lines starting with #, and with genotype denoted as 0,1,2 and -1.
		1;
	} elsif ($plink_tpedfile_format) {
		@header = split (/\s+/, $header);	#separated by space rather than tab
		($name_index, $chr_index, $pos_index) = (1, 0, 3);
		if ($alleleAB) {
			@gt_index = (4..@header-1);
		} else {
			@ac_index = (4..@header-1);
		}
	} elsif ($simple_format) {
		($name_index, $chr_index, $pos_index) = (0, undef, undef);
		if ($alleleAB) {
			@gt_index = (1..@header-1);
		} else {
			@ac_index = (1..@header-1);
		}
	}
	return ($name_index, $chr_index, $pos_index, \@ac_index, \@gt_index, scalar (@header));
}

#this subroutine generate the sample index for all permutation cycles
#during each permutation cycle, the cells in the @perm_index array tells what two groups to compare (for case-control study), or what families to flip transmission/untransmission status (for family-based study), or how to scramble quantitative traits (for QT study)
sub calculatePermIndex {
	my ($cycle, $nf_index) = @_;
	my @perm_index;
	if ($cc) {
		my @allindex = (@{$nf_index->[0]}, @{$nf_index->[1]});
		for my $current_cycle (1 .. $cycle) {
			srand ($current_cycle + ($seed-1)*$cycle);		#this ganrantee that when doing 100 cycles in 10 machines, each of 1000 permutation get a differnt seed
			fisher_yates_shuffle (\@allindex);
			my @case_index = @allindex[0 .. @{$nf_index->[0]}-1];
			my @control_index = @allindex[@{$nf_index->[0]} .. @allindex-1];
			push @perm_index, [\@case_index, \@control_index];
		}
	} elsif ($tdt or $tsp) {
		for my $current_cycle (1 .. $cycle) {
			srand ($current_cycle + ($seed-1)*$cycle);
			my @switch;
			for (1.. @$nf_index) {
				push @switch, (rand () > 0.5) ? 1: 0;
			}
			push @perm_index, \@switch;
		}
	} elsif ($qt) {
		my @allindex = (0 .. @{$nf_index->[0]}-1);			#this is not a random shuffle of first element of nf_index, instead it shuffle the position (so that qt is equally shuffled!)
		for my $current_cycle (1 .. $cycle) {
			srand ($current_cycle + ($seed-1)*$cycle);
			fisher_yates_shuffle (\@allindex);			#the allindex contains shuffled number ranging from 0 to arrayelement(nfindex)-1, but the nf_index position index in the ped_index array
			push @perm_index, [@allindex];
		}
	}
	return \@perm_index;
}

sub generateNFIndex {
	my ($ped) = @_;
	my ($nf_index);
	if ($tdt) {
		$nf_index = $force_tdt ? GWA::identifyTrioFromPed ($ped, 'affected') : GWA::identifyIndepentTrioFromPed ($ped, 'affected');
		@$nf_index = sort {$a->[0] <=> $b->[0]} @$nf_index;			#sort by index of father in pedheaderfile
		@$nf_index or confess "Error: unable to find suitable trios from ped_header file $pedheaderfile";
	} elsif ($tsp) {
		$nf_index = GWA::identifyNuclearFamilyFromPed ($ped, 'affected');
		@$nf_index = sort {$a->[0] <=> $b->[0]} @$nf_index;			#sort by index of father in pedheaderfile
		@$nf_index or confess "Error: unable to find suitable trios or quartets from ped_header file $pedheaderfile";
	} elsif ($cc) {
		$nf_index = GWA::identifyCaseControlFromPed ($ped);			#first element: case; second element: control
		scalar (@{$nf_index->[0]}) or confess "Error: unable to find suitable cases from ped_header file $pedheaderfile";
		scalar (@{$nf_index->[1]}) or confess "Error: unable to find suitable controls from ped_header file $pedheaderfile";
	} elsif ($qt) {
		$nf_index = GWA::identifyQTIndFromPed ($ped);				#first element: posindex; second element: qtvalue
		scalar (@{$nf_index->[0]}) or confess "Error: unable to find suitable individuals for quantitative trait analysis from ped_header file $pedheaderfile";
	}
	return ($nf_index);
}
sub calculateAssociation {
	my ($gtfile, $ped, $ped_index, $ped_stat, $extract_snp, $exclude_snp, $snpprop) = @_;
	my ($name_index, $chr_index, $pos_index, $ac_index, $gt_index, $num_record_per_line);
	my ($nf_index, $perm_index);
	my (@dm_invalidgt, @dm_invalidchr, @dm_gt2allele, @dm_geno, @dm_maf, @dm_mme, @dm_hwe);
	my ($num_analyzed_marker) = (0);
	my (%fam_mendel_error, @ind_mendel_error, @ind_missing);
	my ($counter_line) = (0);
	my ($name_length);				#the length of the marker name in the association output file (sometimes markers such as CNVs have extraordinarily long names)
	
	#read the header line of the genotype file to make sure the number of samples in the file is the same as the ped file
	#determine whether the 'Allele Calls' or the 'GType' records in the gtfile will be used for association analysis
	open (GT, $gtfile) or confess "Error: cannot read from GT file $gtfile: $!";
	
	defined ($_ = <GT>) or confess "Error: cannot read anything from GT file $gtfile";
	chomp;
	$counter_line++;
	($name_index, $chr_index, $pos_index, $ac_index, $gt_index, $num_record_per_line) = analyzeHeaderLine ($_);
	$num_record_per_line >= 2 or confess "Error: the GT file must contain at least 2 space-delimited records per line, but the header line is: $_";
	
	if ($plink_tpedfile_format) {
		open (GT, $gtfile);		#open GT file again to start reading from first line
		$counter_line = 0;
	}
	
	if ($simple_format) {
		defined $name_index or confess "Error: unable to find columns for SNP name in the header line ($_) of GT file $gtfile";
	} else {
		defined $name_index and defined $chr_index and defined $pos_index or confess "Error: unable to find columns for 'Name', 'Chr' and 'Position' in the header line ($_) of GT file $gtfile";
	}
	if ($alleleAB) {
		$gt_index or confess "Error: --alleleAB argument is set but GType column is not found in header line ($_) of GT file $gtfile";
		if ($plink_tpedfile_format) {
			@$gt_index == 2*@$ped_index or confess "Error: Discordant number of samples: ${\(scalar @$ped_index)} in pedheaderfile $pedheaderfile versus ${\(scalar @$gt_index)} in GT file $gtfile";
		} else {
			@$gt_index == @$ped_index or confess "Error: Discordant number of samples: ${\(scalar @$ped_index)} in pedheaderfile $pedheaderfile versus ${\(scalar @$gt_index)} in GT file $gtfile";
		}
	} else {
		if (@$ac_index) {
			if ($plink_tpedfile_format) {
				@$ac_index == 2*@$ped_index or confess "Error: Discordant number of samples: $ped_stat->{num_ind} in pedheaderfile $pedheaderfile versus ${\(scalar @$ac_index)} in genotype file $gtfile";
			} else {
				@$ac_index == @$ped_index or confess "Error: Discordant number of samples: $ped_stat->{num_ind} in pedheaderfile $pedheaderfile versus ${\(scalar @$ac_index)} in genotype file $gtfile";
			}
		} elsif (@$gt_index) {
			confess "Error: Unable to find 'Allele Calls' records in header line of gtfile $gtfile but found 'GType' records (did you forget the --alleleAB argument?)";
		} else {
			confess "Error: Unable to find genotype records in the header line of gtfile $gtfile";
		}
	}

	$nf_index = generateNFIndex ($ped);					#nf_index has different meaning for different study designs
	$cycle and $perm_index = calculatePermIndex ($cycle, $nf_index);	#--cycle argument determines the number of permutation cycles

	#record the nuclear families to the first line of the EVIDENCE file
	if ($tdt || $tsp and $perfam) {					#--evidence is only supported for --tdt and --tsp analysis
		print EVI "Name";
		for my $index (@$nf_index) {
			print EVI "\t", $ped_index->[$index->[0]][0], ":";	#family id
			for my $i (0 .. @$index-1) {
				$i and print EVI ",";
				print EVI $ped_index->[$index->[$i]][1];	#individual id
			}
		}
		print EVI "\n";
	}
	
	#read and process each marker sequentially
	while (<GT>) {
		my (@record, @ac, @gt);					#ac: allele call (such as ACGT), gt: genotype call (AB only)
		my ($name, $chr, $pos);
		my ($marker_bad_flag);					#a flag showing whether marker failed one criteria (GENO, MAF, HWE, MME, etc);
		my ($marker_comment) = '';				#a comment to be printed at the end of each line of the marker info file

		#decide whether to process this marker based on --mstart and --mend arguments, and then read in this line
		$counter_line++;					#current line counter for GT file, which is one more than the counter for markers
		if ($marker_start) {
			if ($plink_tpedfile_format) {
				$counter_line >= $marker_start or next;
			} else {
				$counter_line-1 >= $marker_start or next;
			}
		}
		if ($marker_end) {
			if ($plink_tpedfile_format) {
				$counter_line <= $marker_end or last;
			} else {
				$counter_line-1 <= $marker_end or last;
			}
		}
		s/\s*[\r\n]+$//;					#discard trailing newline or return character



		if ($plink_tpedfile_format) {
			@record = split (/\s+/, $_);				#all records are separated by space
		} else {
			@record = split (/\t/, $_);				#all records are separated by tab
		}
		@record == $num_record_per_line or confess "Error: Discordant number of space/tab-delimited records in $gtfile: $num_record_per_line in header line but ${\(scalar @record)} in line $counter_line: <$_>";
		
		if (defined $chr_index) {
			($name, $chr, $pos) = @record[$name_index, $chr_index, $pos_index];
		} else {
			($name, $chr, $pos) = ($record[$name_index], 99, "NA");		#the chr is treated as a special "autosome"
		}
		$chr eq '23' and $chr = 'X';
		$chr eq '24' and $chr = 'Y';				#for compatibilit with PLINK-formatted TFAM/TPED files
		
		$verbose and print STDERR "NOTICE: Processing marker $name at chr$chr position $pos\n";

		#read genotypes for this marker on all individuals
		%$extract_snp and $extract_snp->{$name} || next;	#this SNP is not specifed in the --extract file so skipped from analysis
		$exclude_snp->{$name} and next;				#this SNP is excluded from analysis
		if ($alleleAB) {
			@ac = @record[@$gt_index];
		} else {
			@ac = @record[@$ac_index];
		}
		if ($plink_tpedfile_format) {
			my @ac1 = @ac;
			@ac = ();
			for my $i (0 .. @ac1/2-1) {
				push @ac, $ac1[$i*2].$ac1[$i*2+1];
			}
		}
		map {s/^NC/00/} @ac;
		map {s/^\-\-/00/} @ac;
		map {s/^0$/00/} @ac;					#all NO_CALL genotypes are converted in the form of 00
		map {s/ //} @ac;					#get rid of spaces in genotype calls (sometimes there is a space between the two alleles in a genotype call)
		
		#check the --geno_threshold and make sure all genotype calls are recognizable
		my @nc_geno = grep {m/^00$/} @ac;
		my @normal_geno = grep {m/^[ABCGT1234][ABCGT1234]$/} @ac;
		my $geno_missing = @nc_geno/@ac;
		my (%ac_count, %allele_count);
		if (@nc_geno+@normal_geno != @ac) {			#unrecognized genotype call (less than or more than 2 letters)
			my @invalid_geno = grep {!m/^[ABCGT01234][ABCGT01234]$/} @ac;
			print STDERR "WARNING: marker $name (chr=$chr pos=$pos) has non-recognized genotype calls (@invalid_geno)\n";
			$marker_comment .= "(invlid_genotype=" . join (",", @invalid_geno) . ")";
			@ac = (('00') x scalar (@ac));			#treat all genotypes as NULL for this marker
			$allmarker or push @dm_invalidgt, $name and $marker_bad_flag++;
		}
		$geno_missing <= $geno_threshold or push (@dm_geno, $name) and $marker_bad_flag++;	#flag this marker as bad marker (does not pass QC criteria)
		
		#convert ACGT allele calls to AB genotypes, identify unique alleles
		my $allac = join ('', @ac);		
		$ac_count{$_}++ for (@ac);
		$allele_count{$_}++ for (split (//, $allac));
		delete $allele_count{0};
		my @uniq_allele = sort keys %allele_count;
		my (%ac2gt, %gt2ac);
		if (not @uniq_allele) {					#completely missing genotype (all individuals have NoCall genotype)
			$marker_comment .= '(all_genotype_missing)';
			@gt = @ac;
		} elsif (@uniq_allele == 1) {				#only one allele found for this marker
			$marker_comment .= "(single_allele=$uniq_allele[0])";
			@gt = (('AA') x scalar (@ac));
		} elsif (@uniq_allele > 2) {
			print STDERR "WARNING: marker $name (chr=$chr pos=$pos) has more than two alleles (@uniq_allele)\n";
			$marker_comment .= "(more_than_two_allele=" . join (",", @uniq_allele) . ")";
			$allmarker or push @dm_gt2allele, $name and $marker_bad_flag++;
			@gt = (('00') x scalar (@ac));			#consider all marker genotypes to be NULL when more than 2 allele is present
		} else {						#the normal scenario
			if ($alleleAB) {
				@gt = @ac;
				map {s/^BA/AB/} @gt;			#make sure A is always listed before B in the genotype calls
			} else {
				$uniq_allele[1] eq 'B' and print STDERR "WARNING: the 'B' allele is found in genotype calls for marker $name (did you forget the --alleleAB argument?)\n";
				@ac2gt{'00', $uniq_allele[0].$uniq_allele[0], $uniq_allele[0].$uniq_allele[1], $uniq_allele[1].$uniq_allele[0], $uniq_allele[1].$uniq_allele[1]} = ('00', 'AA', 'AB', 'AB', 'BB');
				%gt2ac = reverse %ac2gt;
				@gt = map {$ac2gt{$_}} @ac;
			}
		}
		
		#calculate minor allele frequency (for founders only)
		my (%founder_allele, $maf_founder);
		if ($chr eq 'X') {
			for my $indgt (@gt[@{$ped_stat->{male_founder_index}}]) {
				$indgt eq '00' and next;
				$indgt =~ m/^(\w)\1$/ or next;		#invalid genotype calls for male chrX
				$founder_allele{$1}++;
			}
			for my $indgt (@gt[@{$ped_stat->{female_founder_index}}]) {
				$indgt eq '00' and next;
				$founder_allele{$_}++ for (split (//, $indgt));
			}
		} elsif ($chr eq 'Y') {
			for my $indgt (@gt[@{$ped_stat->{male_founder_index}}]) {
				$indgt eq '00' and next;
				$indgt =~ m/^(\w)\1$/ or next;		#invalid genotype calls for male chrY
				$founder_allele{$1}++;
			}
		} else {
			for my $indgt (@gt[@{$ped_stat->{founder_index}}]) {
				$indgt eq '00' and next;
				$founder_allele{$_}++ for (split (//, $indgt));
			}
		}
		if (keys %founder_allele == 0) {
			$maf_founder = 0;
		} elsif (keys %founder_allele == 1) {
			$maf_founder = 0;
		} else {
			$maf_founder = $founder_allele{A}/($founder_allele{A}+$founder_allele{B});
			$maf_founder > 0.5 and $maf_founder = 1-$maf_founder;
		}
		$maf_founder >= $maf_threshold or push (@dm_maf, $name) and $marker_bad_flag++;	#flag this marker as bad marker (does not pass QC criteria)

		#print the header line of the association results (the length of the header line depends on the marker name length)
		if (not defined $name_length) {
			$name_length = length ($name) + 5;			#the idea is that: given the first marker, add its length by 5, it should suffice for all subsequent markers (if not, the full marker name is still printed with mis-alignment)
			$name_length < 12 and $name_length = 12;
			outputAssocHeaderline ($name_length);
		}

		#check the Hardy-Weinberg equilibrium threshold
		my (%founder_gt, $hwe_p_value);
		if ($chr eq 'X') {
			my @founder_gt = @gt[@{$ped_stat->{female_founder_index}}];
			$founder_gt{$_}++ for (@founder_gt);
			if (grep {!m/00/} @founder_gt) {
				$hwe_p_value = kc::SNPHWE ($founder_gt{'AB'}||0, $founder_gt{'AA'}||0, $founder_gt{'BB'}||0);
			} else {
				$hwe_p_value = 'NA';				#no genotype information in founders so cannot calculate HWE
			}
		} elsif ($chr =~ m/^\d+$/) {
			my @founder_gt = @gt[@{$ped_stat->{founder_index}}];
			$founder_gt{$_}++ for (@founder_gt);
			if (grep {!m/00/} @founder_gt) {
				$hwe_p_value = kc::SNPHWE ($founder_gt{'AB'}||0, $founder_gt{'AA'}||0, $founder_gt{'BB'}||0);
			} else {
				$hwe_p_value = 'NA';				#no genotype information in founders
			}
		} else {
			$hwe_p_value = 'NA';					#this statement is reserved for Y chromosome or Mitochondria (non diploid chromosome)
			$marker_comment .= "(non_autosome_or_X)";
			$allmarker or push @dm_invalidchr, $name and $marker_bad_flag++;
		}
		if ($hwe_p_value ne 'NA' and $hwe_p_value < $hwe_threshold) {
			push @dm_hwe, $name;
			$marker_bad_flag++;
		}
		
		#calculate association test statistics
		my ($marker_mendel_error, %fam_mendel_error_recorded) = ('NA');
		if ($tdt or $tsp) {
			my ($chi2, $chi2_p, $t_count, $u_count, $index_mendel_error, $tu_trio, $permstring1, $permstring2);
			my ($hx11, $hx12, $hx22, $hy11, $hy12, $hy22, $sx1122, $sx2211, $sy12, $sy21, $hx_star);
			if ($tdt) {
				($chi2, $chi2_p, $t_count, $u_count, $index_mendel_error, $tu_trio) = GWA::calTDT (\@gt, $nf_index, $ped_index, ($chr eq 'X'), $perfam);
				$cycle and ($permstring1, $permstring2) = GWA::calTDT_perm (\@gt, $nf_index, $ped_index, ($chr eq 'X'), 0, $perm_index, $cycle);			
			
			} elsif ($tsp) {
				($chi2, $chi2_p, $sx1122, $sx2211, $sy12, $sy21, $index_mendel_error, $tu_trio) = GWA::calTSP (\@gt, $nf_index, $ped_index, ($chr eq 'X'), $perfam);
				$cycle and ($permstring1, $permstring2) = GWA::calTSP_perm (\@gt, $nf_index, $ped_index, ($chr eq 'X'), 0, $perm_index, $cycle);			
			}
			
			#check the marker mendel error threshold (the $index_mendel_error contains index for offsprings with mendelian inconsistency)
			$marker_mendel_error = @$index_mendel_error / @$nf_index;
			$marker_mendel_error <= $mme_threshold or push @dm_mme, $name and $marker_bad_flag++;
			
			#mark the family and individual as having mendel error, if this marker is still considered as a good marker
			if (not $marker_bad_flag) {
				for my $index (@$index_mendel_error) {
					$ind_mendel_error[$index]++;				#mark the child has mendelian error (but do not mark the parent as having mendelian error)
					my $fam = $ped_index->[$index][0];
					$fam_mendel_error_recorded{$fam} and next;		#this family is already processed for this marker (but two different ind has mendel error)
					$fam_mendel_error{$fam}++;				#mark the family as having one additional mendel error
					$fam_mendel_error_recorded{$fam}++;			#mark the family as having mendel error and already recorded (several individuals in same family may have mendel error)
				}
			}

			#prepare the output formatting
			if (not $marker_bad_flag) {
				if ($tdt) {
					outputTDTResult ($name_length, $name, $chr, $pos, $chi2, $chi2_p, $t_count, $u_count, \@uniq_allele, $tu_trio, $permstring1, $permstring2, $snpprop->{$name}||'');
				} elsif ($tsp) {
					outputTSPResult ($name_length, $name, $chr, $pos, $chi2, $chi2_p, $sx1122, $sx2211, $sy12, $sy21, \@uniq_allele, $tu_trio, $permstring1, $permstring2, $snpprop->{$name}||'');
				}
			}
		} elsif ($cc) {
			my ($chi2, $chi2_p, $case_af, $control_af, $permstring1, $permstring2);
			($chi2, $chi2_p, $case_af, $control_af) = GWA::calCC (\@gt, $nf_index, ($chr eq 'X'), $cellsize, undef, $ped_stat);
			$cycle and ($permstring1, $permstring2) = GWA::calCC_perm (\@gt, $perm_index, ($chr eq 'X'), $cycle, $perm_method, $ped_stat);

			if (not $marker_bad_flag) {
				outputCCResult ($name_length, $name, $chr, $pos, $chi2, $chi2_p, $case_af, $control_af, \@uniq_allele, $permstring1, $permstring2, $snpprop->{$name}||'');
			}
		} elsif ($qt) {
			my ($a, $b, $f, $f_p, $permstring1, $permstring2);
			if ($marker_comment) {							#something is wrong with the marker
				($a, $b, $f, $f_p) = qw/NA NA NA NA/;
				$cycle and $permstring1 = ',NA'x$cycle and $permstring2 = $permstring1;
			} else {
				($a, $b, $f, $f_p) = GWA::calQT (\@gt, $nf_index, ($chr eq 'X'));
				$cycle and ($permstring1, $permstring2) = GWA::calQT_perm (\@gt, $nf_index, $perm_index, ($chr eq 'X'), $cycle);
			}
			
			if (not $marker_bad_flag) {
				outputQTResult ($name_length, $name, $chr, $pos, $a, $b, $f, $f_p, \@uniq_allele, $permstring1, $permstring2, $snpprop->{$name}||'');
			}
		}
			
		#write marker information regardless of whether marker_bad_flag is set (regardless of the exclusion criteria)
		outputMarkerInfo ($name, $chr, $pos, $geno_missing, $maf_founder, $hwe_p_value, $marker_mendel_error, $marker_comment);
		
		if (not $marker_bad_flag) {
			for my $i (0 .. @ac-1) {
				$ac[$i] eq '00' and $ind_missing[$i]++;			#record individual that has NC genotype at this marker (the marker already passed quality threshold so it is a good marker)
			}
			$num_analyzed_marker++;						#number of analyzed markers (those passing threshold)
		}
	}
	print STDERR "NOTICE: Finished association analysis on $num_analyzed_marker markers that pass inclusion criteria\n";
	outputDiscardMarker (\@dm_invalidgt, \@dm_invalidchr, \@dm_gt2allele, \@dm_geno, \@dm_maf, \@dm_hwe, \@dm_mme);
	
	return (\%fam_mendel_error, \@ind_mendel_error, \@ind_missing, $num_analyzed_marker, $nf_index);
}

#print the header line in association results
sub outputAssocHeaderline {
	my ($name_length) = @_;
	my $name = "Name" . (' ' x $name_length);
	$name = substr ($name, 0, $name_length);
	if ($tdt) {
		print "$name Chr   Position A1:A2       T:U  CHI2_TDT     CHI2_P", $cycle ? "\tCHI2_PERM\tCHI2_P_PERM" : "", $snppropfile ? "\tSNP_property" : "";
	} elsif ($tsp) {
		print "$name Chr   Position       T:U_TRIO       T:U_QURT      CHI2     CHI2_P", $cycle ? "\tCHI2_PERM\tCHI2_P_PERM" : "", $snppropfile ? "\tSNP_property" : "";
	} elsif ($cc) {
		print "$name Chr   Position  A1:A2   case_AF  cont_AF ALLELIC_chi2     A_P  TREND_chi2      T_P  GENO_chi2       G_P  DOM_chi2        D_P  REC_chi2        R_P";
		if ($cycle) {
			my %permcode = (1=>'A', 2=>'T', 3=>'G', 4=>'D', 5=>'R');
			print "\t$permcode{$perm_method}_CHI2_PERM\t$permcode{$perm_method}_CHI2_P_PERM";
		}
		$snppropfile and print "\tSNP_property";
	} elsif ($qt) {
		print "$name Chr   Position  A1:A2 INTERCEPT     SLOPE         F        F_P", $cycle ? "\tF_PERM\tF_P_PERM": "", $snppropfile ? "\tSNP_property":"";
	}
	print "\n";
}

sub outputMarkerInfo {
	my (@info) = @_;
	defined $output or return;							#when --output is not set, do not write any marker information to file
	print MINFO join ("\t", @info), "\n";
}

sub outputDiscardMarker {
	my ($dm_invalidgt, $dm_invalidchr, $dm_gt2allele, $dm_geno, $dm_maf, $dm_hwe, $dm_mme) = @_;
	@$dm_invalidgt and print STDERR "NOTICE: Total of ", scalar (@$dm_invalidgt), " markers are discarded from analysis due to having invalid (unrecognizable) genotype calls\n";
	$output and @$dm_invalidgt and print STDERR "#dm_invalidgt#\n", join ("\n", @$dm_invalidgt), "\n#dm_invalidgt#\n";
	@$dm_invalidchr and print STDERR "NOTICE: Total of ", scalar (@$dm_invalidchr), " markers are discarded from analysis due to being on non-autosome and non-X chromosomes\n";
	$output and @$dm_invalidchr and print STDERR "#dm_invalidchr#\n", join ("\n", @$dm_invalidchr), "\n#dm_invalidchr#\n";
	@$dm_gt2allele and print STDERR "NOTICE: Total of ", scalar (@$dm_gt2allele), " markers are discarded from analysis due to having more than 2 alleles in genotype calls\n";
	$output and @$dm_gt2allele and print STDERR "#dm_gt2allele#\n", join ("\n", @$dm_gt2allele), "\n#dm_gt2allele#\n";
	@$dm_geno and print STDERR "NOTICE: Total of ", scalar (@$dm_geno), " markers are discarded from analysis due to having genotype frequency lower than geno_threshold $geno_threshold\n";
	$output and @$dm_geno and print STDERR "#dm_geno#\n", join ("\n", @$dm_geno), "\n#dm_geno#\n";
	@$dm_maf and print STDERR "NOTICE: Total of ", scalar (@$dm_maf), " markers are discarded from analysis due to having minor allele frequency in founders lower than maf_threshold $maf_threshold\n";
	$output and @$dm_maf and print STDERR "#dm_maf#\n", join ("\n", @$dm_maf), "\n#dm_maf#\n";
	@$dm_hwe and print STDERR "NOTICE: Total of ", scalar (@$dm_hwe), " markers are discarded from analysis due to having Hardy-Weinberg equilibium P value lower than hwe_threshold $hwe_threshold\n";
	$output and @$dm_hwe and print STDERR "#dm_hwe#\n", join ("\n", @$dm_hwe), "\n#dm_hwe#\n";
	@$dm_mme and print STDERR "NOTICE: Total of ", scalar (@$dm_mme), " markers are discarded from analysis due to having mendelian errors more than marker_mendel_error_threshold $mme_threshold\n";
	$output and @$dm_mme and print STDERR "#dm_mme#\n", join ("\n", @$dm_mme), "\n#dm_mme#\n";
}

sub outputQTResult {
	my ($name_length, $name, $chr, $pos, $a, $b, $f, $f_p, $uniq_allele, $permstring1, $permstring2, $property) = @_;
	my ($allele1, $allele2) = @$uniq_allele;
	
	if (length ($name) < $name_length) {
		$name = substr ($name . (' 'x$name_length), 0, $name_length);
	}

	$allele1 ||= '-';
	$allele2 ||= '-';
	$chr = substr ("    $chr", -4, 4);
	$pos = substr ("           $pos", -11, 11);
	print $name, $chr, $pos, "   $allele1:$allele2 ";
	if ($a eq 'NA') {
		print '        NA', ' ', '        NA';
	} else {
		print sprintf ("%10.4g", $a), ' ', sprintf ("%10.4g", $b);		#the b value might be -3.469e-17 (so there won't be space between a and b)
	}
	print ' ';		#some people reported lack of space between SLOPE and F

	if ($f eq 'NA') {
		print '        NA', ' ', '        NA';
	} else {
		print sprintf ("%10.4g", $f), ' ', sprintf ("%10.4g", $f_p);
	}
	$permstring1 and print "\t", $permstring1;
	$permstring2 and print "\t", $permstring2;

	if ($snppropfile) {
		print $property ? "\t$property" : "\tNA";
	}
	print "\n";
}

sub outputCCResult {
	my ($name_length, $marker, $chr, $pos, $chi2, $chi2_p, $case_af, $control_af, $uniq_allele, $permstring1, $permstring2, $property) = @_;
	my ($allele1, $allele2) = @$uniq_allele;
	
	if (length ($marker) < $name_length) {
		$marker = substr ($marker . (' 'x$name_length), 0, $name_length);
	}

	$allele1 ||= '-';
	$allele2 ||= '-';	
	$chr = substr ("    $chr", -4, 4);
	$pos = substr ("           $pos", -11, 11);
	print $marker, $chr, $pos, "   $allele1:$allele2", sprintf ("%10.4g", $case_af), sprintf ("%10.4g", $control_af);
	
	for my $i (0 .. @$chi2-1) {
		if ($chi2->[$i] eq 'NA') {
			print '        NA', ' ', '        NA';
		} else {
			print sprintf ("%10.4g", $chi2->[$i]), ' ', sprintf ("%10.4g", $chi2_p->[$i]);
		}
	}
	$permstring1 and print "\t", $permstring1;
	$permstring2 and print "\t", $permstring2;
	
	if ($snppropfile) {
		print $property ? "\t$property" : "\tNA";
	}
	print "\n";
}

sub outputTDTResult {
	my ($name_length, $name, $chr, $pos, $chi2, $chi2_p, $t_count, $u_count, $uniq_allele, $tu_trio, $permstring1, $permstring2, $property) = @_;
	
	if (length ($name) < $name_length) {
		$name = substr ($name . (' 'x$name_length), 0, $name_length);
	}
	$chr = substr ("    $chr", -4, 4);
	$pos = substr ("           $pos", -11, 11);
	if ($chi2 eq 'NA') {
		($chi2, $chi2_p) = ('        NA', '        NA');
	} else {
		$chi2 = sprintf ("%10.4g", $chi2);
		$chi2_p = sprintf ("%10.4g", $chi2_p);
	}
	
	my ($allele1, $allele2) = @$uniq_allele;
	$allele2 ||= '-';
	my $transmission_status = substr ("          $t_count:$u_count", -10, 10);
	
	print $name, $chr, $pos, "   $allele1:$allele2", $transmission_status, $chi2, ' ', $chi2_p;
	$permstring1 and print "\t", $permstring1;
	$permstring2 and print "\t", $permstring2;

	if (@$tu_trio) {				#when --evidence argument is set, write the detailed output evidence to the file handle
		print EVI "$name\t", join ("\t", @$tu_trio), "\n";
	}
	
	if ($snppropfile) {
		print $property ? "\t$property" : "\tNA";
	}
	print "\n";
}

sub outputTSPResult {
	my ($name_length, $name, $chr, $pos, $chi2, $chi2_p, $sx1122, $sx2211, $sy12, $sy21, $uniq_allele, $tu_trio, $permstring1, $permstring2, $property) = @_;
	my ($trio_tu_ratio, $quartet_tu_ratio);

	if (length ($name) < $name_length) {
		$name = substr ($name . (' 'x$name_length), 0, $name_length);
	}
	$chr = substr ("    $chr", -4, 4);
	$pos = substr ("           $pos", -11, 11);
	
	if ($chi2 eq 'NA') {
		($chi2, $chi2_p) = ('        NA', '        NA');
	} else {
		$chi2 = sprintf ("%10.4g", $chi2);
		$chi2_p = sprintf ("%10.4g", $chi2_p);
	}
	
	my ($allele1, $allele2) = @$uniq_allele;
	$allele2 ||= '-';
	$trio_tu_ratio = substr ("          $allele1($sy12):$allele2($sy21)", -15, 15);
	$quartet_tu_ratio = substr ("          $allele1($sx1122):$allele2($sx2211)", -15, 15);
	
	print $name, $chr, $pos, $trio_tu_ratio, $quartet_tu_ratio, $chi2, ' ', $chi2_p;
	$permstring1 and print "\t", $permstring1;
	$permstring2 and print "\t", $permstring2;
	
	if (@$tu_trio) {				#when --evidence argument is set, write the detailed output evidence to the file handle
		print EVI "$name\t", join ("\t", @$tu_trio), "\n";
	}

	if ($snppropfile) {
		print $property ? "\t$property" : "\tNA";
	}
	print "\n";
}

sub outputSampleStat {
	my ($ped_index, $fam_mendel_error, $ind_mendel_error, $ind_missing, $num_analyzed_marker) = @_;
	if ($output) {
		for my $i (0 .. @$ped_index-1) {
			my ($famid, $indid, $fatid, $motid) = ($ped_index->[$i][0], $ped_index->[$i][1], $ped_index->[$i][2], $ped_index->[$i][3]);
			$ind_mendel_error->[$i] ||= 0;
			$ind_missing->[$i] ||= 0;
			print SINFO "$famid\t$indid\t$num_analyzed_marker\t$ind_missing->[$i]\t", $ind_missing->[$i]/$num_analyzed_marker, "\t$ind_mendel_error->[$i]\t", $ind_mendel_error->[$i]/$num_analyzed_marker, "\n";
		}
		
		for my $famid (keys %$fam_mendel_error) {
			print FINFO $famid, "\t", $fam_mendel_error->{$famid}, "\t", $num_analyzed_marker, "\t", $fam_mendel_error->{$famid} / $num_analyzed_marker, "\n";
		}
	}
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

sub readSNPPropFile {
	my ($snppropfile) = @_;
	my (%snpprop, $maxlength);
	open (SNPPROP, $snppropfile) or confess "Error: cannot read from snppropfile $snppropfile: $!";
	while (<SNPPROP>) {
		m/^(\S+)\t(.+)/ or confess "Error: invalid record found in snppropfile (at least 2 tab-delimited records expected): <$_>";
		$snpprop{$1} = $2;
		$maxlength ||= length ($2);
		$maxlength < length ($2) and $maxlength = length ($2);
	}
	return (\%snpprop, $maxlength);
}

=head1 SYNOPSIS

 calculate_association.pl [arguments] <ped-header-file> <GTfile>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation

 	Analysis Type:
 	    --tdt			use TDT analysis for independent trios
 	    --tsp			use TSP analysis for combined trio and quartet
 	    --pdt			use PDT analysis for general pedigree (not implemented yet)
 	    --cc			case-control study (default analysis option)
 	    --force_tdt			force TDT analysis for non-independent trios in multiplex families
 	    --qt			perform quantitative trait analysis (default is binary trait)

	Output Control:
 	    --output <file>		output root file name (default=STDOUT)
 	    --perfam			record detailed transmission ratio for each family
 	    --snppropfile <file>	a SNP property file to annotate association results
 	    --(no)flush			flush cache after each read/write operation (default=ON)

	Inclusion/Exclusion Criteria for Markers:
	    --allmarker			consider all markers (no exclusion on maf, hwe, geno, mme)
 	    --exclude <file(s)>		specify markers to exclude from analysis
 	    --extract <file(s)>		specify markers to extract (use) in analysis
 	    --geno_threshold <float>	genotype frequency threshold for SNP (default=0.1)
 	    --maf_threshold <float>	minor allele frequency threshold for SNP (default=0.01)
 	    --mme_threshold <float>	marker mendelian error threshold (default=0.1)
 	    --hwe_threshold <float>	Hardy-Weinberg equilibrium P-value threshold (default=0.001)
 	    --mstart <int>		index (1-based) of first marker to include analysis
 	    --mend <int>		index (1-based) of last marker to include in analysis
 	    --cellsize <int>		minimum number (inclusive) in cell in contingency table

	Inclusion/Exclusion Criteria for Samples:
	    --allind			consider all individuals (no exclusion on mind, fme, ime)
 	    --remove <file(s)>		specify individuals to remove from analysis
 	    --keep <file(s)>		specify individuals to keep in analysis
 	    --mind_threshold <float>	missing genotype threshold for individuals (default=0.1, re-run program)
 	    --fme_threshold <float>	family mendelian error threshold (default=1, OPTION NOT DEFINED YET!)
 	    --ime_threshold <float>	individual offspring mendelian error threshold (default=0.02, re-run program)

	File Format Specification
 	    --alleleAB			GT file use AB alphabet for alleles (rather than ACGT)
 	    --ab			GT file use AB alphabet for alleles (rather than ACGT)
 	    --canonical_gtfile_format	data file in canonical format (columns are name, chr, pos, genotypes)
 	    --simple_format		data file in simple format (columns are name, genotypes)
 	    --illumina_gtfile_format	data file in Illumina FullDataTable format
 	    --affymetrix_gtfile_format	data file in Affymetrix Power Tools genotype call format (to be implemented)
 	    --plink_tpedfile_format	data file in TFAM/TPED file format generated by PLINK

	Analysis-specific Options:
 	    --cycle <int>		specify cycles for phenotype permutation
 	    --seed <int>		specify randomization seed for permutation
 	    --perm_method <int>		case-control permutation method (1=allelic as default, 2=trend, 3=genotypic, 4=dom, 5=rec)
 	    --covfile <file>		specify external covariant file for trend test (-cc) or quantitative trait (-qt) (to be implemented)
 	    --covcolumn <int>		specify columns in the ped_header file that contain covariants (to be implemented)

 Function: perform genome-wide association analysis, given a ped_header file and 
 a GT file that can be in various input data format.
 
 Example: calculate_association.pl ped_header gtfile -out result -ab
          calculate_association.pl -plink input.tfam input.tped -out result
          calculate_association.pl ped_header gtfile -out result -cycle 25 -seed 1 -perm 2
          calculate_association.pl ped_header gtfile -out result -remove file1,file2,file3 -exclude file4 -exclude file5

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--tdt>

perform transmission/disequilibrium test on independent trios to test for 
association

=item B<--tsp>

perform the TSP test on trios and quartets to test for association

=item B<--pdt>

perfrom pedigree disequilibrium test on general pedigrees to test for 
association

=item B<--cc>

perform case-control comparisons to test for association

=item B<--force_tdt>

force to TDT on non-independent trios to jointly test linkage and association

=item B<--qt>

perform quantitative trait analysis using simple linear regression models

=item B<--output>

specify output root file name. When this argument is not set, the main 
association results will be written to STDOUT; when this argument is set, 
multiple output files will be generated containing detailed information for each 
marker and each individual.

=item B<--perfam>

print out detailed association evidence (transmission ratio in TDT or TSP tests) 
for each family for each marker (by default, the association test statistic is a 
summary of all families)

=item B<--snppropfile>

a file that specify properties for each marker, such as gene annotation, 
alternative name and so on. The property information will be printed in the 
association results in each line for each SNP marker to help annotate and 
interpret the association results.

=item B<--(no)flush>

this argument is on by default: it forces flushing input/output pipes after each 
read/write operation, so that one can monitor program progress in real-time (for 
example, when writting results to a file, or when running program via rsh or 
ssh). It presumably increase disk overhead so one can turn off this argument by 
--noflush: this is very important to do when running the program in parallel in 
many computers that share the same central storage!

=item B<--allmarker>

use all markers for analysis without exclusion criteria

=item B<--exclude>

specify one or several files that contains markers (one per line) that should be 
excluded from analysis.

=item B<--extract>

specify one or several files that contains markers (one per line) that should be 
used in analysis.

=item B<--geno_threshold>

genotype no-call rate threshold for inclusion of markers in analysis (if the 
marker has no-call rate equal to or less than this threshold)

=item B<--maf_threshold>

minor allele frequency (MAF) threshold for inclusion of markers in analysis (if 
the marker has MAF equal to or more than this threshold)

=item B<--mme_threshold>

marker mendelian error threshold for inclusion of markers in analysis (if for 
this marker, the fraction of trios/quartets with medelian inconsistencies is 
equal to or less than this threshold)

=item B<--hwe_threshold>

Hardy-Weinberg equilibrium threshodl for inclusion of markers in analysis (if 
the marker has HWE P-values more than or equal to this threshold)

=item B<--mstart>

specify the first marker to analyze. Together with the --mend argument, these 
can be used to analyze only a subset of markers in the entire GT file, and is 
suitable to split the genome-wide analysis job into multiple parallele 
computers.

=item B<--mend>

specify the last marker to analyze. Together with the --mstart argument, these 
can be used to analyze only a subset of markers in the entire GT file, and is 
suitable to split the genome-wide analysis job into multiple parallele 
computers.

=item B<--cellsize>

minimum number of counts in each cell in the contingency table for chi2-based 
association test, including genotypic association test, dominant model 
association test, recessive model association test.

=item B<--allind>

use all individuals for analysis without exclusion criteria

=item B<--remove>

specify one or several files containing family and individual identifiers (one 
per line) to be removed from analysis

=item B<--keep>

specify one or several files containing family and individual identifiers (one 
per line) to be used in analysis

=item B<--mind>

specifying the missingness threshold (if the fraction of missing genotypes for 
all markers exceed this threshold, this individual is not uesd in analysis)

=item B<--fme_threshold>

family mendelian error threshold for inclusion of families in association 
analysis (if fraction of markers with Mendelian inconsistency in this family 
exceed this threshold, this family is not used in analysis). THIS ARGUMENT IS 
SET AS 1, AS IT IS NOT WELL-DEFINED YET. IF ONE OF TWO OFFSPRINGS IN A FAMILY 
HAVE EXCESSIVE MENDELIAN INCONSISTENCY, DO WE STILL INCLUDE THE OTHER GENUINE 
OFFSPRING IN ANALYSIS? OR DO WE COMPLETELY DISCARD THIS FAMILY? FME_THRESHOLD 
DOES NOT SEEM TO BE REALLY USEFUL HERE, COMPARED TO IME_THRESHOLD.

=item B<--ime_threshold>

individual mendelian error threshold for inclusion of individuals in association 
analysis (if fraction of markers with Mendelian inconsistency in this individual 
exceed this threshold, this individual is not used in analysis)

=item B<--alleleAB>

specify that the genotype calls in GT file are in AB alphabet, rather than ACGT alphabet.

=item B<--ab>

specify that the genotype calls in GT file are in AB alphabet, rather than ACGT alphabet.

=item B<--canonical_gtfile_format>

specify that the GT file is in canonical gtfile format, where in each line the 
first three tab-delimited columns are marker name, chromosome and position, and 
the rest columns are genotype calls for all individuals.

=item B<--simple_format>

specify that the GT file is in simple format, where the first column of each 
line represent markre name and the rest represent genotype calls for all 
individuals.

=item B<--illumina_gtfile_format>

specify that the GT file is in Illumina FullDataTable format, where the first 
line of the file indicate the column index for Name, Chr, Position and GType for 
all individuals.

=item B<--affymetrix_gtfile_format>

(NOT IMPLEMENTED YET) specify that the GT file is in Affymetrix data format, 
which is the output file from the Affymetrix Power Tools software. It will have 
multiple header lines starting with #, followed by the genotype calls. An 
optional --confidence file will be needed so that only highly confident SNP 
calls are used in analysis.

=item B<--plink_tpedfile_format>

specify that the GT file is in TFAM/TPED file format generated by the PLINK 
software. This options allows PLINK users to quickly use their existing data 
files for analysis by this program. (use the --recode --transpose option in 
PLINK to generate TFAM/TPED files).

=item B<--cycle>

specify the permutation cycle (default cycle=0, meaning no permutation is 
requested)

=item B<--seed>

specify the randomization seed (default seed=1). When running multiple 
permutations in multiple machines (for example, 100 permutations on each of 100 
machines, for a total of 1000 permutations), it is important to specify 
different seeds for differnet machines.

=item B<--perm_method>

permutation method for case-control study (1=allelic, 2=trend, 3=genotypic, 
4=dom, 5=rec): five different methods can be used to calculate association 
statistic for case-control studies, so we have to specify which method is used 
for calculation of permutated chi2 and P values. The default method is allelic 
association test.

=item B<--covfile>

specify an external file that contains covariates for each individual. The file 
contains tab-delimited records, with the first two columns in each line 
corresponding to family id and individual id, and the rest of the columns 
specifying covariates (binary trait or quantitative values). Missing values can 
be represneted as "NA" and individuals with missing values will be ignored in 
analysis.

=item B<--covcolumn>

specify that the covariates are all contained within the ped_header file itself, 
and give the column of these covariates. Multiple columns are separated by 
comma.

=back

=head1 DESCRIPTION

This program is used to perform genome-wide association (GWA) analysis to 
identify the possible association between genetic markers (such as single 
nucleotide polymorphisms and copy number variations) and one or more phenotypic 
measures (such as disease status or quantitative expression values). To learn 
more about GWA, the following papers are good general references:

1. Risch N, Merikangas K (1996) The future of genetic studies of complex human 
diseases. Science 273:1516-1517

2. Hirschhorn J, Daly M (2006), Genome-wide association studies for common diseases 
and complex traits, Nature Reviews Genetics 6, 95-108

3. Wang WYS, Barratt1 BJ, Clayton DG  &  Todd JA (2006) Genome-wide association 
studies: theoretical and practical concerns, Nature Reviews Genetics 6, 109-118

4. McCarthy MI, Abecasis GR, Cardon LR, Goldstein DB, Little J, Ioannidis JP, 
Hirschhorn JN (2008) Genome-wide association studies for complex traits: 
consensus, uncertainty and challenges, Nature Reviews Genetics 9, 356-369

Given a list of individuals, their pedigree relationships, and their 
phenotypes/covariates in a file called ped_header file, as well as genotypes for 
many genetic markers for each individual, this program will test the association 
of each marker with the primary phenotype (such as disease status) adjusting for 
other phenotypes (such as age and sex). Three kinds of testing strategies can be 
analyzed by this program: case-control study (with five different 
association/regression models) and family-based study with 
transmission/disequilibrium test, as well as trio/quartet based study with Tsp 
test (a preliminary form of the Pedigree Disequilibrium Test by Martin et al), 
for which the details were given in I<Tests for linkage and association in 
nuclear families> Am J Hum Genet. 1997 Aug;61(2):439-48 (PUBMED ID: 9311750).

The following sections will give more detailed descriptions on the program:


=head2 Section 1. =======================File formats==========================

=head3 B<ped_header file>

The ped_header file can be considered a superset of the normal PED file used in 
most genetic analysis. It is a tab-delimited text file that records information in one line 
for each individual, and the first six columns of the line are family id, 
individual id, father id, mother id, sex and affection status (or quantitative 
trait). The following columns can be other phenotypes or covariates or genotypes 
or any information (such as race, name, address, etc) one wants to include.

An example of the ped_header file is shown below:

	AU0215  AU021503        AU021502        AU021501        1       2	White	AA
	AU0215  AU021502        0       0       1       1	Mixed	BB
	AU0215  AU021501        0       0       2       1	White	AB
	AU0215  AU021504        AU021502        AU021501        1       1	Black	AA
	AU0285  AU028504        AU028502        AU028501        2       2	Unknown	BB
	AU0285  AU028505        AU028502        AU028501        2       1	Unknown	AB
	AU0285  AU028501        0       0       2       1	White	AA
	AU0285  AU028503        AU028502        AU028501        2       2	White	BB
	AU0285  AU028502        0       0       1       1	White	AB
	AU0121  AU012104        AU012102        AU012101        2       2	White	AB

In addition to the first 6 columns seen in a typical PED file, there is two 
extra columns: one specifying race, and the other specifying the genotype of a 
SNP marker. These two extra columns can be covariates in association analysis. 
Many more additional columns can be added in the file following the first 6 
columns; therefore, it is even possible to use a regular PED file containing all 
SNP genotypes as a ped_header file!

Note that for binary trait, missing phenotypes can be represented as 0; but for 
quantitative trait, it should be represented as B<NA>.

=head3 B<GT file>

The GT file is a genotype file that contains information for one marker per 
line. In a canonical GT file, the first three columns for the file are markerid, 
chromosome and position, while the following columns contain the genotypes for 
each individual. The format for the GT file can be slightly relaxed by the --
illumina_gtfile_format and --affymetrix_gtfile_format arguments. In the case of 
Illumina GT file, one can simply export selected columns in the FullDataTable 
from the BeadStudio software as a GT file, without conforming to the strict 
criteria for a canonical GT file; instead, the definition of columns can be read 
from the first line of the GTfile.

To explain this in more detail, see the first 10 lines of a canonical GT file 
below:

	Name    Chr     Position        sample1        sample2        sample3        sample4        sample5        sample6        sample7
	rs10013542      4       60194938        BB      BB      BB      BB      BB      BB      BB
	rs10013547      4       175988420       BB      BB      BB      BB      BB      BB      BB
	rs10013571      4       28821638        AA      AA      AA      AA      AB      AA      AB
	rs10013576      4       27228651        AB      BB      AB      AB      BB      BB      BB
	rs10013588      4       132717336       AA      AA      AA      AA      AA      AA      AA
	rs10013604      4       43295922        AA      AB      AB      AB      AA      AA      AA
	rs1001362       16      55232359        BB      BB      AB      AB      AB      BB      BB
	rs10013632      4       117373068       AA      AA      AA      AA      AA      AA      AA
	rs10013649      4       162694792       AA      AA      AA      AA      AA      AA      AA
	rs10013727      4       142935271       AB      AB      BB      AB      BB      BB      BB

Each line contains 10 tab-delimited columns, and the first three columns are 
general description of the SNP marker, while the rest columns are the genotypes of 
the markers.

The first 10 lines of a illumina-GTfile is shown below:

	Index	Name    Chr     Position        4028211240_B.GType      4028211240_B.B Allele Freq      4028211240_B.Log R Ratio        4030178116_B.GType      4030178116_B.B Allele Freq      4030178116_B.Log R Ratio
	1	rs389518        2       141691876       AB      0.4122693       -0.3346987      AB      0.5449692       -0.05919398
	2	rs7761056       6       169420041       AB      0.5415538       -0.1324195      AB      0.5364565       -0.08725725
	3	rs2081188       19      21180835        AB      0.5462196       -0.1279078      AB      0.4796815       -0.2688029
	4	rs2951747       12      11518610        BB      1       -0.02597118     BB      0.9963943       0.02994001
	5	rs648102        12      127952571       AA      0.002531445     0.2534947       BB      0.9680079       0.04287543
	6	rs659113        3       42331840        NC      0.2538367       1.197143        AA      0.03036925      -0.2546444
	7	rs2929374       3       19923057        AB      0.5194996       0.1479239       BB      1       -0.08358008
	8	rs2619298       3       103346166       AA      0.004948121     -0.008488706    BB      0.8150933       -0.1196488
	9	rs2959523       18      17567106        BB      1       -0.06788081     BB      0.9987407       0.07631278
	10	rs3104240       18      24521689        AA      0.01723832      -0.2327651      AA      0.01594706      -0.3815838

As you can see, the column position I<Name>, I<Chr> and I<Position> columns are 
determined by the first line (so-called header line of GT file), and the 
genotype for each marker is also determined by the ***.GType column in the 
header line. In other words, the Illumina-GTfile can contain much more 
information than a canonical GTfile, as long as the header line contains enough 
informatoin for identifying the required columns.

Missing genotypes can be specified as B<NC>, B<-->, B<00> or B<0>.

=head3 B<TFAM/TPED file>

To facilitate PLINK users to use this program, the user can alternatively 
generate the TAFM/TPED file in PLINK and feed them into the program.

For example, suppose you have three PLINK files in binary format as input.bim, 
input.fam and input.bed. First it is necessary to transform the files to 
TAFM/TPED text-based format in PLINK:

	plink --bfile input --out output --transpose --recode

Two files including output.tfam and output.tped will be generated.

You can then use the --plink_tpedfile_format argument in this program to read in 
the output.tfam and output.tped file.

=head3 B<Simple-format GT file>

The simple file format is similar to the GT file above, with the exception that 
there is no Chr and Position column in the file. To use simple-format GT file, 
the user needs to specify the --simple_format argument to the program.

This file can be easily generated from Illumina BeadStudio report file format, 
or Affymetrix BRLMM/BirdSeed output file format, or Affymetrix Chiamo fs file 
format, with some very simple one-line Perl programming.

=head3 B<output files>

The main output file will be the association results. When --output argument is 
not set, the association results will be printed to STDOUT (standard output, 
typically terminal screen). Otherwise, the association results will be written 
to the file specified by --output argument.

Many I<NOTICE> and I<WARNING> message will be printed during the program run, 
which gives intermediate notification of the progress of the program, as well as 
possible minor problems encountered. when the --output argument is set, these 
messages will be also written to a log file $output.log.

An additional marker information file will be generated. This file contains 
detailed information for each marker, including allele frequency, Hardy-Weinberg 
equilibrium test and so on. This is a text file too with one marker per line, 
but with many more columns than the association results file, and contains 
information for all markers (including markers not passing MAF, HWE and MME 
thresholds).

For example, the first 10 lines of a minfo file is given below:

	Name    Chr     Position        genotype_missing        minor_allele_frequency  hwe_p_value     marker_mendel_error     comments
	rs3094315       99      NA      0.00284981476204047     0.157470551766894       1.91321160785967e-08    NA
	rs6672353       99      NA      0.00256483328583642     0       1       NA      (single_allele=G)
	rs4040617       99      NA      0.00797948133371331     0.12753036437247        0.0260114518331774      NA
	rs2980300       99      NA      0.00740951838130522     0.156386292834891       1.2885943744048e-06     NA
	rs2905036       99      NA      0.000569962952408094    0.000309885342423303    1       NA
	rs4245756       99      NA      0.00085494442861214     0.000309789343246547    1       NA
	rs4075116       99      NA      0.00740951838130522     0.267857142857143       0.52945617513906        NA
	rs9442385       99      NA      0.00170988885722428     0.0721153846153846      0.290169504144135       NA
	rs10907175      99      NA      0.0042747221430607      0.0930521091811415      0.91712160159049        NA

The first 10 lines of a sinfo file is given below:

	family_id       individual_id   num_analyzed_marker     nocall_count    nocall_rate     mendelian_error_count   mendelian_error_rate
	0       WTCCC66751      391943  643     0.00164054467103635     0       0
	0       WTCCC66761      391943  1076    0.00274529714779955     0       0
	0       WTCCC66771      391943  1614    0.00411794572169933     0       0
	0       WTCCC66723      391943  5272    0.0134509354676573      0       0
	0       WTCCC66731      391943  1577    0.0040235442398512      0       0
	0       WTCCC66742      391943  748     0.00190844076817292     0       0
	0       WTCCC66752      391943  706     0.00180128232931829     0       0
	0       WTCCC66762      391943  1145    0.0029213431544893      0       0
	0       WTCCC66772      391943  1256    0.00320454760003368     0       0

When family information is available, the finfo file will contain the following information:

	family_id       num_analyzed_marker     mendelian_error_count   mendelian_error_rate

When some samples are recommended to be removed (due to call rate, Mendelian 
error, etc), they will be written to a remove file.


=head2 Section 2. ===================Testing strategy==========================

=over 8

=item * B<Case-control association test>

Case-control association test compared two phenotype groups with a binary label, 
such as disease status. There are typically five differnet models to test the 
association of binary phenotypes with bialleleic genotypes, which are briefly 
described below:

1. Genotypic association test (2df association test)

The genotypic association test, also referred to as 2df (2 degrees of freedom) 
associaiton test, uses 2X3 contingency table to test the association between 
binary phenotype label with a trinary genotype label (AA, AB and BB). This is a 
standard chi2 test; however, caution should be exercised where cells in the 
contingency tables are small (less than 5): P-values in this case might be 
inflated or deflated. Typically, when sample size is sufficiently large (>2000) 
and when alleles are relatively common (MAF threshold is set to >1%), this is 
rarely the case. The --cell_size argument can be used to specify the minimum 
cell: the marker will be excluded from association test if any one of the six 
cells in the contingency table is less than 5.

In rare cases, one type of genotypes might be completely absent from the 2X3 
contingency table in both case group and control group. In this case, the 
genotypic association test uses a 2X2 table and has one degree of freedom now.

2. Allelic association test

The allelic association test compares the allele frequency between case groups 
and control groups, using a chi2 test.

3. Cochran-Armitage trend test

The Cochran-Armitage trend association test essentially tries to reduce the 2df 
in genotypic association test to 1df, by assuming an additive mode of 
inheritance. It is conceptually identical to a logistic regression, where 
dependent variables are disease status, and independent variables are genotypes. 
For more information see Armitage, P (1955). Tests for linear trends in 
proportions and frequencies. Biometrics, 11:375-386.

4. Dominant model association test

The dominant model treat AB and BB as the same label (B allele has dominant 
effect over A), and then test the association of binary phenotype labels with 
binary genotype labels. Note that some other software (such as PLINK) may use 
minor/major allele to specify a dominant model, but here we always assume that B 
has dominant effect over A. This test is affected by the --cellsize argument: 
when cells in contingeny table is smaller than the --cellsize argument, the 
association chi2 and P values will be reported as I<NA>.

5. Recessive model association test

The recessive model treat AA and AB as the same label (B allele has recessive 
effect), and then test the association of binary phenotype labels with binary 
genotype label. Note that some other software (such as PLINK) may use 
minor/major allele to specify a recessive model, but here we always assume that 
B has recessive effect over A. This test is affected by the --cellsize argument: 
when cells in contingeny table is smaller than the --cellsize argument, the 
association chi2 and P values will be reported as I<NA>.

=item * B<Transmission/disequilibrium test>

The transmission/disequilibrium test examine the over-transmission of a 
particular allele from heterozygous parents to their offsprings. Strictly 
speaking, for complex diseases, as we always assume that disease phenotype is 
linked with every gene, we are using TDT to test for association only. 
Therefore, the TDT test can be only applied to father-mother-offspring trios in 
GWA studies. (However, some groups do not agree with this point, and use TDT on 
quartet and extended families as well to test for association.)

To know more about the TDT test, see the original publication Spielman RS, 
McGinnis RE, Ewens WJ.Transmission test for linkage disequilibrium: the insulin 
gene region and insulin-dependent diabetes mellitus (IDDM). AJHG 52:506-16 
(1993).

To read more about using TDT in genome-wide association studies, see Ewens WJ, 
Spielman RS. What is the significance of a significant TDT? Hum Hered. 60:206-10 
(2005). Note that some groups apparently do not agree on the main point of this 
paper; and that PLINK does not perform any adjustment for multiplex families. I 
myself have performed simulation study and demonstrated that in the case of 
multiplex families, TDT gives inflated type I error, and results in loss-of- 
power (this is not published and will not be published, it is for my own fun). 
If you have quartets, the TSP (see below) but not TDT, is a valid test, and 
shows correct type I error and improved power. For general pedigrees, there are 
several other programs for association tests, including FBAT, LAMP and so on, 
but they are generally slow and do not scale up to genome-wide association 
studies well.

=item * B<TSP test>

The TSP test is named after Martin et al. I basically used a highly similar 
model as they have used, with one exception (slight modification). This 
exception potentially decrease the power of the method, but I can never 
understand their original pseudo-counting method in their paper. I think it is 
just fine to use equal pseudocounts for A and B alleles when there is ambiguity 
in transmission status, in that the power loss is generally negligible.

=item * B<quantitative trait analysis>

This program also implements a simple regression model for quantitative trait 
analysis. In the regression model, the dependent variable is the QT, while the 
independent variable is the genotypes (AA=0, AB=1 and BB=2). The association 
test is testing the null hypothesis that the slope of the regression is zero. 
Note that in this regression model, we are assuming an additive model for B 
alleles. If one wants to use a multiplicative model, one has to manually convert 
the QT values to logarithm scale in the ped_header file, before running this 
program.

=back

=head2 Section 3. ===========================Quality control===================

Genotype quality control (QC) is an important step for genetic association 
analysis, and may help eliminate spurous association or improve power to detect 
true associations.

Two types of exclusion/inclusion thresholds are implemented in this program: one 
types of QC operates on per-marker level, and the other type of QC operates on 
per-individual level. Unlike other programs (such as PLINK) that does both QC 
steps semi-simultaneously, this program treat these two types of QC differently.

The program scans each line in the GT file sequentially to test for associatoin 
for each marker; therefore, the per-marker level QC is performed during the 
scanning step: markers that fail to meet the --geno_threshold, --maf_threshold, 
--mme_threshold and --hwe_threshold are not printed to association results 
(however, these information for all markers are written to the minfo file so 
that you can check the exact values for these excluded markers).

After the initial scanning step, we would also get a count of the --
ime_threshold, --fme_threshold, --mind and can test whether some indiviudals 
failed to meet these threshold. The individuals that did not meet the per-
individual QC criteria are written to a file as *.indremove.

Therefore, the per-individual level QC is performed after scanning the entire GT 
file. Since all the marker statistics were calculated using all the individuals, 
one must use the --remove argument to add this *.indremove file to the command 
line and re-run the program again to get the correct association results. (Note 
that you can specify the --remove argument in command line multiple times, and 
information from all the --remove arguments will be combined together for 
removing individuals)

The whole process sounds a little weird, but it does make a lot of sense to do 
this in a two-step way. In summary, we use the per-marker scanning to identify 
individuals that fail to meet QC criteria (the criteria is calculated using only 
a set of markers that pass the per-marker QC criteria), and then remove these 
individuals to calculate the final association results.

Several additional notes are described below:
 
=over 8

=item * B<minor allele fequency (MAF)>

The minor allele frequency is calculated on founders only! Founders refers to 
individuals without father and mother annotations in the ped_header file. Note 
that if an individual has father and mother annotation in the ped_header file, 
but such father and mother do not exist in the ped_header file, this individual 
is not treated as founder. (this usually happens when some indivduals in family-
based study do not have DNA samples or are not genotyped.) A warning message 
will be printed out to STDERR during the pedigree structure inspection step.

For male chrX and chrY, if the genotype is heterozygotes, it is treated as 
NoCall instead; otherwise, it contribute one count to the total allele count (as 
opposed to two for female chrX).

Female chrY is ignored in MAF calculation.

For case-control studies, the allele frequency is calculated by pooling cases 
and control together, and stored in the *.minfo file when the --output argument 
is set. However, in the association results, the case allele frequency and 
control allele frequency are also printed.

=item * B<Hardy-Weinberg equilibrium (HWE)>

The HWE calculation uses exact test. See Wigginton JE, Cutler DJ, Abecasis GR, A 
Note on Exact Tests of Hardy-Weinberg Equilibrium. AJHG, 76: 887-93 (2005) for 
more details.

The HWE is calculated only for autosomes, and for female chrX. Other chromosome 
will be assigned a HWE value of NA.

Similar to MAF, the HWE is calculated on founders only! For case-conrol studies, 
again the HWE is calculated using combined samples. This may cause some 
problems, since one may argue that cases should violate HWE in markers that are 
associated with disease status. In the future, an argument will be added that 
change this default behavior and only calculate HWE on controls to eliminate 
such concerns.

=back

=head2 Section 4. ================Advanced usage of the program================

Although the program was originally written for simple association analysis, 
over time it has evolved into a testing platform for implementing various 
analytical strategies beyond simple chi2 based association test. Below I briefly 
describe several advanced usage of the program that are still in experimental 
phase.

=over 8

=item *B<using remove, keep, exclude and extract files>

The program can accept the same types of remove, keep, exclude and extract files 
as used by the PLINK software, so that one can conveniently specify a subset of 
samples or a subset of markers to be used in association analysis. In addition, 
it is possible to specify multiple files separated by comma in the command line, 
which is something that is not available in the PLINK software; this behavior of 
using multiple filter files has some advantage in practical settings.

=item * B<parallel execution of association test>

The program was mainly implemented in Perl: speed-wise, its performance is not 
satisfactory. Although considerable amounts of bits of the program will be 
gradually transferred to C subroutines, sometimes it is still desirable to 
parallelize the program to achieve fast speed. Generally speaking, it takes only 
a few hours to perform association test on 3000 subjects on 550K markers. 
However, if testing association on 10,000 subjects, the time would be much 
longer and it may be desirable to split the job into multiple subset of jobs and 
run them in different computers in parallel.

At this moment, an easy-to-implement strategy for parallelization is through the 
--mstart and --mend argument, which specify the index of the first marker and 
last marker to include in analysis (the index starts from 1 and end at the last 
marker included in GT file).

For example, to analyze 500K markers in 10 machines/CPUs, one can execute the 
following commands in different machines separately:

	calculate_association.pl pedheaderfile GTfile -cc -out part1 -mstart 1 -mend 50000
	calculate_association.pl pedheaderfile GTfile -cc -out part2 -mstart 50001 -mend 100001
	...
	calculate_association.pl pedheaderfile GTfile -cc -out par10 -mstart 450001 -mend 500000

All the output files (excluding LOG file) can be concatenated together easily, 
since they are all line-based file. The above steps can achieve considerable 
speed increase.

Certainly, one could also generate 10 marker files each containing 50K markers 
and then use the --extract argument for each marker file in 10 different 
machines. One could also line-split the GT file to 10 part, but add one extra 
header line to all files except the first one, then run the program on all of 
them. There are many ways for easy parallelization, but the --mstart and --mend 
seem to be the easiest way.

=item * B<parallel permutation procedure>

The permutation procedure takes much longer than than simple association test. 
Therefore, it is generally a good idea to split the total cycle of requested 
permutations (such as 1000) to different computer nodes (for example, use --
cycle 100 for each of 10 computational nodes). Since the default randomization 
seed is set to 1 in the program, do make sure to use a different seed in each 
different machine!

	calculate_association.pl pedheaderfile GTfile -cc -out part1 -seed 1 -cycle 100 -noflush
	calculate_association.pl pedheaderfile GTfile -cc -out part2 -seed 2 -cycle 100 -noflush
	...
	calculate_association.pl pedheaderfile GTfile -cc -out par10 -seed 10 -cycle 100 -noflush

In the above command, we specify a different randomizatoin seed for each 
command, and each command perform 100 permutation cycles. Note that the 
--noflush argument is added to the command line, which ensures that input/output 
are NOT flushed during each read/write operation. This is especially important 
when doing parallel computation, to decrase the system overhead on hard disks.

=item * B<Adjusting for co-variates>

Adjusting for co-variates is NOT IMPLEMENTED due to lack of stimulus/incentive. 
Basically for case-control studies, this implicates that logistic regression 
should be used to incorporate co-variate information. For quantitative trait 
analysis, this implicates that standard multiple linear regression should 
suffice.

=back
                                                                                                                                                                                                                                                                                                                                                                       
