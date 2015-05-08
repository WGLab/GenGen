package GWA;
use strict;
use warnings;
use Carp;
use Exporter qw( import );

eval {
	require "kc.pm";
};
if ($@) {
	require Config;
	my $arch = $Config::Config{archname} || 'UNKNOWN';
	print STDERR "GenGen compilation error: Your system architecture is '$arch', which is not compatible with pre-compiled executables.\n";
	print STDERR "GenGen compilation error: Please download source code from www.openbioinformatics.org/gengen/ and compile executable program.\n";
	exit (100);
}

our $VERSION = 			'$Revision: 310 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-02-06 16:35:51 -0800 (Sat, 06 Feb 2010) $';

our (@EXPORT, @EXPORT_OK, %EXPORT_TAGS);

my @tdt_tags = qw/identifyTrioFromPed identifyIndepentTrioFromPed analyzeTrio/;
my @tsp_tags = qw/identifyNuclearFamilyFromPed/;
my $verbose;

@EXPORT = ();				# symbols to export by default
@EXPORT_OK = ();			# symbols to export on request
%EXPORT_TAGS = (tdt => [@tdt_tags], tsp => [@tsp_tags]);

our %tdt_trans = (	ABAAAA=>1, ABBBAB=>1, ABABAB=>1, AAABAA=>1, BBABAB=>1, ABABAA=>2,
			ABAAAB=>0, ABBBBB=>0, ABABBB=>0, AAABAB=>0, BBABBB=>0,
			AAAAAA=>0, BBBBBB=>0, AABBAB=>0, BBAAAB=>0,
			
			AYABAA=>1, AYABAY=>1, AYABBY=>0, BYABAB=>1, BYABAY=>1, BYABBY=>0,
			AYABAB=>0, BYABBB=>0,
			AYAAAA=>0, AYAAAY=>0, BYBBBB=>0, BYBBBY=>0, AYBBAB=>0, AYBBBY=>0, BYAAAB=>0, BYAAAY=>0,
		);

#tdt_untrans is equal to transmission of alternative allele B
our %tdt_untrans =(	ABAAAB=>1, ABBBBB=>1, ABABAB=>1, AAABAB=>1, BBABBB=>1, ABABBB=>2,
			ABAAAA=>0, ABBBAB=>0, ABABAA=>0, AAABAA=>0, BBABAB=>0,
			AAAAAA=>0, BBBBBB=>0, AABBAB=>0, BBAAAB=>0,
		
			AYABAB=>1, AYABAY=>0, AYABBY=>1, BYABBB=>1, BYABAY=>0, BYABBY=>1,
			AYABAA=>0, BYABAB=>0,
			AYAAAA=>0, AYAAAY=>0, BYBBBB=>0, BYBBBY=>0, AYBBAB=>0, AYBBBY=>0, BYAAAB=>0, BYAAAY=>0,	
		);

1;



#when disease flag is set, only trios that contains affected child is in output
sub identifyTrioFromPed {
	my ($ped, $disease) = @_;
	my %fam_trio;
	my @trio_index;
	my $num_trio;					#total number of trios for analysis
	for my $famid (keys %$ped) {
		for my $indid (keys %{$ped->{$famid}}) {
			$disease and $ped->{$famid}{$indid}[3] == 2 || next;		#when disease is set, the offspring must be affected individual
			my ($fatherid, $motherid) = ($ped->{$famid}{$indid}[0], $ped->{$famid}{$indid}[1]);
			if ($fatherid and $motherid and $ped->{$famid}{$fatherid} and $ped->{$famid}{$motherid}) {
				push @{$fam_trio{$famid}}, [$fatherid, $motherid, $indid];
				push @trio_index, [$ped->{$famid}{$fatherid}[4], $ped->{$famid}{$motherid}[4], $ped->{$famid}{$indid}[4]];
				$num_trio++;
			}
		}
	}
	print STDERR "NOTICE: Identifying trios for TDT analysis: $num_trio trios from ${\(scalar keys %fam_trio)} families\n";
	return (\@trio_index);
}

sub identifyQTIndFromPed {
	my ($ped) = @_;
	my (@posindex, @qtvalue);
	for my $famid (keys %$ped) {
		for my $indid (keys %{$ped->{$famid}}) {
			$ped->{$famid}{$indid}[3] eq 'NA' and next;		#QT is unknown
			push @posindex, $ped->{$famid}{$indid}[4];
			push @qtvalue, $ped->{$famid}{$indid}[3];
		}
	}
	my @nf_index = (\@posindex, \@qtvalue);
	return (\@nf_index);							#nf_index contains not only posindex, but also qt values (this is different from cc or tdt, where phenotype can be grouped into separate arrays)
}

sub identifyCaseControlFromPed {
	my ($ped) = @_;
	my (@case_index, @control_index, @cc_index);
	for my $famid (keys %$ped) {
		for my $indid (keys %{$ped->{$famid}}) {
			if ($ped->{$famid}{$indid}[3] == 2) {
				push @case_index, $ped->{$famid}{$indid}[4];
			} elsif ($ped->{$famid}{$indid}[3] == 1) {
				push @control_index, $ped->{$famid}{$indid}[4];
			}
		}
	}
	@case_index = sort {$a<=>$b} @case_index;
	@control_index = sort {$a<=>$b} @control_index;
	@cc_index = (\@case_index, \@control_index);
	return (\@cc_index);
}

#this is the correct subroutine to use for identifying trios in TDT analysis (if trio are not independent, the statistic is wrong)
#when disease flag is set, only trios that contains affected offspring is kept in output
sub identifyIndepentTrioFromPed {
	my ($ped, $disease) = @_;
	my %fam_trio;
	my @trio_index;
	my $num_trio;					#total number of trios for analysis
	my %trio_found;
	for my $famid (keys %$ped) {
		for my $indid (keys %{$ped->{$famid}}) {
			$disease and $ped->{$famid}{$indid}[3] == 2 || next;	#when disease is set, the individual must be affected offspring in the trio
			my ($fatherid, $motherid) = ($ped->{$famid}{$indid}[0], $ped->{$famid}{$indid}[1]);
			$trio_found{$famid, $fatherid, $motherid} and next;	#this father and mother already contribute to trio, so skip them
			if ($fatherid and $motherid and $ped->{$famid}{$fatherid} and $ped->{$famid}{$motherid}) {
				push @{$fam_trio{$famid}}, [$fatherid, $motherid, $indid];
				push @trio_index, [$ped->{$famid}{$fatherid}[4], $ped->{$famid}{$motherid}[4], $ped->{$famid}{$indid}[4]];
				$num_trio++;
				$trio_found{$famid, $fatherid, $motherid}++;	#mark this father and mother already contribute to trio
			}
		}
	}
	print STDERR "NOTICE: Identifying trios for TDT analysis: $num_trio trios from ${\(scalar keys %fam_trio)} families\n";
	return (\@trio_index);
}

sub identifyNuclearFamilyFromPed {
	my ($ped, $disease) = @_;
	my (%fam_par_ind);
	my @nf_index;
	for my $famid (keys %$ped) {
		for my $indid (keys %{$ped->{$famid}}) {
			$disease and $ped->{$famid}{$indid}[3] == 2 || next;	#when disease is set, the individual must be affected offspring
			my ($fatherid, $motherid) = ($ped->{$famid}{$indid}->[0], $ped->{$famid}{$indid}->[1]);
			if ($fatherid and $motherid and $ped->{$famid}{$fatherid} and $ped->{$famid}{$motherid}) {
				push @{$fam_par_ind{"$famid,$fatherid,$motherid"}}, $indid;
			}
		}
	}
	
	my ($num_trio, $num_quartet) = (0, 0);
	for my $fam_par (keys %fam_par_ind) {
		my ($famid, $fatherid, $motherid) = split (/,/, $fam_par);
		my @indid = @{$fam_par_ind{$fam_par}};
		my $father_index = $ped->{$famid}{$fatherid}[4];
		my $mother_index = $ped->{$famid}{$motherid}[4];
		my @ind_index = map {$ped->{$famid}{$_}[4]} @indid;
		push @nf_index, [$father_index, $mother_index, @ind_index]; #my @temp=map {$_+1} ($father_index, $mother_index, @ind_index); print "@temp\n";
		@ind_index == 1 and $num_trio++;
		@ind_index >= 2 and $num_quartet++;
	}
	print STDERR "NOTICE: Identifying ${\(scalar @nf_index)} nuclear families (including $num_trio trios and $num_quartet quartets)\n";
	return (\@nf_index);
}

sub calQT {
	my ($gt, $qt_index, $chrx_flag) = @_;
	my ($posindex, $qtvalue) = @$qt_index;		#posindex is the position in the ped_header (for example, posindex=<3 1 0 9> qtindex=<1.2 3.1 1.3 1.4>)
	my @validelement;				#validelement is the valid array element that should be selected from 0, 1, 2, 3 only

	my (@x, @y);
	for my $i (0 .. @$posindex-1) {
		if ($gt->[$posindex->[$i]] ne '00' and $qtvalue->[$i] ne 'NA') {	#genotype information is available, phenotype is available (this should not really happen, since the identifyQTIndFromPed () takes care of the NA)
			push @x, $gt->[$posindex->[$i]];
			push @y, $qtvalue->[$i];
		}
	}


	my ($a, $b, $f, $p) = qw/0 0 0 0/;
	map {s/AA/0/ or s/AB/1/ or s/BB/2/} @x;
	my $ndata = scalar (@x);
	
	if (not $ndata) {
		return (qw/NA NA NA NA/);				#all genotype are no call
	}
	
	#for my $i (0 .. @x-1) {
	#	print "$x[$i]\t$y[$i]\n";
	#}
	
	kc::reg_linear (\@x, \@y, $ndata, \$a, \$b, \$f, \$p);
	return ($a, $b, $f, $p);
}

sub calQT_perm {
	my ($gt, $qt_index, $perm_index, $chrx_flag, $cycle) = @_;
	my ($posindex, $qtvalue) = @$qt_index;
	my ($string1, $string2);
	for my $current_cycle (0 .. $cycle-1) {
		my $newindex = $perm_index->[$current_cycle];
		my @newqtvalue = @$qtvalue[@$newindex];
		my ($v1, $v2, $v3, $v4) = calQT ($gt, [$posindex, \@newqtvalue], $chrx_flag);
		if ($v3 ne 'NA') {
			$v3 = sprintf ("%.3g", $v3);
			$v4 = sprintf ("%.3g", $v4);
		}
		$string1 .= ",$v3";
		$string2 .= ",$v4";
	}
	return ($string1, $string2);
}

#calculate case-control statistics by permutation, the perm_flag is used to specify which of the five models to use for permutated chi2 and P values
sub calCC_perm {
	my ($gt, $perm_index, $chrx_flag, $cycle, $perm_flag, $ped_stat) = @_;
	my ($chi2, $chi2_p, $v1, $v2, $string1, $string2);
	$perm_flag ||= 1;						#default is to use allelic association test for permutation
	for my $current_cycle (0 .. $cycle-1) {
		($chi2, $chi2_p) = calCC ($gt, $perm_index->[$current_cycle], $chrx_flag, 0, $perm_flag, $ped_stat);
		($v1, $v2) = ($chi2->[0], $chi2_p->[0]);
		if ($v2 ne 'NA') {
			$v1 = sprintf ("%.3g", $v1);			#reformat the numeric values for easier representation as string
			$v2 = sprintf ("%.3g", $v2);			#reformat the numeric values for easier representation as string
		}
		$string1 .= ",$v1";
		$string2 .= ",$v2";
	}
	return ($string1, $string2);
}

sub calCC {
	my ($gt, $cc_index, $chrx_flag, $cellsize, $perm_flag, $ped_stat) = @_;
	my (%case_gt_count, %control_gt_count);
	my $case_index = $cc_index->[0];
	my $control_index = $cc_index->[1];
	my ($chi2, $chi2_p) = qw/0 0/;
	my (@chi2, @chi2_p);
	my ($case_af, $control_af);				#A allele frequency in cases and controls
	my @table;						#contingency table holding genotype counts
	
	@case_gt_count{'AA', 'AB', 'BB'} = (0, 0, 0);
	@control_gt_count{'AA', 'AB', 'BB'} = (0, 0, 0);
	
	#follow the same rule as used in PLINK: http://pngu.mgh.harvard.edu/~purcell/plink/faq.shtml#faq9
	#in PLINK: For the --model test and Hardy-Weinberg calculations, male X chromosome genotypes are excluded. 
	
	if ($chrx_flag) {
		my %fi = map {$_, 1} @{$ped_stat->{female_founder_index}};		#female index (only female are used in individual-based association; later on some changes are applied for allele-based association
		for my $ci (@$case_index) {
			if ($fi{$ci}) {
				$case_gt_count{$gt->[$ci]}++;
			}
		}
		for my $ci (@$control_index) {
			if ($fi{$ci}) {
				$control_gt_count{$gt->[$ci]}++;
			}
		}
	} else {
		$case_gt_count{$_}++ for (@$gt[@$case_index]);
		$control_gt_count{$_}++ for (@$gt[@$control_index]);
	}
	@table = ($case_gt_count{AA}, $case_gt_count{AB}, $case_gt_count{BB}, $control_gt_count{AA}, $control_gt_count{AB}, $control_gt_count{BB});


	#the @table_allele array is useful later on to calculate allele frequency for chrX
	my @table_allele = (2*$table[0]+$table[1], 2*$table[2]+$table[1], 2*$table[3]+$table[4], 2*$table[5]+$table[4]);

	#step 1: calculate allelic association test (chi2-based); SPECIAL CHROMOSOME X HANDLING HERE, SINCE THIS IS AN ALLELE-BASED ASSOCIATION TEST
	if (not $perm_flag or $perm_flag == 1) {
		
		if ($chrx_flag) {
			my %mi = map {$_, 1} @{$ped_stat->{male_founder_index}};
			for my $ci (@$case_index) {
				if ($mi{$ci}) {
					if ($gt->[$ci] eq 'AA') {
						$table_allele[0] += 2;		#per PLINK convention, males are treated as females, as if they have two alleles, in allelic association test
					} elsif ($gt->[$ci] eq 'BB') {
						$table_allele[1] += 2;
					}
				}
			}
			for my $ci (@$control_index) {
				if ($mi{$ci}) {
					if ($gt->[$ci] eq 'AA') {
						$table_allele[2] += 2;
					} elsif ($gt->[$ci] eq 'BB') {
						$table_allele[3] += 2;
					}
				}
			}
		}
		
		kc::chi2test_2by2table (\@table_allele, \$chi2, \$chi2_p);
		if ($chi2 < 0) {					#this means something is going wrong
			$chi2 = 'NA';
			$chi2_p = 'NA';
		}
		push @chi2, $chi2;
		push @chi2_p, $chi2_p;
		($chi2, $chi2_p) = qw/0 0/;
	}
	
	if ($perm_flag and $perm_flag == 1) {
		return (\@chi2, \@chi2_p);
	}
	
	#step 2: calculate Cochran-Armitage trend association test (chi2-based)
	#this test essentially tries to reduce the 2df in genotypic association test to 1df, by assuming an additive mode of inheritance. see Armitage, P (1955). Tests for linear trends in proportions and frequencies. Biometrics, 11:375-386
	if (not $perm_flag or $perm_flag == 2) {
		kc::chi2test_trend_2by3table (\@table, \$chi2, \$chi2_p);
		if ($chi2 < 0) {					#this means something is going wrong
			$chi2 = 'NA';
			$chi2_p = 'NA';
		}
		push @chi2, $chi2;
		push @chi2_p, $chi2_p;
		($chi2, $chi2_p) = qw/0 0/;
	}
	
	if ($perm_flag and $perm_flag == 2) {
		return (\@chi2, \@chi2_p);
	}

	#step 3. calculate genotypic association test (2df chi2-based association test) by 2X3 contingency table	
	if (not $perm_flag or $perm_flag == 3) {
		if ($cellsize and scalar (grep {$_<$cellsize} @table)) {
			$chi2 = 'NA';
			$chi2_p = 'NA';
		} else {
			kc::chi2test_2by3table (\@table, \$chi2, \$chi2_p);
			if ($chi2 < 0) {				#this means something is going wrong
				$chi2 = 'NA';
				$chi2_p = 'NA';
			}
		}
		push @chi2, $chi2;
		push @chi2_p, $chi2_p;
		($chi2, $chi2_p) = qw/0 0/;
	}
	
	if ($perm_flag and $perm_flag == 3) {
		return (\@chi2, \@chi2_p);
	}
	
	#step 4: calculate dominant model association test (B has dominant effect over A, regardless of minor or major allele!)
	if (not $perm_flag or $perm_flag == 4) {
		my @table_dom = ($table[0], $table[1]+$table[2], $table[3], $table[4]+$table[5]);
		if ($cellsize and scalar (grep {$_<$cellsize} @table_dom)) {
			$chi2 = 'NA';
			$chi2_p = 'NA';
		} else {
			kc::chi2test_2by2table (\@table_dom, \$chi2, \$chi2_p);
			if ($chi2 < 0) {					#this means something is going wrong
				$chi2 = 'NA';
				$chi2_p = 'NA';
			}
		}
		push @chi2, $chi2;
		push @chi2_p, $chi2_p;
		($chi2, $chi2_p) = qw/0 0/;
	}
	
	if ($perm_flag and $perm_flag == 4) {
		return (\@chi2, \@chi2_p);
	}
	
	#step 5: calculate recessive model association test (B has recessive effect over A, regardless of minor or major allele!)
	if (not $perm_flag or $perm_flag == 5) {
		my @table_rec = ($table[0]+$table[1], $table[2], $table[3]+$table[4], $table[5]);
		if ($cellsize and scalar (grep {$_<$cellsize} @table_rec)) {
			$chi2 = 'NA';
			$chi2_p = 'NA';
		} else {
			kc::chi2test_2by2table (\@table_rec, \$chi2, \$chi2_p);
			if ($chi2 < 0) {					#this means something is going wrong
				$chi2 = 'NA';
				$chi2_p = 'NA';
			}
		}
		push @chi2, $chi2;
		push @chi2_p, $chi2_p;
	}
	
	if ($perm_flag and $perm_flag == 5) {
		return (\@chi2, \@chi2_p);
	}

=head1	
	#calculate case and control allele frequency
	if ($table[0]+$table[1]+$table[2]) {
		$case_af = ($table[0]+$table[1]/2)/($table[0]+$table[1]+$table[2]);
	} else {
		$case_af = 'NA';
	}
	if ($table[3]+$table[4]+$table[5]) {
		$control_af = ($table[3]+$table[4]/2)/($table[3]+$table[4]+$table[5]);
	} else {
		$control_af = 'NA';
	}
=cut

	if ($table_allele[0]+$table_allele[1]) {
		$case_af = $table_allele[0]/($table_allele[0]+$table_allele[1]);
	} else {
		$case_af = 'NA';
	}
	if ($table_allele[2]+$table_allele[3]) {
		$control_af = $table_allele[2]/($table_allele[2]+$table_allele[3]);
	} else {
		$control_af = 'NA';
	}
	
	return (\@chi2, \@chi2_p, $case_af, $control_af);
}

#this subroutine calculates the TDT statistic, given a list of genotypes from many individuals and a list of trios
#$gt: an array containing genotypes for multiple individuals (some of them form trios), where AB alphabet for genotype is used and 0 means NoCall genotype
#$trio_index: an array containing indexes (positions within the gt array) for multiple trios
#$chrx_flag: indicate whether the GT is for chrX
#$tu_flag: indicate whether detailed tu_ratio for each trio is stored and returned
#$perm_trio_index: specify an array containing the flipping status for each trio in $trio_index
#$current_cycle: specify the current permutation cycle to get the permutation indicator from perm_index
sub calTDT {
	my ($gt, $trio_index, $ped_index, $chrx_flag, $tu_flag, $perm_trio_index) = @_;
	my ($t_count, $u_count, $chi2, $chi2_p) = qw/0 0 NA NA/;
	my (@index_mendel_error, @tu_trio);
	my ($index, $trio_gt, $fgt, $mgt, $ogt);					#father, mother, offspring genotypes
	
	for my $i (0 .. @$trio_index-1) {
		$index = $trio_index->[$i];
		($fgt, $mgt, $ogt) = @$gt[@$index];
		if (grep {m/^0/} ($fgt, $mgt, $ogt)) {					#NoCall genotype is specified as 00
			$tu_flag and push @tu_trio, 'NC';
			next;
		}
		
		if ($chrx_flag) {							#father is heterozygous (possibly wrong genotype call and skipped)
			if ($fgt eq 'AB') {
				$tu_flag and push @tu_trio, 'NC';
				next;
			}
			$fgt =~ s/.$/Y/;						#father: the other allele is Y allele
			if ($ped_index->[$index->[2]][2] == 1) {			#if offspring sex is male
				if ($ogt eq 'AB') {					#male is heterozygous (possibly wrong genotype call and skipped)
					$tu_flag and push @tu_trio, 'NC';
					next;
				}
				$ogt =~ s/.$/Y/;
			}
		}

		$trio_gt = $fgt . $mgt . $ogt;
		my ($t, $u) = ($tdt_trans{$trio_gt}, $tdt_untrans{$trio_gt});		#make two variables here, representing the transmission/untransmission status
		if (not defined $t) {
			push @index_mendel_error, $index->[2];
			$tu_flag and push @tu_trio, 'ME';
			next;
		}
		
		if ($perm_trio_index) {
			$perm_trio_index->[$i] and ($t, $u) = ($u, $t);		#permutate transmission/untransmission status
		}
		
		$t_count += $t;							#total transmission count
		$u_count += $u;							#total untransmission count
		$tu_flag and push @tu_trio, $t - $u;
	}

	if ($t_count+$u_count) {							#if there is count at all (some markers are heterozygous in parents, and has no mendelian inconsistency)
		$chi2 = ($t_count-$u_count)**2 / ($t_count+$u_count);
		$chi2_p = 1 - kc::cdf_chi2 (1, $chi2);
	}
	return ($chi2, $chi2_p, $t_count, $u_count, \@index_mendel_error, \@tu_trio);
}

sub calTDT_perm {
	my ($gt, $trio_index, $ped_index, $chrx_flag, $tu_flag, $perm_index, $total_cycle) = @_;
	my ($chi2_string, $chi2_p_string);						#a string of chi2 and chi2_P values for all the permutation cycles
	for my $current_cycle (0 .. $total_cycle-1) {
		my ($chi2, $chi2_p) = calTDT ($gt, $trio_index, $ped_index, $chrx_flag, 0, $perm_index->[$current_cycle]);
		if ($chi2 ne 'NA') {
			$chi2 = sprintf ("%.3g", $chi2);					#reformat the numeric values for easier representation as string
			$chi2_p = sprintf ("%.3g", $chi2_p);					#reformat the numeric values for easier representation as string
		}
		$chi2_string .= ",$chi2";
		$chi2_p_string .= ",$chi2_p";
	}
	return ($chi2_string, $chi2_p_string);
}
	

#calculate TDT statistic for quartets using Martin et al (1997, AJHG, Tests for linkage and association in nuclear families) formula for combinations of data from affected sib pairs and singletons
#hx12: parents with marker genotype 12 and with two affected children
#hy12: parrents with marker genotype 12 and with one affected child
#hx: number of heterozygous parents with two affected children
#hx_star: number of heterozygous parents to transmit same allele to both children
#sx1122: number of parents with genotype 12 who have two affected children and give allele 1 to both children
#sy12: number of parents with genotype 12 who have a sinle child and give allele 1 to that child
sub calTSP {
	my ($gt, $nf_index, $ped_index, $chrx_flag, $tu_flag, $perm_nf_index) = @_;
	my ($hx11, $hx12, $hx22, $hy11, $hy12, $hy22) = qw/0 0 0 0 0/;
	my ($hx_star, $hx, $sx1122, $sx2211, $sy12, $sy21) = qw/0 0 0 0 0 0/;
	my (@index_mendel_error, @tu_trio, $chi2, $chi2_p);
	my ($index, $trio_gt, $fgt, $mgt, $o1gt, $o2gt);

	for my $i (0 .. @$nf_index-1) {
		$index = $nf_index->[$i];
		($fgt, $mgt, $o1gt, $o2gt) = @$gt[@$index];			#TSP test only consider the trios and quartets and delete other additional siblings
		if (grep {m/^0/} ($fgt, $mgt, $o1gt, $o2gt||'')) {		#missing genotype for the marker for this nuclear family
			$tu_flag and push @tu_trio, 'NC';
			next;
		}
		
		if ($chrx_flag) {
			if ($fgt eq 'AB') {					#genotype error: father cannot be heterozygotes
				$tu_flag and push @tu_trio, 'NC';
				next;
			}
			$fgt =~ s/.$/Y/;
			if ($ped_index->[$index->[2]][2] == 1) {		#first male offspring cannot be heterozygotes
				if ($o1gt eq 'AB') {
					$tu_flag and push @tu_trio, 'NC';
					next;
				}
				$o1gt =~ s/.$/Y/;
			}
			if ($o2gt and $ped_index->[$index->[3]][2] == 1) {	#second offspring cannot be heterozygotes
				if ($o2gt eq 'AB') {
					$tu_flag and push @tu_trio, 'NC';
					next;
				}
				$o2gt =~ s/.$/Y/;
			}
		}
		
		if (not defined $o2gt) {					#trio (the o2 genotype is not defined)
			$trio_gt = $fgt . $mgt . $o1gt;
			my ($t, $u) = ($tdt_trans{$trio_gt}, $tdt_untrans{$trio_gt});
			if (not defined $t) {
				push @index_mendel_error, $index->[2];
				$tu_flag and push @tu_trio, 'ME';
				next;
			}

			$perm_nf_index->[$i] and ($t, $u) = ($u, $t);		#permutate transmission/untransmission status

			$sy12 += $t;
			$sy21 += $u;
			$hy12 += ($t+$u);					#number of het parents with one child
			$tu_flag and push @tu_trio, $t-$u;
		} else {
			my ($t1, $t2) = ($tdt_trans{$fgt . $mgt . $o1gt}, $tdt_trans{$fgt . $mgt . $o2gt});
			my ($sx1122_add, $sx2211_add) = qw/0 0/;		#rather than recording sx1122 and sx2211 directly, we record their potential additions (so that we can permutate the additions)
			if (not defined $t1) {
				push @index_mendel_error, $index->[2];
			}
			if (not defined $t2) {
				push @index_mendel_error, $index->[3];
			}
			if (not defined $t1 or not defined $t2) {		#both markers should be clear of mendel error to proceed
				$tu_flag and push @tu_trio, 'ME';
				next;
			}
			
			my $t12 = $t1 + $t2;
			if ($fgt eq 'AB' and $mgt eq 'AB') {			#both father and mother are heterozygous
				$hx12 += 2;
				if ($t12 == 4) {				#ABABAAAA
					$sx1122_add = 2;			#$sx1122 += 2;
				} elsif ($t12 == 3) {				#ABABAAAB ABABABAA
					$sx1122_add = 1;			#$sx1122++;
				} elsif ($t12 == 2) {				#ABABABAB
					$hx_star++;				#either transmit 0 or 2, so take their mean
				} elsif ($t12 == 1) {				#ABABBBAB ABABABBB
					$sx2211_add = 1;			#$sx2211++;
				} elsif ($t12 == 0) {				#ABABBBBB
					$sx2211_add = 2;			#$sx2211 += 2;
				}
				$tu_flag and push @tu_trio, $t12 * 2 - 4;
			} elsif ($fgt eq 'AB' or $mgt eq 'AB') {		#either father or mother are heterozygous but not both
				$hx12++;
				if ($t12 == 2) {				#AAABAAAA ABAAAAAA ABBBABAB BBABABAB
					$sx1122_add = 1;			#$sx1122++;
				} elsif ($t12 == 0) {				#AAABABAB ABAAABAB ABBBBBBB BBABBBBB
					$sx2211_add = 1;			#$sx2211++;
				}
				$tu_flag and push @tu_trio, $t12 * 2 - 2;
			} else {
				$tu_flag and push @tu_trio, 0;
			}
			
			$perm_nf_index->[$i] and ($sx1122_add, $sx2211_add) = ($sx2211_add, $sx1122_add);	#permutate entire quartet transmission status
			$sx1122 += $sx1122_add;
			$sx2211 += $sx2211_add;					#add desired counts (after permutation) to the sx1122 and sx2211
		}
	}
	$hx_star += ($sx1122+$sx2211);

	if ($hx12+$hy12==0 or $hx_star+$hy12==0) {
		($chi2, $chi2_p) = ("NA", "NA");
	} else {
		$chi2 = (2 * $hx12 + $hy12) / (4 * $hx_star + $hy12) * (2 * $sx1122 + $sy12 - 2 * $sx2211 - $sy21)**2 / (2 * $hx12 + $hy12);
		$chi2_p = 1 - kc::cdf_chi2 (1, $chi2);
	}@tu_trio and print STDERR "ERROR!!! @tu_trio\n";
	return ($chi2, $chi2_p, $sx1122, $sx2211, $sy12, $sy21, \@index_mendel_error, \@tu_trio);
}

sub calTSP_perm {
	my ($gt, $nf_index, $ped_index, $chrx_flag, $tu_flag, $perm_index, $total_cycle) = @_;
	my ($chi2_string, $chi2_p_string);						#a string of chi2 and chi2_P values for all the permutation cycles
	for my $current_cycle (0 .. $total_cycle-1) {
		my ($chi2, $chi2_p) = calTSP ($gt, $nf_index, $ped_index, $chrx_flag, 0, $perm_index->[$current_cycle]);
		if ($chi2 ne 'NA') {
			$chi2 = sprintf ("%.3g", $chi2);					#reformat the numeric values for easier representation as string
			$chi2_p = sprintf ("%.3g", $chi2_p);					#reformat the numeric values for easier representation as string
		}
		$chi2_string .= ",$chi2";
		$chi2_p_string .= ",$chi2_p";
	}
	return ($chi2_string, $chi2_p_string);
}

sub readPedHeaderFile {
	my ($pedheaderfile, $qt) = @_;
	my ($index, $ped, @ped_index) = (0);			#the starting index is 0, which is the first column after marker, chr, position
	open (PEDHEADER, $pedheaderfile) or confess "Error: cannot read from pedheaderfile $pedheaderfile: $!";
	while (<PEDHEADER>) {
		s/\s*[\r\n]+$//;				#discard trailing spaces and return characters
		my @record = split (/\s+/, $_);			#the six fields are family id, individual id, father id, mother id, sex, affection status (phenotype)
		@record >= 6 or confess "Error: invalid record found in pedheaderfile $pedheaderfile (at least six space-delimited fields expected): <$_>";
		my ($famid, $indid, $fatherid, $motherid, $sexid, $phenotype) = @record;

		lc $sexid eq 'm' || lc $sexid eq 'male' and $sexid = 1;
		lc $sexid eq 'f' || lc $sexid eq 'female' and $sexid = 2;
		$sexid eq '0' || lc $sexid eq 'unknown' and $sexid = 0;
		$sexid =~ m/^[012]$/ or confess "Error: family $famid individual $indid does not have a valid sex identifier (1 for male, 2 for female and 0 for unknown): <$sexid>";
		if ($qt) {
			if ($phenotype ne 'NA') {
				$phenotype =~ m/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ or confess "Error: family $famid individual $indid does not have a valid phenotype as floating point number for quantitiatve trait: <$phenotype>";
			}
			if ($phenotype eq '-9') {
				print STDERR "WARNING: The phenotype for family $famid individual $indid is labeled as -9 in $pedheaderfile, and will be treated as unknown\n";
				$phenotype = 'NA';
			}
		} else {
			$phenotype eq '-9' and $phenotype = 0;
			$phenotype =~ m/^[012]$/ or confess "Error: family $famid individual $indid does not have a valid binary phenotype indicator (1 for unaffected, 2 for affected or 0 for unknown): <$phenotype>";
		}
		if ($ped->{$famid}{$indid}) {			#same individual occur more than once in pedheaderfile
			confess "Error: family $famid individual $indid occur more than once in pedheaderfile $pedheaderfile";
		} else {
			$ped->{$famid}{$indid} = [$fatherid, $motherid, $sexid, $phenotype, $index];	#given famid and indid, find the information for the individual
			push @ped_index, [$famid, $indid, $sexid, $phenotype];				#given index in file, find the information for the individual
		}
		$index++;
	}
	close (PEDHEADER);
	print STDERR "NOTICE: Finished reading pedigree information for $index individuals from header file $pedheaderfile\n";
	return ($ped, \@ped_index);
}

#check errors in pedigree file, including sex errors for father and mother and so on
sub inspectPedigree {
	my ($ped, $qt) = @_;
	my ($num_fam, $num_ind);
	my @sex = (0, 0, 0);
	my @phenotype = (0, 0, 0);
	my $ped_stat = {male_founder_index=>[], female_founder_index=>[], founder_index=>[], unknownsex_founder_index=>[]};
	$qt and @phenotype = ();									#for quantitative trait, we calculate their means and standard deviations
	my ($fat_not_found, $mot_not_found) = qw/0 0/;

	for my $famid (keys %$ped) {
		$num_fam++;
		for my $indid (keys %{$ped->{$famid}}) {
			$num_ind++;
			$sex[$ped->{$famid}{$indid}[2]]++;						#the [2] element is the sexid information
			if ($qt) {
				push @phenotype, $ped->{$famid}{$indid}[3];				#the [3] element is the phenotype (dependent variable) information
			} else {
				$phenotype[$ped->{$famid}{$indid}[3]]++;
			}
			if ($ped->{$famid}{$indid}[0] eq '0' and $ped->{$famid}{$indid}[1] eq '0') {	#both father and mother are unknown (treated as "founder")
				push @{$ped_stat->{founder_index}}, $ped->{$famid}{$indid}[4];		#the [4] element is the index in pedheaderfile (or GT file)
				if ($ped->{$famid}{$indid}[2] == 1) {
					push @{$ped_stat->{male_founder_index}}, $ped->{$famid}{$indid}[4];
				} elsif ($ped->{$famid}{$indid}[2] == 2) {
					push @{$ped_stat->{female_founder_index}}, $ped->{$famid}{$indid}[4];
				} elsif ($ped->{$famid}{$indid}[2] == 0) {
					push @{$ped_stat->{unknownsex_founder_index}}, $ped->{$famid}{$indid}[4];
				} else {
					confess "Error: unknown sex: sex code can be 1, 2 and 0 only!";
				}
			}
			
			if ($ped->{$famid}{$indid}[2] == 1) {						#this paragraph get the sex information for all samples
				push @{$ped_stat->{male_index}}, $ped->{$famid}{$indid}[4];
			} elsif ($ped->{$famid}{$indid}[2] == 2) {
				push @{$ped_stat->{female_index}}, $ped->{$famid}{$indid}[4];
			} elsif ($ped->{$famid}{$indid}[2] == 0) {
				push @{$ped_stat->{unknownsex_index}}, $ped->{$famid}{$indid}[4];
			} else {
				confess "Error: unknown sex: sex code can be 1, 2 and 0 only!";
			}
			
			if ($ped->{$famid}{$indid}[0]) {
				$ped->{$famid}{$ped->{$famid}{$indid}[0]} or $fat_not_found++;		#father is annotated in PED annotation of offspring, but father PED is not found in PED file
			}
			if ($ped->{$famid}{$indid}[1]) {
				$ped->{$famid}{$ped->{$famid}{$indid}[1]} or $mot_not_found++;
			}
		}
	}
	$ped_stat->{num_fam} = $num_fam;
	$ped_stat->{num_ind} = $num_ind;
	$ped_stat->{num_male} = $sex[1];
	$ped_stat->{num_female} = $sex[2];
	$ped_stat->{num_unknown_sex} = $sex[0];
	$ped_stat->{num_founder} = scalar (@{$ped_stat->{founder_index}});
	print STDERR "NOTICE: Current pedigree contains $num_fam families with $num_ind individuals\n";
	print STDERR "        Sex: $sex[1] are male, $sex[2] are female, $sex[0] are of unknown sex\n";
	print STDERR "        Founder: $ped_stat->{num_founder} are founders, ", $num_ind-$ped_stat->{num_founder}, " are non-founders\n";
	if ($qt) {
		print STDERR "        Phenotype (quantitative trait): mean is ", sprintf ("%.4g", mean (\@phenotype)), ", standard deviation is ", sprintf ("%.4g", sd (\@phenotype)), "\n\n";
	} else {
		$ped_stat->{num_unaffected} = $phenotype[1];
		$ped_stat->{num_affected} = $phenotype[2];
		$ped_stat->{num_unknown_affection} = $phenotype[0];
		print STDERR "        Phenotype (binary trait): $phenotype[1] are unaffected, $phenotype[2] are affected, and $phenotype[0] are unknown\n\n";
	}
	$fat_not_found and print STDERR "WARNING: $fat_not_found offspring have annotations for father in ped_header file, but father is not present in ped_header file\n";
	$mot_not_found and print STDERR "WARNING: $mot_not_found offspring have annotations for mother in ped_header file, but mother is not present in ped_header file\n";
	return $ped_stat;
}

sub mean {
	my ($score) = @_;
	@$score or confess "Error: NO VALUES for calculating mean";
	my $sum;
	for (@$score) {
		$sum += $_;
	}
	return $sum/@$score;
}

sub sd {
	my ($score) = @_;
	my $mean = mean ($score);
	my $sum;
	for my $i (0 .. @$score-1) {
		$sum += ($score->[$i]-$mean)*($score->[$i]-$mean);
	}
	$sum /= (@$score-1);
	return sqrt ($sum);
}
