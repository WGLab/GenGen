#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 87a93fbdbc7fb54dadc30cc526acab51d337dfce $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) $';

our ($verbose, $help, $man);
our ($bimfile, $snptablefile);
our ($outfile, $intype, $outtype, $replacezero, $force, $autoflip, $strandfile);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'intype=s'=>\$intype, 'outtype=s'=>\$outtype, 
	'replacezero'=>\$replacezero, 'force'=>\$force, 'autoflip'=>\$autoflip, 'strandfile=s'=>\$strandfile) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($bimfile, $snptablefile) = @ARGV;
$intype or print STDERR "NOTICE: The default --intype of 'ilmn12' is assumed for bimfile $bimfile\n" and $intype = 'ilmn12';
$outtype or print STDERR "NOTICE: The default --outtype of 'dbsnp' is assumed as output format\n" and $outtype = 'dbsnp';

$intype eq 'top12' and $intype = 'ilmn12';
$intype eq $outtype and pod2usage ("Error in argument: the --intype $intype cannot be identical to --outtype $outtype");
$intype =~ m/^(ilmn12|ilmnab|top|dbsnp|forward)$/ or pod2usage ("Error in argument: the --intype $intype cannot be procesed (supported intype are ilmn12, ilmnab, top, dbsnp, forward)!");
$autoflip and $intype eq 'dbsnp' || pod2usage ("Error in argument: the --autoflip argument can be used only when --intype is 'dbsnp'");

$outtype eq 'forward' and defined $strandfile || pod2usage ("Error in argument: please specify --strandfile when '-outtype' is 'forward'");
$intype eq 'forward' and defined $strandfile || pod2usage ("Error in argument: please specify --strandfile when '-intype' is 'forward'");

defined $strandfile and $outtype eq 'forward' || $intype eq 'forward' || pod2usage ("Error in argument: the --strandfile is supported only when '-intype' or '-outtype' is 'forward'");

my ($ilmn_strand, $cust_strand, $indel, $no_anno) = readSNPTableFile ($snptablefile);
my $strandinfo;
$strandfile and $strandinfo = readStrandfile ($strandfile, $ilmn_strand);

writeBim ($bimfile, $ilmn_strand, $cust_strand, $outfile, $indel, $no_anno, $strandinfo);




sub readStrandfile {
	my ($strandfile, $ilmn_strand) = @_;
	my (%strandinfo);
	my ($pos, $neg) = (0, 0);
	open (STRAND, $strandfile) or die "Error: cannot read strandfile $strandfile: $!\n";
	while (<STRAND>) {
		s/[\r\n]+$//;
		my @field = split (/\t/, $_);
		$field[1] eq '+' or $field[1] eq '-' or die "Error: invalid record found in strandfile $strandfile: <$_>\n";
		$ilmn_strand->{$field[0]} or next;		#skip markers not useful for the BIM file
		$strandinfo{$field[0]} = $field[1];
		$field[1] eq '+' and $pos++;
		$field[1] eq '-' and $neg++;
	}
	print STDERR "NOTICE: Finished reading $pos positive strand and $neg negative strand\n";
	return (\%strandinfo);
}

sub writeBim {
	my ($bimfile, $ilmn_strand, $cust_strand, $outfile, $indel, $no_anno, $strandinfo) = @_;
	my (@indel);
	my %rev = ('A'=>'T', 'C'=>'G', 'G'=>'C', 'T'=>'A');
	my ($without_strand) = (0);
	
	open (BIM, $bimfile) or confess "Error: cannot read from bimfile $bimfile: $!";
	if ($outfile) {
		open (STDOUT, ">$outfile") or confess "Error: cannot write to output file $outfile: $!";
		print STDERR "NOTICE: The new bim file will be written to $outfile ... ";
	}
	
	my $flipcount;
	while (<BIM>) {
		s/[\r\n]+$//;
		my @record = split (/\s+/, $_);
		@record == 6 or confess "Error: invalid record found in bimfile $bimfile (6 tab/space-delimited records expected): <$_>\n";
		if ($indel->{$record[1]}) {
			$verbose and print STDERR "WARNING: the SNP $record[1] in bimfile $bimfile is annotated as indel in SNPTable file (removal from BIM file is recommended)\n";
			push @indel, $record[1];
		}
		if ($no_anno->{$record[1]}) {
			confess "Error: the SNP $record[1] in bimfile $bimfile is annotated as ALLELE UNKNOWN (CNV marker?) in SNPTable file (remove it from BIM file before conversion)\n";
		}
		$ilmn_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNPTable file\n";
		
		my ($allele1, $allele2) = @record[4,5];			#same the allele identity since we'll modify record[4,5] below
		$allele1 and not $allele2 and confess "FATAL ERROR: the minor allele for SNP $record[1] is $allele1 but major allele is a zero allele in BIM file $bimfile\n";
		
		if ($intype eq 'ilmn12' or $intype eq 'ilmnab') {
			if ($intype eq 'ilmnab') {
				$record[4] =~ tr/AB/12/;
				$record[5] =~ tr/AB/12/;
				$record[4] =~ m/^0|1|2$/ and $record[5] =~ m/^0|1|2$/ or confess "Error: invalid allele designation (0/A/B expected based on --intype of ilmnab): <$_>\n";
			} else {
				$record[4] =~ m/^0|1|2$/ and $record[5] =~ m/^0|1|2$/ or confess "Error: invalid allele designation (0/1/2 expected based on --intype of ilmn12): <$_>\n";
			}
			
			if ($record[4] == 0) {
				$record[5] == 1 and $record[4] = 2;
				$record[5] == 2 and $record[4] = 1;
				$record[5] == 0 and ($record[4], $record[5]) = (1, 2);		#randomly assign 1 and 2, since no genotype call is generated for this SNP
			}
			if ($outtype eq 'dbsnp' or $outtype eq 'top') {
				my $a12 = ($outtype eq 'dbsnp' ? $cust_strand->{$record[1]} : $ilmn_strand->{$record[1]}) or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file $snptablefile\n";
				my @a = split (//, $a12);
				if ($record[4] == 1) {
					@record[4,5] = @a[0,1];
				} else {
					@record[4,5] = @a[1,0];
				}
			} elsif ($outtype eq 'forward') {
				my $a12 = $cust_strand->{$record[1]};
				my @a = split (//, $a12);
				
				if (not defined $strandinfo->{$record[1]}) {
					@a = (0, 0);
					$without_strand++;
				} elsif ($strandinfo->{$record[1]} eq '-') {
					@a = ($rev{$a[0]}, $rev{$a[1]});
				} elsif ($strandinfo->{$record[1]} eq '+') {
					1;
				} else {				#for those without strand information, make the genotype as zero as it is undetermined
					@a = (0, 0);
				}
				
				if ($record[4] == 1) {
					@record[4,5] = @a[0,1];
				} else {
					@record[4,5] = @a[1,0];
				}
				
			} elsif ($outtype eq 'ilmnab') {
				$record[4] =~ tr/12/AB/;
				$record[5] =~ tr/12/AB/;
			} elsif ($outtype eq 'ilmn12') {
				$record[4] =~ tr/AB/12/;
				$record[5] =~ tr/AB/12/;
			}
		} elsif ($intype eq 'top') {
			$record[4] =~ m/^[ACTG0]$/ and $record[5] =~ m/^[ACTG0]$/ or confess "Error: invalid alleles found in BIM file (the -intype of 'top' indicates possible alleles as ACGT0): <$_>\n";
			if ($outtype eq 'ilmn12' or $outtype eq 'ilmnab') {
				my $a12 = $ilmn_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file\n";
				my @a = split (//, $a12);
				if ($record[5] eq '0') {
					@record[4,5] = (1, 2);				#randomly assign 1 and 2, since no genotype call is generated for this SNP
				} else {
					if ($record[5] eq $a[1]) {
						$record[4] and $record[4] eq $a[0] || confess "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = (1,2);
					} elsif ($record[5] eq $a[0]) {
						$record[4] and $record[4] eq $a[1] || confess "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = (2,1);
					} else {
						if ($force) {
							warn "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file. Treated as zero allele in output\n";
							@record[4,5] = (0, 0);
						} else {
							confess "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						}
					}
				}
				
				if ($outtype eq 'ilmnab') {
					$record[4] =~ tr/12/AB/;
					$record[5] =~ tr/12/AB/;
				}
			} elsif ($outtype eq 'dbsnp') {
				my $a12 = $ilmn_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file\n";
				my @a = split (//, $a12);
				my $b12 = $cust_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file\n";
				my @b = split (//, $b12);
				
				if ($record[5] eq '0') {
					$record[4] and confess "Error: the SNP $record[1] major allele is annotated as 0 but minor allele is $record[4] in BIM file $bimfile\n";
					@record[4,5] = @b;				#randomly assign alleles, since no genotype call is generated for this SNP
				} else {
					if ($record[5] eq $a[1]) {
						$record[4] and $record[4] eq $a[0] || confess "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = @b[0,1];
					} elsif ($record[5] eq $a[0]) {
						$record[4] and $record[4] eq $a[1] || confess "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = @b[1,0];
					} else {
						if ($force) {
							warn "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file; Treated as zero allele in output\n";
							@record[4,5] = (0, 0);
						} else {
							confess "Error: the SNP $record[1] top allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						}
					}
				}
			}
		} elsif ($intype eq 'dbsnp' or $intype eq 'forward') {
			$record[4] =~ m/[ACGT0]$/ and $record[5] =~ m/^[ACGT0]$/ or confess "Error: invalid alleles found in BIM file (the -intype of 'dbsnp' indicates possible alleles as ACGT0): <$_>\n";

			if ($intype eq 'forward') {		#convert forward to dbSNP first
				$strandinfo->{$record[1]} or die "Error: strand information for $record[1] is not found in strandfile\n";
				if ($strandinfo->{$record[1]} eq '-') {
					$record[4] = complement($record[4]);
					$record[5] = complement($record[5]);
				} elsif ($strandinfo->{$record[1]} eq '+') {
					1;
				}
			}

			if ($outtype eq 'ilmn12' or $outtype eq 'ilmnab') {
				my $a12 = $cust_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file\n";
				my @a = split (//, $a12);
				if ($record[5] eq '0') {
					$record[4] and confess "Error: the SNP $record[1] major allele is annotated as 0 but minor allele is $record[4] in BIM file $bimfile\n";
					@record[4,5] = (1,2);			#this SNP has no genotyping call in all subjects
				} else {
					
					#when intype is dbSNP, check the allele correspondence and flip alleles when necessary
					if ($autoflip and not pairnt ($record[4], $record[5])) {
						if ($record[4] eq '0') {
							if ($record[5] eq complement($a[0]) or $record[5] eq complement($a[1])) {
								$record[5] = complement ($record[5]);
								$flipcount++;
							}
						} elsif ($record[4] eq complement($a[0]) and $record[5] eq complement($a[1]) or $record[4] eq complement($a[1]) and $record[5] eq complement($a[0]) ) {
							@record[4,5] = (complement($record[4]), complement ($record[5]));
							$flipcount++;
						}
					}
					
					if ($record[5] eq $a[0]) {
						$record[4] and $record[4] eq $a[1] || confess "Error: the SNP $record[1] forward strand allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = (2,1);
					} elsif ($record[5] eq $a[1]) {
						$record[4] and $record[4] eq $a[0] || confess "Error: the SNP $record[1] forward strand allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = (1,2);
					} else {
						if ($force) {
							warn "Error: the SNP $record[1] forward strand allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file; Treated as zero allele in output\n";
							@record[4,5] = (0, 0);
						} else {
							die "Error: the SNP $record[1] forward strand allele is annotated as '@a' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						}
					}
				}
				if ($outtype eq 'ilmnab') {
					$record[4] =~ tr/12/AB/;
					$record[5] =~ tr/12/AB/;
				}	
			} elsif ($outtype eq 'top') {
				my $a12 = $ilmn_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file\n";
				my @a = split (//, $a12);
				my $b12 = $cust_strand->{$record[1]} or confess "Error: the SNP $record[1] in bimfile $bimfile is not annotated in SNP Table file\n";
				my @b = split (//, $b12);
				if ($record[5] eq '0') {
					$record[4] and confess "Error: the SNP $record[1] major allele is annotated as 0 but minor allele is $record[4] in BIM file $bimfile\n";
					@record[4,5] = @a;
				} else {
					if ($record[5] eq $b[0]) {
						$record[4] and $record[4] eq $b[1] || confess "Error: the SNP $record[1] forward strand allele is annotated as '@b' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = @a[1,0];
					} elsif ($record[5] eq $b[1]) {
						$record[4] and $record[4] eq $b[0] || confess "Error: the SNP $record[1] forward strand allele is annotated as '@b' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						@record[4,5] = @a[0,1];
					} else {
						if ($force) {
							warn "Error: the SNP $record[1] forward strand allele is annotated as '@b' in snptable file, but is shown as '@record[4,5]' in BIM file; Treated as zero allele in output\n";
							@record[4,5] = (0, 0);
						} else {
							confess "Error: the SNP $record[1] forward strand allele is annotated as '@b' in snptable file, but is shown as '@record[4,5]' in BIM file\n";
						}
					}
				}
			}
		} else {
			pod2usage ("Error in argument: the -intype of $intype is not supported by the program yet!");
		}
		
		if (not $replacezero) {					#if --replacezero is not set, then do not automatically change the unobserved allele to the actual allele
			$allele1 eq '0' and $record[4] = 0;
			$allele2 eq '0' and $record[5] = 0;
		}
		
		print join("\t", @record), "\n";
	}
	close (BIM);
	print STDERR "Done\n";
	$flipcount and print STDERR "NOTICE: the auto-flipping function was performed on $flipcount SNPs in the input BIM file\n";
	@indel and print STDERR "WARNING: ${\(scalar @indel)} indels (examples:", join (',', splice (@indel, 0, 3)), ") are included in the output file\n";
	$without_strand and print STDERR "WARNING: $without_strand markers have no strand information and their genotypes are treated as zero\n";
}

sub pairnt {
	@_ == 2 or die "Error: subroutine pairnt() requires two arguments: <@_>\n";
	my ($nt1, $nt2) = @_;
	if ($nt1 eq 'A' and $nt2 eq 'T' or $nt1 eq 'T' and $nt2 eq 'A' or $nt1 eq 'C' and $nt2 eq 'G' or $nt1 eq 'G' and $nt2 eq 'C') {
		return 1;
	} else {
		return 0;
	}
}

sub complement {
	my ($nt) = @_;
	$nt =~ tr/ACGT/TGCA/;
	return $nt;
}

sub readSNPTableFile {
	my ($snptablefile) = @_;
	my ($name_index, $snp_index, $ilmn_index, $cust_index);
	my @record;
	my %rev = ('A'=>'T', 'C'=>'G', 'G'=>'C', 'T'=>'A');
	my (%ilmn_strand, %cust_strand, @indel, @no_anno);
	
	open (SNPTABLE, $snptablefile) or confess "Error: cannot read from SNP Table file $snptablefile: $!";
	print STDERR "NOTICE: Reading SNP Table file $snptablefile ...";
	
	$_ = <SNPTABLE>;
	s/[\r\n]+$//;
	@record = split (/\t/, $_);
	@record >= 4 or confess "Error: invalid record found in SNP Table file $snptablefile (at least 4 tab-delimited fields expected): <$_>\n";
	for my $i (0 .. @record-1) {
		$record[$i] eq "Name" and $name_index = $i;
		$record[$i] eq "SNP" and $snp_index = $i;
		$record[$i] eq "ILMN Strand" and $ilmn_index = $i;
		$record[$i] eq "Customer Strand" and $cust_index = $i;
	}
	defined $name_index and defined $snp_index and defined $ilmn_index and defined $cust_index or confess "Error: cannot find Name, SNP, ILMN Strand or Customer Strand column in the first line of SNP Table file $snptablefile: <$_>\n";
	
	while (<SNPTABLE>) {
		s/[\r\n]+$//;
		@record = split (/\t/, $_);
		@record >= 4 or confess "Error: invalid record found in SNP Table file $snptablefile (at least 4 tab-delimited fields expected): <$_>\n";
		my ($name, $snp, $ilmn, $cust) = @record[$name_index, $snp_index, $ilmn_index, $cust_index];
		$ilmn = uc $ilmn;		#sometimes is is annotated as Bot and Top
		$cust = uc $cust;
		if ($snp eq '[D/I]' or $snp eq '[I/D]') {		#for example, rs17838100 is an insertion/deletion polymorphism
			$ilmn =~ m/^[PM]$/ and $cust =~ m/^[PM]$/ or confess "Error: invalid ILMN/Customer Strand found for insertiion/deletion ('P' or 'M' expected): <$ilmn> or <$cust> in <$_>\n";
			push @indel, $name;
			$snp =~ s/[\[\]\/]//g;
			$ilmn_strand{$name} = $snp;
			$cust_strand{$name} = $snp;
			next;
		}
		if ($snp eq '[N/A]') {
			push @no_anno, $name;
			next;
		}
		if ($ilmn eq 'P' or $cust eq 'P') {
			push @no_anno, $name;
			next;
		}
	
		$ilmn =~ m/^TOP|BOT$/ or confess "Error: invalid ILMN Strand found ('TOP' or 'BOT' expected): <$ilmn> in <$_>\n";
		$cust =~ m/^TOP|BOT$/ or confess "Error: invalid Customer Strand found ('TOP' or 'BOT' expected): <$cust> in <$_>\n";
		$snp =~ m#^\[(\w)/(\w)\]$# or confess "Error: Invalid SNP allele designation found: <$snp>\n";
		my ($a1, $a2) = ($1, $2);
		
		if ($ilmn eq 'TOP') {
			if ($cust eq 'TOP') {
				$ilmn_strand{$name} = $a1 . $a2;
				$cust_strand{$name} = $a1 . $a2;
			} elsif ($cust eq 'BOT') {
				$ilmn_strand{$name} = $a1 . $a2;
				$cust_strand{$name} = $rev{$a1} . $rev{$a2};
			} else {
				confess "Error: Neither TOP nor BOT is specified for the Customer strand: <$cust>\n";
			}
		} elsif ($ilmn eq 'BOT') {
			if ($cust eq 'TOP') {
				$ilmn_strand{$name} = $rev{$a1} . $rev{$a2};
				$cust_strand{$name} = $rev{$a1} . $rev{$a2};
			} elsif ($cust eq 'BOT') {
				$ilmn_strand{$name} = $rev{$a1} . $rev{$a2};
				$cust_strand{$name} = $a1 . $a2;
			} else {
				confess "Error: Neither TOP nor BOT is specified for the Customer strand: <$cust>\n";
			}
		} else {
			confess "Error: Neither TOP nor BOT is specified for the ILMN strand: <$ilmn>\n";
		}
	}
	close (SNPTABLE);
	print STDERR " Done with ${\(scalar keys %ilmn_strand)} SNPs\n";
	@indel and print STDERR "NOITCE: ${\(scalar @indel)} insertion/deletion polymorphism are annotated in $snptablefile (examples:", join (',', splice (@indel, 0, 3)), ")\n";
	@no_anno and print STDERR "WARNING: ${\(scalar @no_anno)} ALLELE UNKNOWN markers are annotated in $snptablefile (example:", join (',', splice (@no_anno, 0, 3)), ")\n";
	
	my %indel = map {$_, 1} @indel;
	my %no_anno = map {$_, 1} @no_anno;
	return (\%ilmn_strand, \%cust_strand, \%indel, \%no_anno);
}


=head1 SYNOPSIS

 convert_bim_allele.pl [arguments] <bimfile> <snptablefile>

 Optional arguments:
        -h, --help                   		print help message
        -m, --man                    		print complete documentation
        -v, --verbose                 		use verbose output
            --intype <ilmn12|ilmnab|top|dbsnp|forward>	specify input type
            --outtype <ilmn12|ilmnab|top|dbsnp|forward>	specify output type
            --force				force execution when allele mismatch is detected (output zero allele)
            --replacezero			replace the zero (unobserved) allele by known allele (default: output zero allele)
            --outfile <file>			specify output file name (default: STDOUT)
            --autoflip				automatically flip strand where appropriate
            --strandfile <file>			specify a file with dbSNP strand information

 Function: convert a BIM file so that the allele names are changed to the 
 specified outtype.
 
 Notes: ilmn12=Illumina AB allele designation coded as 1 and 2
        ilmnab=Illumina AB allele designation
        top   =Illumina Top Allele call
        dbsnp =dbSNP (not necessarily forward) allele designation
        forward=forward strand in reference genome

 Example: convert_bim_allele.pl old.bim hh550_610.snptable -outfile new.bim -intype ilmnab -outtype top
          convert_bim_allele.pl old.bim hh550_610.snptable -outfile new.bim -intype top -outtype dbsnp

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--intype>

specify the intype of allele names in the BIM file. The type can be ilmn12 (1/2 
designation), ilmnab (A/B designation), top (TOP allele) and dbsnp (allele in 
forward strand), with the ilmn12 being the default option.

=item B<--outtype>

specify the outtype of allele names in the generated new BIM file. The dbsnp is 
the default output option.

==item B<--force>

force the program to continue (but print out a warning message) when allele 
mismatch (between input BIM file and the SNPTABLE file) is detected, so that 
users can inspect all the mismatches based on the warning message. The alleles 
with strand mismatch will be printed out as 0. Normally, when this occurs, the 
program will throw out an error and stop; this argument helps user to manually 
inspect all SNPs showing these problems (usually caused by the use of different 
genome assemblies).

=item B<--replacezero>

replace unobserved alleles in BIM file by the predicted actual allele from the 
SNPTABLE file. When both alleles are unobserved (SNP has no genotyping call at 
all), two alleles are randomly assigned. When one allele is unobserved 
(annotated as zero in the BIM file), it will be filled in based on the other 
allele.

=item B<--outfile>

specify the output file. When the argument is not specified, the output will be 
printed to the STDOUT.

=item B<--autoflip>

automaticallly flip the allele strands, when dbSNP is the input type, and when 
the allele designation for input differ from those in the SNPTableFile.

=item B<--strandfile>

specify a file with dbSNP strand information. This file can be generated by the 
user, for example, first download 
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp138.txt.gz then "cut 
-f 5,7 snp138.txt > strandfile".

=back

=head1 DESCRIPTION

This program is used to convert allele names in a BIM file, based on specified 
-intype and -outtype arguments in the command line. BIM file is a file generated 
by the PLINK software when creating a binary PED file with --make-bed argument. 
The BIM file contains information, including allele names, for each SNP.

One common use of the allele name conversion is to do imputation analysis, since 
such analysis typically requires that allele names follow forward strand 
designation of the dbSNP data.


=head2 Annotation mistake in HapMap or dbSNP

Since human genome assembly is not perfect, there is no doubt that annotation 
mistakes exist in databases. Generally speaking, only dozens of mistakes should 
be present for each chromoome.

The mistakes can be several types:

1. The forward allele is annotated errananeouly in HapMap haplotype phasing file. 

2. Discordant strand designation in different sequencing builds, or different 
SNP submission batches. 

=back

=cut
