#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-12-17 11:15:14 -0800 (Thu, 17 Dec 2009) $';

our ($verbose, $help, $man);
our ($gen_suffix, $sample_suffix, $keepfile, $prefix);
our (@prefix);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'gen_suffix=s'=>\$gen_suffix, 'sample_suffix=s'=>\$sample_suffix, 'keepfile=s'=>\$keepfile, 'prefix=s'=>\$prefix) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV >= 1 or pod2usage ("Syntax error: <inputfile> missing");

(@prefix) = @ARGV;
$gen_suffix ||= 'gen';
$sample_suffix ||= 'sample';
if (not $prefix) {
	print STDERR "NOTICE: The output prefix is automatically set as $prefix[0].combine (use --prefix to override)\n";
	$prefix = "$prefix[0].combine";
}



my (@fh1, @allsample, $keep);
$keepfile and $keep = readKeepFile ($keepfile);

for my $i (0 .. @prefix-1) {
	my ($fh1, $fh2);
	if ($gen_suffix =~ m/\.gz$/) {
		open ($fh1, "gunzip -c $prefix[$i].$gen_suffix |") or confess "Error: cannot read from gen file $prefix[$i].$gen_suffix: $!\n";
	} else {
		open ($fh1, "$prefix[$i].$gen_suffix") or confess "Error: cannot read from gen file $prefix[$i].$gen_suffix: $!\n";
	}
	push @fh1, $fh1;

	if ($sample_suffix =~ m/\.gz$/) {
		open (FH, "gunzip -c $prefix[$i].$sample_suffix |") or confess "Error: cannot read from gen file $prefix[$i].$gen_suffix: $!\n";
	} else {
		open (FH, "$prefix[$i].$sample_suffix") or confess "Error: cannot read from gen file $prefix[$i].$gen_suffix: $!\n";
	}
	$_ = <FH>;
	m/^ID_1\s+ID_2\s+missing/ or confess "Error: invalid record found (expected: ID_1 ID_2 missing)in header line of sample file $prefix[$i].$sample_suffix: <$_>\n";
	$_ = <FH>;
	m/^0 0 0/ or confess "Error: invalid record found in second line of sample file (0 0 0 expected): <$_>\n";
	while (<FH>) {
		s/[\r\n]+$//;
		m/^(\S+)\s+(\S+)/ or confess "Error: the sample file $prefix[$i].sample should contain at least two columns in each line: <$_>\n";
		push @{$allsample[$i]}, [split (/\s+/, $_)];
	}
	close (FH);
}

print STDERR "NOTICE: Writting genotype data to $prefix.gen\n";
my $countline = 0;
open (GEN, ">$prefix.gen") or confess "Error: cannot write to combined gen file $prefix.gen: $!\n";
MAIN: while (1) {
	my ($line1, @line1, $snpid, $rsid, $position, $a1, $a2, @prob);
	my ($refsnpid, $refrsid, $refposition, $refa1, $refa2);
	for my $i (0 .. @prefix-1) {
		my $switch;
		$line1 = readline ($fh1[$i]);
		$line1 or last MAIN;				#finished reading all lines
		$line1 =~ s/[\r\n]+$//;
		($snpid, $rsid, $position, $a1, $a2, @prob) = split (/\s+/, $line1);
		$i or ($refsnpid, $refrsid, $refposition, $refa1, $refa2) = ($snpid, $rsid, $position, $a1, $a2);
		@prob % 3 == 0 or confess "Error: invalid record found in $prefix[$i].gen: prob scores should be in multiples of 3\n";
		if ($snpid ne $refsnpid or $rsid ne $refrsid or $position ne $refposition) {
			confess "Error: reference ($prefix[0]) has SNP annotation of ($refsnpid, $refrsid, $refposition, $refa1, $refa2) but current file ($prefix[$i]) has annotation of ($snpid, $rsid, $position, $a1, $a2)\n";
		}
		if ($a1 ne $refa1 or $a2 ne $refa2) {
			$a1 eq $refa2 and $a2 eq $refa1 or confess "Error: reference ($prefix[0]) has SNP allele annotation ($refa1, $refa2) that is different from current file $prefix[$i] ($a1, $a2)\n";
			$switch++;
		}
		$i or print GEN join (" ", $snpid, $rsid, $position, $a1, $a2);
		for my $j (0 .. @prob/3-1) {
			if ($keep) {
				$keep->{$allsample[$i]->[$j][0]}{$allsample[$i]->[$j][1]} or next;
			}
			if ($switch) {
				print GEN " ", $prob[$j*3+2], " ", $prob[$j*3+1], " ", $prob[$j*3];
			} else {
				print GEN " ", $prob[$j*3], " ", $prob[$j*3+1], " ", $prob[$j*3+2];
			}
		}
		$i == @prefix-1 and print GEN "\n";
		$i or $countline++;
		$i or $countline % 10000 == 0 and print "NOTICE: Processing line $countline\n";
	}
}
print STDERR "NOTICE: Finished writting genotype data to $prefix.gen with $countline markers\n";

$countline = 0;
print STDERR "NOTICE: Writting sample data to $prefix.sample\n";
open (SAMPLE, ">$prefix.sample") or confess "Error: cannot write to combined sample file $prefix.sample: $!\n";
print SAMPLE "ID_1 ID_2 missing\n";
print SAMPLE "0 0 0\n";
for my $i (0 .. @allsample-1) {
	for my $j (0 .. @{$allsample[$i]}-1) {
		if ($keep) {
			$keep->{$allsample[$i]->[$j][0]}{$allsample[$i]->[$j][1]} or next;
		}
		print SAMPLE join (" ", @{$allsample[$i]->[$j]}), "\n";
		$countline++;
	}
}
print STDERR "NOTICE: Finished writting sample data to $prefix.sample with $countline samples\n";






sub readKeepFile {
	my ($keepfile) = @_;
	my $keep;
	my $countline = 0;
	open (KEEP, $keepfile) or confess "Error: cannot read from keepfile $keepfile: $!\n";
	print STDERR "NOTICE: Reading keepfile $keepfile ...";
	while (<KEEP>) {
		$countline++;
		m/^(\S+)\s+(\S+)$/ or confess "Error: invalid record found in keepfile (2 fields expected): <$_>\n";
		$keep->{$1}{$2}++;
	}
	print " Done with $countline subjects to keep in output file\n";
	close (KEEP);
	return $keep;
}

=head1 SYNOPSIS

 combine_snptest.pl [arguments] <prefix1 | ...>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --prefix <string>		output prefix (default:<combine_prefix1>)
 	    --keepfile <file>		a two-column file specifying subjects to output
 	    --gen_suffix <string>	gen file suffix (default: gen)
 	    --sample_suffix <string>	sample file suffix (default: sample)
 	
 Function: combine multiple SNPTEST files (*.gen and *.sample) into one single file
 
 Example: combine_snptest.pl file1 file2 -keep caseid.keep -prefix combine
          #the above command combines file1.gen, file1,sample, file2.gen, file2.sample, and output combine.gen and combine.sample based on subject in caseid.keep file
          
          combine_snptest.pl file1 -keep caseid.keep -prefix caseonly
          #the above command retrieve a subset of subjects (based on caseid.keep file) from file1.gen and file2.sample, then write to new file caseonly.gen and caseonly.sample

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--prefix>

specify the prefix of output file names

=item B<--keepfile>

specify a file containing the subjects to be written to output files

==item B<--gen_suffix>

the suffix of the genotype file (default: gen)

=item B<--sample_suffix>

the suffix of the sample file (default: sample)

=back

=head1 DESCRIPTION

This program is used to combine multiple files in SNPTEST formats into one 
single file, to be analyzed by the SNPTEST software for association tests. 
Additionall, the program can also extract a specific list of subjects and write 
them to new files.

=back
                                                                                                                                                                                                                                                                                                                                                                       
