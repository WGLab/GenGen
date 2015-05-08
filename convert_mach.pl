#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2010-08-11 14:05:50 -0700 (Wed, 11 Aug 2010) $';

our ($verbose, $help, $man);
our ($genofile, $infofile, $qcfile, $prefix, $pedfile, $legendfile, $chr, $rsq, $mapfile, $qc_threshold);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'prefix=s'=>\$prefix, 'pedfile=s'=>\$pedfile, 'legendfile=s'=>\$legendfile, 'chr=s'=>\$chr, 'rsq=f'=>\$rsq, 'mapfile=s'=>\$mapfile,
	'qc_threshold=f'=>\$qc_threshold) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 3 or pod2usage ("Syntax error: <inputfile> missing");

($genofile, $infofile, $qcfile) = @ARGV;
$prefix ||= $genofile;
$legendfile or pod2usage ("Error in argument: please specify --legendfile argument");

$rsq ||= 0.3;
$qc_threshold ||= 0.9;
if (not $chr) {
	if ($legendfile =~ m/chr(\d+)/) {
		print STDERR "NOTICE: The -chr argument is automatically set as $1 based on the -legendfile argument\n";
		$chr = $1;
	} else {
		pod2usage ("Error in argument: please specify --chr argument");
	}
}


my ($legend, $add_legend, $pedinfo, @marker);



$legend = readLegend ($legendfile);
$mapfile and $add_legend = readMap ($mapfile);
$pedinfo = readPedInfo ($pedfile);

writeMap ($prefix, $infofile, $legend, $add_legend, $chr);
writePed ("$prefix.ped", $genofile, $qcfile, $qc_threshold);


sub writePed {
	my ($outfile, $genofile, $qcfile, $qc_threshold) = @_;
	open (PED, ">$outfile") or confess "Error: cannot write to PED file $outfile";
	
	my ($countline, $count1, $count2, $genofh, $qcfh);
	open ($genofh, $genofile) or confess "Error: cannot read from genofile $genofile: $!\n";
	open ($qcfh, $qcfile) or confess "Error: connot read from qcfile $qcfile: $!\n";
	while (1) {
		$_ = readline ($genofh);
		$_ or last;
		
		$countline++;
		my ($famid, $indid, $geno, @geno);
		s/[\r\n]+$//;
		m/^(\S+)\->(\S+)\s+ML_GENO\s+(.+)/ or confess "Error: invalid record (famid, indid, ML_GENO expected): <$_>\n";
		($famid, $indid, $geno) = ($1, $2, $3);
		#$geno =~ s#/##g;
		@geno = split (/\s+/, $geno);
		
		my ($qc, $qcfamid, $qcindid, $qcscore, @qcscore);
		$qc = readline ($qcfh);
		$qc or confess "Error: Unable to read QC line from $qcfile after reading genotypes from $genofile\n";
		$qc =~ s/[\r\n]+$//;
		$qc =~ m/^(\S+)\->(\S+) ML_QC\s+(.+)/ or confess "Error: invalid record (famid, indid, ML_QC expected) at line $countline: <$qc>\n";
		($qcfamid, $qcindid, $qcscore) = ($1, $2, $3);
		$famid eq $qcfamid or confess "Error: genofile $genofile has famid=$famid but qcfile $qcfile has famid=$qcfamid\n";
		$indid eq $qcindid or confess "Error: genofile $genofile has indid=$indid but qcfile $qcfile has indid=$qcindid\n";
		@qcscore = split (/\s+/, $qcscore);
		@geno == @qcscore or confess "Error: unequal number of genotype calls for family $famid individual $indid in $genofile (" . scalar (@geno) . ") versus $qcfile (" . scalar (@qcscore) . ")\n";
		
		$pedinfo->{$famid}{$indid} or confess "Error: the family $famid individual $indid is listed in GENO file $genofile, but not annotated in PED file $pedfile\n";
		print PED $pedinfo->{$famid}{$indid};
		for my $i (0 .. @geno-1) {
			if ($qcscore[$i] >= $qc_threshold) {
				my @allele = split (/\//, $geno[$i]);
				print PED " ", $allele[0], " ", $allele[1];
				$count1++;
			} else {
				print PED " 0 0";
			}
			$count2++;
		}
		print PED "\n";
	}
	close (PED);
	print STDERR "NOTICE: Finished writting the $prefix.ped (total genotype call rate: ", $count1/$count2, ", use --exclude argument in PLINK to filter bad SNPs)\n";
}




sub readMap {
	my ($mapfile) = @_;
	my %add_legend;
	open (MAP, $mapfile) or confess "Error: cannot read from map file $mapfile: $!\n";
	while (<MAP>) {
		s/[\r\n]+$//;
		my @record = split (/\s+/, $_);
		@record == 4 or confess "Error: invalid record found in map file $mapfile (4 fields expected): <$_>\n";
		$add_legend{$record[1]} = $record[3];
	}
	return (\%add_legend);
}

sub writeMap {
	my ($prefix, $infofile, $legend, $add_legend, $chr) = @_;
	open (MAP, ">$prefix.map") or confess "Error: cannot write to MAP file $prefix.map";
	open (EXCLUDE, ">$prefix.exclude") or confess "Error: cannot write to exclude file $prefix.exclude";
	open (INFO, $infofile) or confess "Error: cannot read from infofile $infofile: $!\n";
	$_ = <INFO>;
	$_ =~ m/^SNP\s+Al1\s+Al2\s+Freq1\s+MAF\s+Quality\s+Rsq/ or confess "Error: invalid first line in info file: <$_>\n";
	
	my $add_marker;
	my $countline;
	while (<INFO>) {
		m/^(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/ or confess "Error: invalid record found in info file: <$_>\n";
		if ($legend->{$1}) {
			print MAP "$chr\t$1\t0\t", $legend->{$1}, "\n";
			$2<$rsq and print EXCLUDE $1, "\n";
		} elsif ($add_legend->{$1}) {
			$add_marker++;
			print MAP "$chr\t$1\t0\t", $add_legend->{$1}, "\n";
			$2<$rsq and print EXCLUDE $1, "\n";
		} else {
			confess "Error: the SNP $1 is listed in the INFO file $infofile, but not annotated in the LEGEND file $legendfile ", $mapfile?" or mapfile $mapfile":"", "\n";
		}
		$countline++;
	}
	
	$add_marker and print STDERR "NOTICE: Total of $add_marker markers are not annotated in LEGEND file $legendfile but in MAP file $mapfile\n";
	print STDERR "NOTICE: Finished writting the mapfile $prefix.map (with $countline markers) and the exclude file $prefix.exclude (which contains low-quality SNPs not passing --rsq threshold)\n";
	close (INFO);
	close (EXCLUDE);
	close (MAP);
}


sub readPedInfo {
	my ($pedfile) = @_;
	my $pedinfo = {};
	open (PED, $pedfile) or confess "Error: cannot read from PED file $pedfile: $!\n";
	while (<PED>) {
		m/^((\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+)/ or confess "Error: invalid record found: <$_>\n";
		$pedinfo->{$2}{$3} = $1;
	}
	close (PED);
	return ($pedinfo);
}

sub readLegend {
	my ($legendfile) = @_;
	my (%legend);
	open (LEGEND, $legendfile) or confess "Error: cannot read from LEGEND file $legendfile: $!\n";
	while (<LEGEND>) {
		m/^(\S+)\s+(\S+)/ and $legend{$1} = $2;
	}
	return (\%legend);
}



	
	

=head1 SYNOPSIS

 convert_mach.pl [arguments] <geno-file> <info-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --legendfile <file>		HapMap legend file for the chromosome
 	    --pedinfo <file>		ped information file or PLINK fam file (first 6 columns will be used)
 	    --prefix <string>		output prefix (default:<geno-file>
 	    --rsq <float>		rsq threshold to define high-qulaity SNP and write to extract file (default: 0.3)
 	    --qc_threshold <float>	genotype posterior probability threshold to output (default: 0.9)
 	    --mapfile <file>		optional file that provides coordinates for SNPs (if not in legend file)
 	
 Function: convert MACH imputation results to standard PED and MAP files for association analysis.
 
 Example: convert_mach.pl chr22_step2.mlgeno chr22_step2.mlinfo chr22_step2.mlqc -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt -ped chr22.ped -prefix chr22_step2

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--legendfile>

the HapMap legend file that were used in the MACH imputation. This file is 
useful to infer the information on SNPs, such as chromosome and position 
information.

=item B<--mapfile>

The map file can be used optionally to provide chromosome and position 
information for SNPs. Depending on program version of MACH and the actual 
parameter settings, some markers in the original input MAP file to MACH may 
be included in the final MACH output, but these markers may not be included in the 
legendfile. This argument addresses this issue, to make sure that Chr and 
position information are available for all markers in MACH output.

=item B<--pedinfo>

specify a file whose first six columns describe each sample, so that the same 
columns will be copied to the output files. A PED file used in imputation, or a 
FAM file generated by PLINK, can be used here.

=item B<--prefix>

specify the output prefix. The suffix of the output files will be PED and MAP.

=item B<--rsq>

the Rsq threshold to write SNPs to an exclude file. The threshold is reported by 
MACH to flag poorly imputed SNPs, and a default threshold of 0.3 is generally 
recommanded.

=item B<--qc_threshold>

the genotype posterior probability threshold to output a genotype call. By 
default the threshold is 0.9, and a value less than that will be written as null 
genotype.

=back

=head1 DESCRIPTION

This program is used to convert MACH-generated imputation files into standard 
PED/MAP format for subsequent associatoin tests.

=back
                                                                                                                                                                                                                                                                                                                                                                       
