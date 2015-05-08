#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2013-04-22 21:03:06 -0700 (Mon, 22 Apr 2013) $';

our ($verbose, $help, $man);
our ($probfile, $infofile);
our ($prefix, $rsq, $mapfile, $legendfile, $keepfile, $blocksize);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'prefix=s'=>\$prefix, 'rsq=f'=>\$rsq, 'mapfile=s'=>\$mapfile, 'legendfile=s'=>\$legendfile, 'keepfile=s'=>\$keepfile, 'blocksize=i'=>\$blocksize) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error: <inputfile> missing");

($probfile, $infofile) = @ARGV;

#add the path of the executable to enviroment variable, since later on we'll need to call combine_snptest.pl program
my $path = $0;
$path =~ s/[^\\\/]+$//;
$ENV{PATH} = "$path:$ENV{PATH}";


defined $rsq or $rsq = 0.3;
$prefix ||= $probfile;
$legendfile or pod2usage ("Error in argument: please specify the --legendfile argument");
$blocksize ||= 100_000_000;


my ($keep, $legend, $add_legend, $allprob, $allfam, $allind, $info);
$keepfile and $keep = readKeepfile ($keepfile);
$mapfile and $add_legend = readMap ($mapfile);
$legend = readLegend ($legendfile);


$info = readInfoFile ($infofile, "$prefix.exclude", $rsq, $legend, $add_legend);

my ($block, $startline) = (1, 1);
while (1) {
	my ($allprob, $allfam, $allind, $countline, $maxmem) = readProbFile ($probfile, $keep, $block, $startline);
	if (not $allprob) {			#last round has reached maximum memory, but the file has already been finished reading
		$block--;
		last;
	}
	print STDERR "NOTICE: Finished prob for ", scalar (@$allprob), " subjects\n";
	
	writeGenFile ("$prefix.block$block.gen", $allprob, $info, $legend, $add_legend, $keep);
	writeSampleFile ("$prefix.block$block.sample", $allfam, $allind, $keep);

	$maxmem or last;			#did not reach maximum memory, indicating the prob file has already been finished reading

	$block++;
	$startline = $countline+1;
	
}

my $command = "combine_snptest.pl";
for my $i (1 .. $block) {
	$command .= " $prefix.block$i";
}
$command .= " -prefix $prefix";
print STDERR "NOTICE: Executing system command $command\n";
system ($command);

for my $i (1 .. $block) {
	unlink ("$prefix.block$i.gen", "$prefix.block$i.sample");
}


sub writeGenFile {
	my ($outfile, $allprob, $info, $legend, $add_legend, $keep) = @_;
	
	open (GEN, ">$outfile") or confess "Error: cannot write to output file $outfile\n";
	for my $i (0 .. @$info-1) {
		print GEN join (" ", $info->[$i][0], @{$info->[$i]});
		for my $j (0 .. @$allprob-1) {
			defined $allprob->[$j][$i*2] or warn "Probability not defined for i=$i j=$j position=2*i\n";
			defined $allprob->[$j][$i*2+1] or warn "Probability not defined for i=$i j=$j position=2*i+1\n";
			print GEN " ", $allprob->[$j][$i*2], " ", $allprob->[$j][$i*2+1], " ", sprintf("%.3f", 1-$allprob->[$j][$i*2]-$allprob->[$j][$i*2+1]);
		}
		print GEN "\n";
	}
	close (GEN);
}

sub readInfoFile {
	my ($infofile, $excludefile, $rsq, $legend, $add_legend) = @_;
	
	my ($skip_rsq, $countline, $position, @info) = (0 , 0);
	open (INFO, $infofile) or confess "Error: cannot read from infofile $infofile: $!\n";
	open (EXCLUDE, ">$excludefile") or confess "Error: cannot write to excludefile $excludefile: $!\n";
	$_ = <INFO>;
	$_ =~ m/^SNP\s+Al1\s+Al2\s+Freq1\s+MAF\s+Quality\s+Rsq/ or confess "Error: invalid first line in info file: <$_>\n";
	while (<INFO>) {
		$countline++;
		m/^(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)/ or confess "Error: invalid record found in info file: <$_>\n";
		if ($legend->{$1}) {
			$position = $legend->{$1};
		} elsif ($add_legend and $add_legend->{$1}) {
			$position = $add_legend->{$1};
		} else {
			confess "Error: cannot find position information in legendfile for marker $1 in infofile $infofile\n";
		}
		push @info, [$1, $position, $2, $3];
		if ($4 < $rsq) {
			$skip_rsq++;
			print EXCLUDE $1, "\n";
		}
	}
	close (EXCLUDE);
	close (INFO);
	
	print STDERR "NOTICE: Finished reading $countline markers from infofile $infofile\n";
	print STDERR "NOTICE: Total of $skip_rsq markers not passing --rsq threshold were written to $excludefile\n";
	return (\@info);
}

sub readProbFile {
	my ($probfile, $keep, $block, $startline) = @_;
	my (@allprob, @allfam, @allind);
	my ($countline, $memsize) = (0, 0);
	my ($maxmem);
	
	print STDERR "NOTICE: Reading prob file $probfile for block $block with starting line $startline ...";
	open (PROB, $probfile) or confess "Error: cannot read from probfile $probfile: $!\n";
	while (<PROB>) {
		$countline++;
		$countline < $startline and next;
		s/[\r\n]+$//;
		m/^(\S+)\->(\S+)\s+ML_PROB\s+(.+)/ or confess "Error: invalid record in probfile: <$_>\n";
		my ($famid, $indid, $prob) = ($1, $2, $3);
		if ($keep) {
			$keep->{$famid}{$indid} or next;
		}
		push @allprob, [split(/\s+/, $prob)];
		push @allfam, $famid;
		push @allind, $indid;
		$memsize += length ($prob);
		
		if ($memsize >= $blocksize) {		#100MB data requires 2GB memory after splitting
			$maxmem++;
			last;
		}
		$verbose and $countline =~ m/00$/ and print STDERR "NOTICE: countline=$countline, memsize=$memsize\n";
	}
	close (PROB);
	if (not @allprob) {				#nothing to read, already finished reading the file
		print STDERR " Nothing read\n";
		return (undef);
	}
	print STDERR "Done with ", scalar (@allprob), " subjects and ", @{$allprob[0]}/2, " markers\n";
	return (\@allprob, \@allfam, \@allind, $countline, $maxmem);
}

sub writeSampleFile {
	my ($outfile, $allfam, $allind, $keep) = @_;
	my ($countsample) = (0);
	open (SAMPLE, ">$outfile") or confess "Error: cannot write to output file $outfile: $!\n";
	print SAMPLE "ID_1 ID_2 missing\n";
	print SAMPLE "0 0 0\n";
	for my $i (0 .. @$allfam-1) {
		if ($keep) {
			$keep->{$allfam->[$i]}{$allind->[$i]} or next;
		}
		print SAMPLE $allfam->[$i], " ", $allind->[$i], " ", 0, "\n";
		$countsample++;
	}
	close (SAMPLE);
	print STDERR "NOTICE: Total of $countsample samples were written to samplefile $outfile\n";
}

sub readKeepfile {
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

sub readLegend {
	my ($legendfile) = @_;
	my (%legend);
	my $countline = 0;
	open (LEGEND, $legendfile) or confess "Error: cannot read from LEGEND file $legendfile: $!\n";
	print STDERR "NOTICE: Reading legend file $legendfile ...";
	$_ = <LEGEND>;
	m/^rs\s+position\s+0\s+1/ or confess "Error: invalid header line found (expected: rs position 0 1) in LEGEND file $legendfile: <$_>\n";
	while (<LEGEND>) {
		$countline++;
		m/^(\S+)\s+(\S+)/ and $legend{$1} = $2;
	}
	print STDERR " Done with legend information for $countline markers\n";
	close (LEGEND);
	return (\%legend);
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
	close (MAP);
	return (\%add_legend);
}


=head1 SYNOPSIS

 convert_mach2snptest.pl [arguments] <prob-file> <info-file>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	    --prefix <string>		output prefix (default:<geno-file>)
 	    --rsq <float>		rsq threshold to write SNPs to exclude file (default: 0.3)
 	    --legendfile <file>		HapMap legend file used in MACH imputation
 	    --mapfile <file>		MAP file (only useful when imputation results contains markers not annotated in legend file)
 	    --keepfile <file>		a two-column file specifying subjects to output
 	    --blocksize <int>		size of each reading block (default=100,000,000, reduce if out-of-memory)
 	
 Function: convert MACH imputation results to SNPTEST input files (a gen file, a sample file, and an exclude file).
 
 Example: convert_mach2snptest.pl chr22.mlprob chr22_step2.mlinfo -prefix output -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt
          convert_mach2snptest.pl chr22.mlprob chr22_step2.mlinfo -prefix output -legend genotypes_chr22_CEU_r22_nr.b36_fwd_legend.txt -rsq 0.6 -keep caseid

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=item B<--prefix>

specify the output prefix. The suffix of the output files will be gen and 
sample, for the genotype data and the same data, respectively.

=item B<--rsq>

the Rsq threshold to write SNPs to an exclude file. The threshold is reported by 
MACH to flag poorly imputed SNPs, and a default threshold of 0.3 is generally 
recommanded.

=item B<--legendfile>

the HapMap legend file that were used in the MACH imputation. This file is 
useful to infer the information on SNPs, such as chromosome and position 
information.

=item B<--mapfile>

the MAP file used in the imputation. Normally this file is not necessary. 
However, in older versions of MACH, some SNPs that are not in HapMap haplotypes 
are still written to the output by MACH, so the MAP file is requires to figure 
out the chromosome and position information for these SNPs.

=item B<--keepfile>

specify a file containing family and individual identifiers that should be kept 
in the output. This argument is useful to select cases and controls from the 
combined data set.

=item B<--blocksize>

the size of a memory block to process by the program internally. Default value 
is 100000000 (or 100MB). If there is an out-of-memory error reported by the 
system, try to decrease this value.

=back

=head1 DESCRIPTION

This program is used to convert MACH-generated imputation files into SNPTEST 
formats, to be analyzed for association by the SNPTEST software.

=back
                                                                                                                                                                                                                                                                                                                                                                       
