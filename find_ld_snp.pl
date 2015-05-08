#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2011-03-18 01:20:57 -0700 (Fri, 18 Mar 2011) $';

our ($verbose, $help, $man);
our ($queryfile, $candfile, $genofile);

GetOptions ('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, ) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 3 or pod2usage ("Syntax error: <inputfile> missing");

($queryfile, $candfile, $genofile) = @ARGV;

my $system_command;
my (@query, %cand);
my $randfile = "rand" . substr (rand()*1000, 0, 3);

open (QUERY, $queryfile) or confess "Error: cannot read from query list file $queryfile: $!\n";
while (<QUERY>) {
	s/[\r\n]+$//;
	m/^(\S+)/ and push @query, $1;
}
close (QUERY);
@query or confess "Error: no query SNP can be read from query file $queryfile\n";
print STDERR "NOTICE: Finished reading ", scalar (@query), " SNPs from query file $queryfile\n";


open (CAND, $candfile) or confess "Error: cannot read from cand list file $candfile: $!\n";
while (<CAND>) {
	s/[\r\n]+$//;
	m/^(\S+)/ and $cand{$1}++;
}
close (CAND);
%cand and print STDERR "NOTICE: Finished reading ", scalar (keys %cand), " SNPs from cand file $candfile\n";

my $command = "plink --noweb --bfile $genofile --ld-snp-list $queryfile --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out $randfile --silent";
print STDERR "NOTICE: Executing command <$command> ... ";
system ($command) and confess "Error: cannot execute PLINK to calculate LD\n";
print STDERR "Done!\n";

my (%ldsnp, %ldval, %lddist);
open (LD, "$randfile.ld") or confess "Error: cannot open the LD file $randfile.ld";
$_ = <LD>;
defined $_ and m/^\s+CHR_A\s+BP_A/ or confess "Error: invalid LD information format from file $randfile.ld\n";
while (<LD>) {
	s/^\s+|\s+$//g;
	my @field = split (/\s+/, $_);
	@field == 7 or confess "Error: invalid format in PLINK LD output (7 columns expected): <$_>\n";
	$cand{$field[5]} or next;		#only examine SNPs in candidate list
	
	defined $ldval{$field[2]} or $ldval{$field[2]} = $field[6];
	defined $ldsnp{$field[2]} or $ldsnp{$field[2]} = $field[5];
	defined $lddist{$field[2]} or $lddist{$field[2]} = abs ($field[1]-$field[4]);
	
	if ($field[6] > $ldval{$field[2]} or $field[6] == $ldval{$field[2]} and abs ($field[1]-$field[4]) < $lddist{$field[2]}) {
		if (length ($field[6]) > 4) {
			$ldval{$field[2]} = sprintf ("%.2f", $field[6]);
		} else {
			$ldval{$field[2]} = $field[6];
		}
		$ldsnp{$field[2]} = $field[5];
		$lddist{$field[2]} = abs ($field[1]-$field[4]);
	}
}

for my $nextsnp (@query) {
	print $nextsnp, "\t";
	if ($ldsnp{$nextsnp}) {
		print $ldsnp{$nextsnp}, "\t", $ldval{$nextsnp}, "\t", $lddist{$nextsnp}, "\n";
	} else {
		print "NONE\tNA\tNA\n";
	}
}
			
print STDERR "NOTICE: Deleting temporary files $randfile.log and $randfile.ld\n";
unlink ("$randfile.log", "$randfile.ld");





	
	

=head1 SYNOPSIS

 find_ld_snp.pl [arguments] <query-SNP-list-file> <candidate-SNP-list-file> <PLINK-binary-prefix>

 Optional arguments:
 	-v, --verbose			use verbose output
 	-h, --help			print help message
 	-m, --man			print complete documentation
 	
 Function: find SNPs in candidate list that are best proxy for SNPs in query list
 
 Example: find_ld_snp.pl querylist humanhap550.snplist hapmap_ceu_r23a
 
 Version: $LastChangedDate: 2011-03-18 01:20:57 -0700 (Fri, 18 Mar 2011) $
 

=head1 OPTIONS

=over 8

=item B<--help>

print a brief help message and exit

=item B<--man>

print the complete manual of how to use the program

=item B<--verbose>

use verbose output

=back

=head1 DESCRIPTION

This program is used to identify SNPs in a particular array (for example, 
Illumina array with 550K markers), that are best proxy for a list of candidate 
SNPs (for example, dozens of SNPs reported in a meta-analysis paper using 
imputation data). Users need to provide a PLINK-formatted file for the purpose 
of calculating LDs. The files need to corresopnd to the population under study, 
for example, if hapmap_ceu_r23a is used, then the program can used only for 
identifying proxy SNPs in subjects of European ancestry.

The query-snp-list-file and candidate-snp-list-file should both be simple text 
files, with one SNP per line. SNPs that are not in the PLINK-binary file will 
obviously not be interrogated. Most users probably want to use HapMap 
populations as reference populations to calculate LD (linkage disequilibrium), 
and these files can be downloaded from PLINK website. For example, CEU genotype 
data is at http://pngu.mgh.harvard.edu/~purcell/plink/dist/hapmap_CEU_r23a.zip. 
Uncompress the downloaded files will generate three files, with the bed, bim and 
fam suffix.

One example is demonstrated below:

	[kaiwang@cc ~/]$ cat snplist.candidate 
	rs9939609
	rs6548238
	rs17782313
	rs10938397
	rs7498665
	rs10838738
	rs11084753
	rs2815752

The candidate file contains eight SNPs, and given a Illumina GWAS data, we want 
to see whether these eight SNPs (or their best proxys) show any significance. Do 
this:

	[kaiwang@cc ~/]$ find_ld_snp.pl snplist.candidate snplist.array hapmap_CEU_r23a
	NOTICE: Finished reading 8 SNPs from query file snplist.candidate
	NOTICE: Finished reading 536636 SNPs from cand file snplist.array
	NOTICE: Executing command <plink --noweb --bfile /home/kaiwang/lib/hapmap/plink1/hapmap_CEU_r23a --ld-snp-list snplist.candidate --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out rand21. --silent> ... Done!
	rs9939609       rs3751812       1       2067
	rs6548238       rs4854344       0.844156        3239
	rs17782313      rs10871777      1       666
	rs10938397      rs12641981      1       2644
	rs7498665       rs7498665       1       0
	rs10838738      rs10838738      1       0
	rs11084753      rs29941 0.638822        12605
	rs2815752       rs2568958       1       47324

So the best proxy for each candidate SNP, as well as their r2 measures are 
printed in the output. The last column is the distance between the index SNP and 
the best proxy SNP.

For questions, comments or bug reports, please contact me at 
kai@openbioinformatics.org.


=cut                                                                                                                                                                                                                                                                                                                                                 
