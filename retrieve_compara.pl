#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2007-11-01 17:11:51 -0700 (Thu, 01 Nov 2007) $';

our ($verbose, $help, $man);
our ($file, $fh, $keep_blank_line, $single_line_output, $prefix_output, $suffix_output, $append_query, $affix_query, $query_line);
our ($species, $outspecies, $latin_species, $intype, $outtype, $host, $user, $port, $password, $comment, $ev);
our ($homotype);

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'file=s'=>\$file, 'species=s'=>\$species,
	'intype=s'=>\$intype, 'outtype=s'=>\$outtype, 'keep_blank_line'=>\$keep_blank_line, 'single_line_output'=>\$single_line_output,
	'prefix_output|pre'=>\$prefix_output, 'suffix_output'=>\$suffix_output, 'append_query'=>\$append_query, 'affix_query'=>\$affix_query, 'host=s'=>\$host, 'user=s'=>\$user, 
	'port=s'=>\$port, 'password=s'=>\$password,
	'homotype=s'=>\$homotype, 'comment_flag=s'=>\$comment, 'ev'=>\$ev) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

#checking the validity of arguments
$file or @ARGV or pod2usage ("Error in argument: please either specify the --file argument or provide identifiers in command line\n");
$file and @ARGV and pod2usage ("Error in argument: please do not specify the --file argument and provide identifiers in command line simultaneously\n");
$species or pod2usage ("Error in argument: please specify the --species argument\n");
$homotype and $homotype eq 'paralog' || $homotype eq 'ortholog' || pod2usage ("Error in argument: the --homotype argument can be only 'homolog' or 'ortholog'\n");

#setting up default arguments
$intype ||= 'geneid';
$outtype ||= 'geneid';
if ($ev) {
	#the --host, --user, --port and --password argument take precedence over the --ev argument
	$host ||= $ENV{ENSEMBLDB} || '';
	$user ||= $ENV{ENSEMBLUSER} || '';
	$port ||= $ENV{ENSEMBLPORT} || '';
	$password ||= $ENV{ENSEMBLPASS} || '';
}
$host ||= 'ensembldb.ensembl.org';
$user ||= 'anonymous';
$port ||= 3306;
$password ||= '';
$| = 1;

#continue checking the validity of arguments
$intype eq 'geneid' and $outtype eq 'geneid' or pod2usage ("Error in argument: --intype and --outtype can be only geneid at this moment");

#loading registry from remote database
$verbose and print STDERR "NOTICE: Loading registry from database\n";
Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host, -user => $user, -port => $port, -pass=>$password);
if (Bio::EnsEMBL::Registry->alias_exists ($species)) {
	$outspecies = Bio::EnsEMBL::Registry->get_alias ($species);
	$latin_species = ucfirst $outspecies;
	$latin_species =~ s/_/ /;
} else {
	confess "Error in argument: the species name <$species> cannot be recognized";
}

#getting database adaptors
#my $db1 = Bio::EnsEMBL::Registry->get_DBAdaptor($inspecies, "core");
my $db2 = Bio::EnsEMBL::Registry->get_DBAdaptor("Multi", "compara");
my $member_adaptor = $db2->get_MemberAdaptor();
my $homology_adaptor = $db2->get_HomologyAdaptor();

#preparing input/output operation to read identifiers
if (defined $file) {
	if ($file eq 'stdin') {
		$fh = \*STDIN;
	} else {
		open ($fh, $file) or confess "Error: cannot read from file $file: $!";
	}
}

while (defined (my $identifier = nextIdentifier ($fh))) {
	$identifier eq '' and writeResult ($identifier) and next;
	$comment and $identifier =~ m/^$comment/ and print "$identifier\n" and next;
	$verbose and print STDERR "NOTICE: Processing identifier <$identifier>\n";
	
	if ($intype eq 'geneid') {
		writeResult ($identifier, geneid2output ($identifier));
	}
}


sub geneid2output {
	my ($identifier) = @_;
	# then you get the homologies where the member is involved
	my $member1 = $member_adaptor->fetch_by_source_stable_id ('ENSEMBLGENE', $identifier);
	not $member1 and print STDERR "WARNING: Unable to retrieve gene <$identifier> in Ensembl-Compara database\n" and return;

	my @output;
	my @homology = @{$homology_adaptor->fetch_all_by_Member ($member1)};	#fetch_by_Member is depreciated and is not an alias to fetch_all_by_Member in version 42
	
	foreach my $homology (@homology) {
	  # You will find different kind of description
	  # UBRH, MBRH, RHS, YoungParalogues
	  # see ensembl-compara/docs/docs/schema_doc.html for more details
	
		#$verbose and print STDERR $homology->description, " ", $homology->subtype,"\n";
		#$homology->print_homology; next;
		$homotype and not $homology->description =~ m/$homotype/i and next;

		my @member_attribute = @{$homology->get_all_Member_Attribute};
		@member_attribute == 2 or confess "Error retrieving paired species members";	#first element is always the query gene, second element is its homolog
		if ($member_attribute[1]->[0]->genome_db->name eq $latin_species) {
			push @output, $member_attribute[1]->[0]->stable_id;
		}
	}
	return @output;
}


#the following subroutine is used to find homologs for an ensembl gene, however it stops working after I upgrade to compara-42 for many genes
sub geneid2output_old {
	my ($identifier) = @_;
	my $member1 = $member_adaptor->fetch_by_source_stable_id ('ENSEMBLGENE', $identifier);
	not $member1 and print STDERR "WARNING: Unable to retrieve gene <$identifier> in Ensembl-Compara database\n" and return;
	
	my @output;
	my (@homology) = @{$homology_adaptor->fetch_all_by_Member_paired_species ($member1, $outspecies)};	#original version is fetch_by_Member_paired_species, which is depreciated as of ENSEMBL-COMPARA 42.2

	@homology or print STDERR "WARNING: Unable to retrieve homolog for identifier <$identifier>\n";

	for my $homology (@homology) {
		my $description = $homology->description;
		if ($homotype) {
			if ($description =~ m/$homotype/i) {
				my @member_attribute = @{$homology->get_all_Member_Attribute};
				@member_attribute == 2 or confess "Error retrieving paired species members";
				push @output, $member_attribute[1]->[0]->stable_id;
			} else {
				$verbose and print STDERR "WARNING: Skipping hit because it is of type <$description>\n";
				next;
			}
		} else {
			my @member_attribute = @{$homology->get_all_Member_Attribute};
			@member_attribute == 2 or confess "Error retrieving paired species members";
			push @output, $member_attribute[1]->[0]->stable_id;
		}
	}
	return @output;
}


#retrieve the next query identifier from the command line or a file or STDIN
sub nextIdentifier {
	my ($fh) = @_;
	my $identifier;
	if ($fh) {
		$identifier = <$fh>;
	} else {
		$identifier = shift @ARGV;
	}
	defined $identifier or return undef;
	$identifier =~ s/[\r\n]+$//;
	$query_line = $identifier;			#query_line is a global variable that store the current line (it is useful for the --append_query argument)
	$comment and $identifier =~ m/^$comment/ and return $identifier;	#this is a comment line, so identifier is returned verbatim
	$identifier =~ m/^\s*(\w[\w:\-\.\/]*)/ or return '';	#identifier is the first word (including colon and dash and dot)
	return $1;
}

#this subroutine print out the database retrieval results. It uses several arguments provided to the program to decide how to format the outputs
sub writeResult {
	my ($identifier, @result) = @_;
	if ($identifier eq '') {		#there is no identifier (maybe the input is an empty line or a line starting with non-alphabetic characters)
		$keep_blank_line and print "\n";
	} elsif (not @result) {			#there is no results for the identifier
		if ($keep_blank_line) {		#still print the prefix and suffix, even if there is no result
			$prefix_output and print $identifier, "\t";
			$suffix_output and print "\t", $identifier;
			$append_query and print "\t", $query_line;
			$affix_query and print $query_line, "\t";
			print "\n";
		}
	} else {				#this is the common situation
		if ($single_line_output) {
			$affix_query and print $query_line, "\t";
			$prefix_output and print $identifier, "\t";
			print join (";", @result);
			$suffix_output and print "\t", $identifier;
			$append_query and print "\t", $query_line;
			print "\n";
		} else {
			for my $result (@result) {
				$affix_query and print $query_line, "\t";
				$prefix_output and print $identifier, "\t";
				print $result;
				$suffix_output and print "\t", $identifier;
				$append_query and print "\t", $query_line;
				print "\n";
			}
		}
	}
	return 1;				#must return true so that main program can use "writeResult() and next" statement
}


=head1 SYNOPSIS

 retrieve_compara.pl [arguments] <identifiers | idfiles>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        -f, --file <string>		file containing identifiers
            --intype <string>		input type
            --outtype <string>		output type
            --species <string>		species for output
            --homotype <string>		type of homology

            --prefix_output		prefix the output by some identifiers
            --suffix_output		suffix output by identifier
            --affix_query		prefix output with the query line
            --append_query		append output with the query line
        -s, --single_line_output	single line for each identifier
        -k, --keep_blank_line		output blank line for blank input line
        
            --host <string>		host URL for the database
            --user <string>		user name for the database
            --port <string>		port number for the database
            --password <string>		password for the database
            --ev			environment variable contains database info 

 Function: retrieve information from the ENSEMBL 
 database based on a given list of identifiers

 Example: 
 retrieve_compara.pl --inspecies human --species mouse ENSG00000164674
 retrieve_compara.pl -ins human -outs mouse -f genefile -homo ortholog

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--intype>

the type of input identifier. It can be 'geneid' only at this moment.

=item B<--outtype>

the output type depends on the input type. It can be 'geneid' only at 
this moment.

=item B<--file>

a file that contains the identifiers to use. When --file is 'stdin', 
it will read from the standard input.

=item B<--species>

homologous species database to retrieve data from. You can use any 
alias that can be recognized by ENSEMBL as valid species here.

=item B<--prefix_output>

add a prefix in every output line. This is useful when one input can 
generate several outputs. For example, one geneid may generate many 
rsid.

=item B<--suffix_output>

add a suffix in every output line. This is useful when one input can 
generate several outputs. For example, one geneid may generate many 
ortholog geneid.

=item B<--append_query>

append every output line by the entire input (query) line. This is useful when 
one input can generate several outputs, and when the input line contains many 
more information than the identifier.

=item B<--single_line_output>

each input identifier only generate one line of output, by concatenating 
the outputs together for each input identifier. This is useful for copying 
the output to Excel sheet, or for matching the input to output conveniently.

=item B<--keep_blank_line>

keep the blank lines in the output. Blank lines are generated when 
there is no database match, or if the input identifier itself is 
blank. This is useful for copying the output to Excel files as a new 
column that corresponds to an input column.

=item B<--host>

the host of the database. Usually it is the ENSEMBL database, but if 
you install a local copy of ensembl, you can use a different URL for 
faster database operation.

=item B<--user>

the user name for the ENSEMBL database, and default is 'anonymous'. If 
you install a local copy of ENSEMBL, you may set up a different user 
name and password.

=item B<--password>

the password for the --user to the ENSEMBL database. By default there 
is no password for user 'anonymous'.

=item B<--port>

the port to use to connect to the ENSEMBL database. If you install a 
local copy of ENSEMBL, you will probably need to change the port 
number to accormodate your local configuration of the MySQL server.

=item B<--homotype>

the type of homology can be 'ortholog' or 'paralog'. When this option 
is not set, all homologs will be retrieved.

=item B<--ev>

specify that certain environment variables contains database information. The 
user can set up the ENSEMBLDB, ENSEMBLPORT, ENSEMBLUSER, ENSEMBLPASS environment 
variables so that a different database (other than the ensembldb.ensembl.org) 
can be queried for faster access.

=back

=head1 DESCRIPTION

This program is used to retrieve information on genomic regions from 
the ENSEMBL MySQL database. Users supply a list of identifiers (in 
command line or in files), and are supplied with characteristics of 
genomic regions corresponding to these identifiers.

Some examples are explained below:

 1. retrieve_compara.pl --species mouse ENSG00000164674

The above command retrieve the mouse homolog corresponding to the 
human gene. About six mouse genes will be retrieved.

 2. retrieve_compara.pl --species mouse -homotype ortholog ENSG00000164674

The above command only retrieve orthologs and discard the paralogs in 
ENSEMBL. Only one mouse gene will be retrieved. You can then use 
programs such as retrieve_ensembl.pl to retrieve other information 
related to this gene.

=cut
                                                                                                                                                                       
