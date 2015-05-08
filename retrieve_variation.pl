#!/usr/bin/env perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-01-31 20:05:16 -0800 (Sat, 31 Jan 2009) $';

our ($verbose, $help, $man);
our ($file, $fh, $keep_blank_line, $single_line_output, $prefix_output, $suffix_output, $append_query, $affix_query, $query_line);
our ($intype, $outtype, $species, $consequence);
our ($host, $user, $password, $port);
our ($flanking_region, $upstream_region, $downstream_region, $population_name, $strain_name, $comment, $ev, $r2_threshold);

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'file|f=s'=>\$file, 'keep_blank_line'=>\$keep_blank_line, 'species=s'=>\$species,
	'single_line_output'=>\$single_line_output, 'prefix_output|pre'=>\$prefix_output, 'suffix_output'=>\$suffix_output, 'append_query'=>\$append_query,
	'affix_query'=>\$affix_query, 'intype=s'=>\$intype, 'outtype=s'=>\$outtype,
	'consequence=s'=>\$consequence, 'host=s'=>\$host, 'user=s'=>\$user, 'password=s'=>\$password, 'port=s'=>\$port,
	'flanking_region=i'=>\$flanking_region, 'upstream_region=i'=>\$upstream_region, 'downstream_region=i'=>\$downstream_region,
	'population_name=s'=>\$population_name, 'strain_name=s'=>\$strain_name, 'comment_flag=s'=>\$comment, 'ev'=>\$ev, 'r2_threshold=f'=>\$r2_threshold) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

#checking the validity of arguments
$file or @ARGV or pod2usage ("Error in argument: please either specify the --file argument or provide identifiers in command line\n");
$file and @ARGV and pod2usage ("Error in argument: please do not specify the --file argument and provide identifiers in command line simultaneously\n");
$outtype and $outtype eq 'seq' and $single_line_output and pod2usage ("Error in argument: when --outtype is 'seq', --single_line_output cannot be applied\n");
$consequence and $consequence =~ m/^(NSC|NSP|NS|SC|C|NC|IMPNC|NIMPNC|FRAMESHIFT_CODING|STOP_GAINED|STOP_LOST|NON_SYNONYMOUS_CODING|5PRIME_UTR|3PRIME_UTR|INTRONIC|UPSTREAM|DOWNSTREAM|INTERGENIC|ESSENTIAL_SPLICE_SITE|SPLICE_SITE|REGULATORY_REGION)$/i || pod2usage ("Error in argument: --consequence <$consequence> is not a valid variation type. Use --help to see all valid types");

#setting up default arguments
$| = 1;										#we may want to see immediate output if the output is redirected to a file handle
$species ||= 'human';
if ($ev) {
	#the --host, --user, --port and --password argument take precedence over the --ev argument
	$host ||= $ENV{ENSEMBLDB} || '';
	$user ||= $ENV{ENSEMBLUSER} || '';
	$port ||= $ENV{ENSEMBLPORT} || '';
	$password ||= $ENV{ENSEMBLPASS} || '';
}
$intype ||= 'rsid';
$outtype ||= 'rsid';
$host ||= 'ensembldb.ensembl.org';
$user ||= 'anonymous';
$password ||= '';
$port ||= 3306;
$r2_threshold ||= 0.8;
$outtype eq 'seq' and $consequence ||= 'CODING';				#set up a obvious default constraint so less warning message will be printed out
uc $consequence eq 'NSC' and $consequence = 'NON_SYNONYMOUS_CODING';
uc $consequence eq 'NSP' and $consequence = 'FRAMESHIFT_CODING|STOP_GAINED|STOP_LOST';
uc $consequence eq 'NS' and $consequence = 'FRAMESHIFT_CODING|STOP_GAINED|STOP_LOST|NON_SYNONYMOUS_CODING';
uc $consequence eq 'SC' and $consequence = 'SYNONYMOUS_CODING';
uc $consequence eq 'C' and $consequence = 'FRAMESHIFT_CODING|STOP_GAINED|STOP_LOST|NON_SYNONYMOUS_CODING|SYNONYMOUS_CODING';
uc $consequence eq 'NC' and $consequence = '5PRIME_UTR|3PRIME_UTR|INTRONIC|UPSTREAM|DOWNSTREAM|INTERGENIC|ESSENTIAL_SPLICE_SITE|SPLICE_SITE|REGULATORY_REGION';
uc $consequence eq 'IMPNC' and $consequence = '5PRIME_UTR|3PRIME_UTR|ESSENTIAL_SPLICE_SITE|SPLICE_SITE|REGULATORY_REGION';
uc $consequence eq 'NIMPNC' and $consequence = 'INTRONIC|UPSTREAM|DOWNSTREAM|INTERGENIC';

#loading registry from the database, so the various database adaptors can be retrieved automatically, without worrying about table names or versions.
$verbose and print STDERR "NOTICE: Loading registry from database\n";
Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host, -user => $user, -port => $port, -pass=>$password);
Bio::EnsEMBL::Registry->alias_exists ($species) or pod2usage ("Error in argument: the --species <$species> is not recognized by ensembl registry as valid species");

our ($slice_adaptor, $gene_adaptor, $transcript_adaptor, $v_adaptor, $vf_adaptor, $tv_adaptor, $pop_adaptor, $pop_gen_adaptor, $af_adaptor, $ldFeatureContainerAdaptor) = initializeDatabaseHandler ($species);

#preparing input/output operation to read identifiers
if (defined $file) {
	if ($file eq 'stdin') {
		$fh = \*STDIN;
	} else {
		open ($fh, $file) or confess "Error: cannot read from file $file: $!";
	}
}

#process identifiers one-by-one
my $num_processed_identifiers;
while (defined (my $identifier = nextIdentifier ($fh))) {
	$identifier eq '' and writeResult ($identifier) and next;
	$comment and $identifier =~ m/^$comment/ and print "$identifier\n" and next;
	$verbose and print STDERR "NOTICE: Processing identifier <$identifier>\n";
	if ($intype eq 'geneid') {
		my $gene = $gene_adaptor->fetch_by_stable_id ($identifier);
		$gene or print STDERR "WARNING: No gene found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, gene2output ($gene, $identifier));
	} elsif ($intype =~ m/^(transcriptid|exonid|chr_region)$/) {
		my $slice;
		if ($intype eq 'transcriptid') {
			$slice = $slice_adaptor->fetch_by_transcript_stable_id ($identifier) or print STDERR "WARNING: Unable to retrieve slice for transcript <$identifier>\n" and writeResult ($identifier) and next;
		} elsif ($intype eq 'exonid') {
			$slice = $slice_adaptor->fetch_by_exon_stable_id ($identifier) or print STDERR "WARNING: Unable to retrieve slice for exon <$identifier>\n" and writeResult ($identifier) and next;
		} elsif ($intype eq 'chr_region') {
			my ($chr, $chr_start, $chr_end);
			$identifier =~ m/^(?:chr)?(\w+):(\d+)\-(\d+)$/ or print STDERR "WARNING: Unable to recognize chromosome location <$identifier>\n" and writeResult ($identifier) and next;
			$slice = $slice_adaptor->fetch_by_region('chromosome', $1, $2, $3) or print STDERR "WARNING: Unable to retrieve slice for chromosome region <$identifier>\n" and writeResult ($identifier) and next;;
		}
		$flanking_region and $slice = $slice->expand ($flanking_region, $flanking_region);
		$upstream_region and $slice = $slice->expand ($upstream_region);
		$downstream_region and $slice = $slice->expand (0, $downstream_region);
		writeResult ($identifier, slice2output ($slice, $identifier));
	} elsif ($intype eq 'rsid') {
		my $var = $v_adaptor->fetch_by_name($identifier) or print STDERR "WARNING: Unable to retrieve variation <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, var2output ($var, $identifier));
	} elsif ($intype eq 'population_name') {
		#population table contains (sample_id and is_strain), while sample table contains (sample_id, name, size and other information)
		my $pop = $pop_adaptor->fetch_by_name ($identifier);
		$pop or print STDERR "WARNING: Unable to retrieve population $identifier\n" and writeResult ($identifier) and next;
		writeResult ($identifier, population2output ($pop, $identifier));
	} else {
		pod2usage ("Error in argument: --intype must be one of geneid, chr_region, rsid");
	}
	
	$num_processed_identifiers++;
	#re-initialize database handlers (after querying remote database for a certain period of time, there may be errors such as "DBD::mysql::st execute failed: Lost connection to MySQL server during query"
	if ($num_processed_identifiers =~ m/0000$/) {
		print STDERR "NOTICE: Re-initializing remote database handlers\n";
		($slice_adaptor, $gene_adaptor, $transcript_adaptor, $v_adaptor, $vf_adaptor, $tv_adaptor, $pop_adaptor, $pop_gen_adaptor, $af_adaptor, $ldFeatureContainerAdaptor) = initializeDatabaseHandler ($species);
	}
	#after about 38725 rsid->chr_region queries, the DB usually fails by losing connection. Therefore, I decide to re-initialize database every 10000 queries
}

sub initializeDatabaseHandler {
	my ($species) = @_;
	
	#first disconnect from all databases if already connected
	Bio::EnsEMBL::Registry->disconnect_all ();
	
	#preparing database adaptors and check whether they function correctly
	Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host, -user => $user, -port => $port, -pass=>$password);
	my $db1 = Bio::EnsEMBL::Registry->get_DBAdaptor ($species, 'core');
	my $db2 = Bio::EnsEMBL::Registry->get_DBAdaptor ($species, 'variation');
	$db1 or confess "ERROR: Unable to get the core database adaptor for species $species\n";
	$db2 or confess "ERROR: Unable to get the variation database adaptor for species $species\n";
	
	my $slice_adaptor = $db1->get_SliceAdaptor ();
	my $gene_adaptor = $db1->get_GeneAdaptor ();
	my $transcript_adaptor = $db1->get_TranscriptAdaptor ();
	
	my $v_adaptor = $db2->get_VariationAdaptor;
	my $vf_adaptor = $db2->get_VariationFeatureAdaptor ();
	my $tv_adaptor = $db2->get_TranscriptVariationAdaptor ();
	my $pop_adaptor = $db2->get_PopulationAdaptor ();
	my $pop_gen_adaptor = $db2->get_PopulationGenotypeAdaptor ();
	my $af_adaptor = $db2->get_AlleleFeatureAdaptor ();
	my $ldFeatureContainerAdaptor;
	eval {
		require IPC::Run;
		$ldFeatureContainerAdaptor = $db2->get_LDFeatureContainerAdaptor;
	};
	$@ and print STDERR "WARNING: Database handler ldFeatureContainerAdaptor cannot be loaded, some components of the program (those involving LD calculation) may not work\n";
	return ($slice_adaptor, $gene_adaptor, $transcript_adaptor, $v_adaptor, $vf_adaptor, $tv_adaptor, $pop_adaptor, $pop_gen_adaptor, $af_adaptor, $ldFeatureContainerAdaptor);
}

sub population2output {
	my ($pop, $identifier) = @_;
	my @output;
	if ($outtype eq 'info') {
		push @output, join (',', $pop->name, $pop->size||'', $pop->description||'', $pop->is_strain);
	}
	return @output;
}

#this subroutine generate outputs for a given variation object.
sub var2output {
	my ($var, $identifier) = @_;
	my @output;
	if ($outtype eq 'transcriptid') {
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		my $tvs = $tv_adaptor->fetch_all_by_VariationFeatures ($vfs);
		@$tvs or print STDERR "WARNING: No transcript variation is associated with variation <$identifier>\n" and return;
		for my $tv (@$tvs) {
			$consequence and join (',', @{$tv->consequence_type}) =~ m/^($consequence)$/i || next;
			push @output, $tv->transcript->stable_id;
		}
		my %output = map {$_, 1} @output;			#the same transcriptid might be associated with multipel variation features
		@output = sort keys %output;
	} elsif ($outtype eq 'geneid') {
		
		#I am a little reluctant to implement this --outtype option, because it may cause a lot of confusion. However, many times it is convenient to have an option like this
		#just remember, the 'geneid' here only corresponds to transcriptid that can be associated with the rsid
		#in other word, even if the rsid is located inside the gene, it may not be retrieved if it is located in an intron (non-transcript region)
		
		#the following is a verbatim copy of the above paragraph on 'transcriptid'
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		my $tvs = $tv_adaptor->fetch_all_by_VariationFeatures ($vfs);
		@$tvs or print STDERR "WARNING: No transcript variation is associated with variation <$identifier>\n" and return;
		for my $tv (@$tvs) {
			$consequence and join (',', @{$tv->consequence_type}) =~ m/^($consequence)$/i || next;
			push @output, $tv->transcript->stable_id;
		}
		my %output = map {$_, 1} @output;			#the same transcriptid might be associated with multipel variation features
		@output = sort keys %output;
		
		#the following is used to get the geneid from the list of transcriptid
		%output = ();
		for my $transcriptid (@output) {
			my $gene = $gene_adaptor->fetch_by_transcript_stable_id ($transcriptid) or print STDERR "WARNING: No gene associated with transcript id $transcriptid\n" and next;
			$output {$gene->stable_id} = 1;
		}
		@output = sort keys %output;	
	} elsif ($outtype eq 'seq') {
		#to get the sequence, one has to first get vf via vf_adaptor, then get tvf via tvf_adaptor
		#the ENSEMBL API documentation mentions that tra can directly retrieve var by fetch_all_by_Variation, which does not exist
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		
		my ($max_vf, $max_tv, $max_length);
		for my $vf (@$vfs) {
			my $tvs = $tv_adaptor->fetch_all_by_VariationFeatures ([$vf]);
		
			@$tvs or print STDERR "WARNING: No transcript is associated with variation <$identifier>\n" and next;
			$verbose and @$tvs == 1 || print STDERR "NOTICE: variation <identifier> has ${\(scalar @$tvs)} transcript variations\n";
			for my $tv (@$tvs) {
				$consequence and join (',', @{$tv->consequence_type}) =~ m/^($consequence)$/i || next;
				$max_vf ||= $vf;
				$max_tv ||= $tv;
				$max_length ||= $max_tv->transcript->length;

				if ($tv->transcript->length > $max_length) {
					$max_vf = $vf;
					$max_tv = $tv;
					$max_length = $max_tv->transcript->length;
					$verbose and print STDERR "NOTICE: A longer transcript ", $max_tv->transcript->stable_id, " found for variation ", $max_vf->variation_name, "\n";
				}
			}
		}
		defined $max_tv or print STDERR "NOTICE: variation <$identifier> does not have transcript associated with it\n" and return;
		my $fasta_head = '>' . $max_tv->transcript->stable_id . ' ' . $max_vf->variation_name . '|' . $max_tv->pep_allele_string . '|' . $max_tv->translation_start . '-' . $max_tv->translation_end;
		my $fasta_body = $max_tv->transcript->translation->seq;
		$fasta_body =~ s/(.{80})/$1\n/g;
		@output = ("$fasta_head\n$fasta_body");
	} elsif ($outtype eq 'test') {
		foreach my $vf (@{$vf_adaptor->fetch_all_by_Variation($var)}){
			print $vf->variation_name(),","; # print rsID
			print $vf->allele_string(),","; # print alleles
			print $vf->get_consequence_type(),","; # print consequenceType
			print substr($var->five_prime_flanking_seq,-10) , "[",$vf->allele_string,"]"; #print the allele string
			print substr($var->three_prime_flanking_seq,0,10), ","; # print RefSeq
			print $vf->seq_region_name, ":", $vf->start,"-",$vf->end, "\n"; # print position in Ref in format Chr:start-end
			#&get_TranscriptVariations($vf); # get Transcript information
		}
	} elsif ($outtype eq 'rsid') {
		#if the only required output is rsid without other constraints, we can simply print out the name of the variation.
		#but if user specified a special consequence type, we have to check whether at least one of the variation_feature for the variation meets the criteria
		if ($consequence) {
			for my $vf (@{$vf_adaptor->fetch_all_by_Variation ($var)}){
				join (',', @{$vf->get_consequence_type}) =~ m/\b($consequence)\b/i and push @output, $var->name and last;
			}
		} else {
			push @output, $var->name;
		}
	} elsif ($outtype eq 'chr_region') {
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		for my $vf (@$vfs) {
			$consequence and join (',', @{$vf->get_consequence_type}) =~ m/\b$consequence\b/i || next;
			push @output, $vf->seq_region_name . ':' . $vf->seq_region_start . '-' . $vf->seq_region_end;
		}
		
	} elsif ($outtype eq 'tagged_population') {
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		for my $vf (@$vfs) {
			#we can use either pop_adaptor or $vf->is_tagged here. (tagging means r2 > 0.99 in a population)
			my $pops = $pop_adaptor->fetch_tagged_Population ($vf);
			@$pops or print STDERR "WARNING: variation feature ", $vf->name, " has not been tagged in any population\n" and next;
			push @output, map {$_->name} @$pops;
		}
	} elsif ($outtype eq 'ldrsid') {
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		defined $population_name or pod2usage ("Error in argument: you must specify a population by --population_name argument for the --outtype of af (you can use --outtype of allaf if you do not want to specify a particular population)");
		my $population = $pop_adaptor->fetch_by_name ($population_name);
		for my $vf (@$vfs) {
			my $ldFeatureContainer = $ldFeatureContainerAdaptor->fetch_by_VariationFeature ($vf, $population);

			foreach my $r_square (@{$ldFeatureContainer->get_all_r_square_values}){
				if ($r_square->{r2} > $r2_threshold) {
					push @output, (($r_square->{variation1}->variation_name eq $var->name) ? $r_square->{variation2}->variation_name : $r_square->{variation1}->variation_name);
				}
			}

		}
	} elsif ($outtype eq 'ldinfo') {
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		defined $population_name or pod2usage ("Error in argument: you must specify a population by --population_name argument for the --outtype of af (you can use --outtype of allaf if you do not want to specify a particular population)");
		my $population = $pop_adaptor->fetch_by_name ($population_name);
		for my $vf (@$vfs) {
			my $ldFeatureContainer = $ldFeatureContainerAdaptor->fetch_by_VariationFeature ($vf, $population);

			foreach my $r_square (@{$ldFeatureContainer->get_all_r_square_values}){
				if ($r_square->{r2} > $r2_threshold) {
					push @output, join (',', (($r_square->{variation1}->variation_name eq $var->name) ? $r_square->{variation2}->variation_name : $r_square->{variation1}->variation_name), $r_square->{r2});
				}
			}

		}
	} elsif ($outtype eq 'genotype_freq') {
		my $pop_gens = $pop_gen_adaptor->fetch_all_by_Variation ($var);
		@$pop_gens or print STDERR "WARNING: Unable to retrieve population genotype for <$identifier>\n" and return;
		for my $pop_gen (@$pop_gens) {
			if (defined $population_name) {
				$pop_gen->population->name eq $population_name || next;
				push @output, $pop_gen->allele1 . $pop_gen->allele2 . ',' . $pop_gen->frequency;
			} else {
				push @output, $pop_gen->population->name . ',' . $pop_gen->allele1 . $pop_gen->allele2 . ',' . $pop_gen->frequency;
			}
		}
	} elsif ($outtype eq 'af') {
		my $alleles = $var->get_all_Alleles;
		defined $population_name or pod2usage ("Error in argument: you must specify a population by --population_name argument for the --outtype of af (you can use --outtype of allaf if you do not want to specify a particular population)");
		@$alleles or print STDERR "WARNING: Unable to retrieve allele for variation ", $var->name, "\n" and return;
		for my $allele (@$alleles) {
			defined $allele->population or next;
			$population_name eq $allele->population->name or next;
			push @output, join (',', $allele->allele, $allele->frequency||'');
		}
	} elsif ($outtype eq 'maf') {
		my $alleles = $var->get_all_Alleles;
		defined $population_name or pod2usage ("Error in argument: you must specify a population by --population_name argument for the --outtype of af (you can use --outtype of allaf if you do not want to specify a particular population)");
		@$alleles or print STDERR "WARNING: Unable to retrieve allele for variation ", $var->name, "\n" and return;
		my @maf;		#an array of allele frequencies
		for my $allele (@$alleles) {
			defined $allele->population or next;
			$population_name eq $allele->population->name or next;
			defined $allele->frequency or next;
			push @maf, $allele->frequency;
		}
		if (@maf) {
			@maf = sort {$a <=> $b} @maf;
			$maf[0] == 1 and unshift @maf, 0;
			push @output, $maf[0];
		}
	} elsif ($outtype eq 'allaf') {
		my $alleles = $var->get_all_Alleles;
		@$alleles or print STDERR "WARNING: Unable to retrieve allele for variation ", $var->name, "\n" and return;
		for my $allele (@$alleles) {
			defined $allele->population or next;
			push @output, join (',', $allele->population->name, $allele->allele, $allele->frequency||'');
		}
	} elsif ($outtype eq 'allmaf') {
		my $alleles = $var->get_all_Alleles;
		@$alleles or print STDERR "WARNING: Unable to retrieve allele for variation ", $var->name, "\n" and return;
		my %maf;		#key is population, value is array of allele frequencies
		for my $allele (@$alleles) {
			defined $allele->population or next;
			defined $allele->frequency or next;
			push @{$maf{$allele->population->name}}, $allele->frequency;
		}
		if (%maf) {
			for my $key (keys %maf) {
				my @maf = @{$maf{$key}};
				@maf = sort {$a <=> $b} @maf;
				$maf[0] == 1 and unshift @maf, 0;
				push @output, join (',', $key, $maf[0]);
			}
		}
	} elsif ($outtype eq 'consequence') {
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		for my $vf (@$vfs) {
			push @output, $vf->display_consequence;
		}
	} elsif ($outtype eq 'genotype') {
		#this option makes no much sense for mouse, because all genotyped mouse are actually homozygous! the population for mouse is a mixture of different strains!
		my $pop_gens = $pop_gen_adaptor->fetch_all_by_Variation ($var);
		@$pop_gens or print STDERR "WARNING: Unable to retrieve pop_gen for variation, ", $var->name, "\n" and return;
		for my $pop_gen (@$pop_gens) {
			$population_name and $population_name ne $pop_gen->population->name and next;
			push @output, join (',', $pop_gen->population->name, $pop_gen->variation->name, $pop_gen->frequency, $pop_gen->allele1 . $pop_gen->allele2);
		}
	} elsif ($outtype eq 'allele') {
		my $alleles = $var->get_all_Alleles;
		@$alleles or print STDERR "WARNING: Unable to retrieve allele for variation ", $var->name, "\n" and return;
		for my $allele (@$alleles) {
			defined $allele->population or print STDERR "WARNING: $identifier allele cannot be associated with a population\n" and next;
			$population_name and $population_name ne $allele->population->name and next;
			push @output, join (',', $allele->population->name, $allele->allele, $allele->frequency||'');
		}
	} elsif ($outtype eq 'strainallele') {
		$strain_name or pod2usage ("Error in argument: you must specify a strain name by --strain_name argument for the --outtype of strainallele");
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		for my $vf (@$vfs) {
			my @strain_name = split (/,/, $strain_name);
			if (@strain_name == 1) {
				my $strain_slice = $vf->feature_Slice->get_by_strain ($strain_name);
				push @output, $strain_slice->seq;
			} else {
				my $output = '';
				for my $nextstrain (@strain_name) {
					my $strain_slice = $vf->feature_Slice->get_by_strain ($nextstrain);
					$output .= $nextstrain . ':' . $strain_slice->seq . ',';
				}
				$output =~ s/,$//;
				push @output, $output;
			}
		}
	} elsif ($outtype eq 'straindif') {
		$strain_name or pod2usage ("Error in argument: you must specify a strain name by --strain_name argument for the --outtype of strainallele");
		my $vfs = $vf_adaptor->fetch_all_by_Variation ($var);
		@$vfs or print STDERR "WARNING: Unable to retrieve variation feature for <$identifier>\n" and return;
		$verbose and @$vfs == 1 || print STDERR "NOTICE: variation <$identifier> has ${\(scalar @$vfs)} variation features\n";
		for my $vf (@$vfs) {
			my $strain_slice = $vf->feature_Slice->get_by_strain ($strain_name);
			my $differences = $strain_slice->get_all_differences_Slice;			#get differences between reference and strain sequence
			if ($differences) {
				foreach my $diff (@$differences){
					push @output, 'DIF';
				}
			} else {
				push @output, 'SAME';
			}
		}
	} elsif ($outtype eq 'validationstatus') {
		my $val_status = $var->get_all_validation_states;
		@$val_status or print STDERR "WARNING: There is no validation status for $identifier\n" and return;
		@output = @$val_status;
	} else {
		pod2usage ("Error in argument: when --intype is rsid, --outtype can be only seq, all, rsid and others");
	}
	return @output;
}

sub gene2output {
	my ($gene, $identifier) = @_;
	my @output;
	if ($outtype eq 'rsid') {
		my $slice = $slice_adaptor->fetch_by_gene_stable_id ($gene->stable_id);
		@output = slice2output ($slice, $identifier);
	} elsif ($outtype eq 'test') {
		my $transcripts = $transcript_adaptor->fetch_all_by_Gene ($gene);
		@$transcripts or print STDERR "WARNING: No transcripts found for gene <$identifier>\n" and return;
		$verbose and print STDERR "NOTICE: For identifier <$identifier>, a total of ${\(scalar @$transcripts)} transcripts are found\n";
		my $tvs = $tv_adaptor->fetch_all_by_Transcripts ($transcripts);
		@$tvs or print STDERR "WARNING: No variation found in transcripts for gene <$identifier>\n" and return;
		$verbose and print STDERR "NOTICE: For identifier <$identifier>, a total of ${\(scalar @$tvs)} transcript variations are found\n";
	} elsif ($outtype eq 'ld' or $outtype eq 'compseq') {
		my $slice = $slice_adaptor->fetch_by_gene_stable_id ($gene->stable_id);
		@output = slice2output ($slice, $identifier);
	} else {
		pod2usage ("Error in argument: when --intype is $intype, the --outtype $outtype is not supported");
	}
	return @output;
}

#the subroutine generate the appropriate outputs for a given slice object that may contain variations.
sub slice2output {
	my ($slice, $identifier) = @_;
	my @output;
	if ($outtype eq 'rsid') {
		my $vfs = $vf_adaptor->fetch_all_by_Slice ($slice);		#return ALL variations defined in $slice
		@$vfs or print STDERR "WARNING: Unable to find variation in slice for identifier <$identifier>\n";
		foreach my $vf (@{$vfs}){
			$consequence and join (',', @{$vf->get_consequence_type}) =~ m/\b$consequence\b/i || next;
			$verbose and print STDERR "NOTICE: Variation ", $vf->variation_name, " with consequence ", join (',', @{$vf->get_consequence_type}), " and alleles ", $vf->allele_string, " in slice ", $slice->seq_region_name, " and position ", $vf->start,"-",$vf->end,"\n";
			push @output, $vf->variation_name;
		}
	} elsif ($outtype eq 'test') {
		confess;
		my $transcripts = $transcript_adaptor->fetch_all_by_Slice ($slice);
	} elsif ($outtype eq 'ld') {
		$population_name or pod2usage ("Error in argument: please specify --population_name for --outtype ld");
		my $population = $pop_adaptor->fetch_by_name($population_name);
		$population or print STDERR "WARNING: No population found with identifier <$identifier>\n" and writeResult ($identifier) and next;

		defined $ldFeatureContainerAdaptor or confess "ERROR: ldFeatureContainerAdaptor is not loaded successfully so LD calculation cannot be performed. Check whether IPC::Run is installed in the system";
		my $ldFeatureContainer = $ldFeatureContainerAdaptor->fetch_by_Slice($slice, $population);

		foreach my $r_square (@{$ldFeatureContainer->get_all_r_square_values}){
			push @output, join (',', $r_square->{r2}, $r_square->{variation1}->variation_name, $r_square->{variation2}->variation_name, $population->name, $population->size||'UNKNOWN');
		}
	} elsif ($outtype eq 'compseq') {
		$strain_name or pod2usage ("Error in argument: please specify --strain_name for --outtype of compseq");
		my $strainSlice = $slice->get_by_strain($strain_name);
		$strainSlice or print STDERR "WARNING: Unable to retrieve slice for strain $strain_name\n" and return;
		my $differences = $strainSlice->get_all_differences_Slice();
		foreach my $diff (@{$differences}){
			push @output, join (',', $diff->start, $diff->end, $diff->allele_string);
		}
	} elsif ($outtype eq 'allele') {
		my $afs;
		if ($population_name) {
			my $pop = $pop_adaptor->fetch_by_name ($population_name);
			$pop or print STDERR "WARNING: Unable to retrieve population $population_name for identifier <$identifier>\n" and return;
			$afs = $af_adaptor->fetch_all_by_Slice_Population ($slice, $pop);
		} else {
			$afs = $af_adaptor->fetch_all_IndividualSlice ($slice);
		}
		@$afs or print STDERR "WARNING: Unable to retrieve allele features from slice for identifier <$identifier>\n" and return;
		for my $af (@$afs) {
			push @output, join (',', $af->population->name, $af->variation_name, $af->allele_string);
		}
	} else {
		pod2usage ("Error in argument: for --intype $intype, valid --outtype is rsid");
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

 retrieve_varation.pl [arguments] [identifiers]

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        -i, --intype <string>		input type
        -o, --outtype <string>		output type
        -f, --file <string>		file containing identifiers
        -s, --species <string>		species for the variations
        -c, --consequence <string>	functional consequence of variation
        
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

 Function: retrieve information from the ENSEMBL-VARIATION database 
 based on a given list of identifiers such as Gene ID or dbSNP ID.

 Example:
 retrieve_variation.pl -i geneid ENSG00000182473
 retrieve_variation.pl -i geneid -o rsid -c nsc ENSG00000182473
 retrieve_variation.pl -i geneid -o seq -c nsc ENSG00000182473
 retrieve_variation.pl -intype rsid -c ns rs12944187 rs2303116 rs11881259 rs2303115
 retrieve_variation.pl -intype rsid -c ns -file snp.file
 retrieve_variation.pl -i transcriptid -o seq ENST00000335675 -c nsc -v
 retrieve_variation.pl -i geneid ENSG00000182473 -o ld -pop CSHL-HAPMAP:HapMap-CEU
 retrieve_variation.pl -i geneid ENSMUSG00000026196 -o compseq -strain 129X1/SvJ -spe mouse
 retrieve_variation.pl -i population_name -o info -spe mouse C57BL/6J
 retrieve_variation.pl -i rsid rs3022801 -spe mouse -o allele
 retrieve_variation.pl -i rsid rs3828309 -o ldinfo -pop CSHL-HAPMAP:HapMap-CEU

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--intype>

the type of input identifier can be 'geneid', 'rsid' and 'transcriptid'.

=item B<--outtype>

the output type can be 'geneid', 'rsid', 'transcriptid' and 'seq' (in 
FASTA format).

=item B<--file>

a file that contains the identifiers to use. When --file is 'stdin', 
it will read from the standard input.

=item B<--species>

species database to retrieve variations from. By default it will use 
'homo_sapeins' database. You can use any alias that can be recognized 
by ENSEMBL as valid species here.

=item B<--consequence>

the functional consequence of the variation. It could be FRAMESHIFT_CODING, 
STOP_GAINED, STOP_LOST, NON_SYNONYMOUS_CODING, SYNONYMOUS_CODING, 5PRIME_UTR, 
3PRIME_UTR, INTRONIC, UPSTREAM, DOWNSTREAM, INTERGENIC. In addition, consequence 
could be 'NSC' (non-synonymous coding), 'NSP' (non-synonymous peptide change), 
'NS' (NSC+NSP), 'SC' (synonymous coding), 'CODING' (NS+SC), 'NONAA' (non-
coding). I speculate that 'UPSTREAM' and 'DOWNSTREAM' are defined in ENSEMBL as 
within 5kbp of the gene region and confirmed this by some examples. I do not 
know how 'REGULATORY_REGION' is defined, but it is probably defined as 1kbp 
with gene (I do not really know!!!).

=item B<--prefix_output>

add a prefix in every output line. This is useful when one input can 
generate several outputs. For example, one geneid may generate many 
rsid.

=item B<--suffix_output>

add a suffix in every output line. This is useful when one input can 
generate several outputs. For example, one geneid may generate many 
rsid.

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

=item B<--ev>

specify that certain environment variables contains database information. The 
user can set up the ENSEMBLDB, ENSEMBLPORT, ENSEMBLUSER, ENSEMBLPASS environment 
variables so that a different database (other than the ensembldb.ensembl.org) 
can be queried for faster access.

=back

=head1 DESCRIPTION

This program is used to retrieve information on genomic regions from 
the ENSEMBL-VARIATION MySQL database. Users supply a list of 
identifiers (in command line or in files), and will be able to 
retrieve characteristics of corresponding genomic regions or new types 
of identifiers corresponding to the list of input identifiers. I have 
also written two other programs that can be used to retrieve 
information from the EMSEMBL-COMPARA and EMSEMBL database, 
respectively.

Currently these types of --intype are supported: geneid, transcriptid, 
rsid and chr_region. I am actively expanding the capacity of this 
program to include more input and output options.

Some specific examples for using this program are explained below:

 1. retrieve_variation.pl -i geneid ENSG00000076944
 
The above command retrieve the rsid (default --outtype) from the 
ENSEMBL-VARIATION database that are located in the GENOMIC REGION of 
the gene ENSG00000076944. Because it is totally based on genomic 
region, there is still a slight possiblity that the retrieved SNP 
belongs to another gene, if two genes share a certain genomics region.

 2. retrieve_variation.pl -i geneid -o rsid -c nsc ENSG00000076944

The above command adds a constraint '--consequence nsc' here, so that 
the output only contains rsid that is associated with non-synonymous 
coding chagnes.

 3. retrieve_variation.pl -i geneid -o seq -c nsc ENSG00000076944

The above command output the FASTA protein sequence (rather than the 
rsid as in the previous example). The sequences are TRANSCRIPT-BASED, 
which means that if the variation is associated with multiple 
transcripts, then multiple sequences will be printed out.

 4. retrieve_variation.pl -intype rsid -c ns rs4241318 rs2303116 rs11881259 rs2303115

The above command retrieves the rsid that has non-synonymous amino 
acid changes from four rsids.

 5. retrieve_variation.pl -intype rsid -c ns -file snp.file

The above command reads rsid from the snp.file file, then perform the 
same task as the previous command.

 6. retrieve_variation.pl -i transcriptid -o seq ENST00000335675 -c nsc -v

The above command prints out the FASTA sequence that corresponds to 
the transcript ENST00000335675 with verbose output. Since there might 
be multiple SNP in this transcript, the description line of the FASTA 
file will contain all the relevant information.

Additional Note:

A note on consequence type: upstream and downstream refers to 5kb from the start 
or end of transcript. splice_site refers to within 5bp from the exon-intron 
boundary (or 10 bp surrounding the boundary). I do not know what is essential 
splice site, since many SNPs (rs3926203) overlap with exon-intron boundary.

=cut                                  
