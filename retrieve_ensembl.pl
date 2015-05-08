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
our ($intype, $file, $outtype, $prefix_output, $suffix_output, $single_line_output, $keep_blank_line, $species, $host, $user, $port, $password, $nice_output, $keyword, 
	$xreftype, $xrefdb, $auto_expand, $comment, $ev, $affix_query, $append_query, $expand_max, $expand3, $expand5, $uni_direction, $biotype, $extend_slice,
	$forward_strand, $reverse_strand, $unirev_direction);
our ($fh, $query_line);

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'intype=s'=>\$intype, 'file|f=s'=>\$file, 'outtype=s'=>\$outtype,
	'prefix_output|pre'=>\$prefix_output, 'suffix_output'=>\$suffix_output, 'single_line_output'=>\$single_line_output, 'keep_blank_line'=>\$keep_blank_line,
	'species=s'=>\$species, 'host=s'=>\$host, 'user=s'=>\$user, 'port=s'=>\$port, 'password=s'=>\$password, 
	'nice_output'=>\$nice_output, 'keyword=s'=>\$keyword, 'xreftype=s'=>\$xreftype, 'xrefdb=s'=>\$xrefdb,
	'comment_flag=s'=>\$comment, 'ev'=>\$ev, 'append_query'=>\$append_query, 'affix_query'=>\$affix_query, 'extend_slice=i'=>\$extend_slice,
	'auto_expand'=>\$auto_expand, 'expand_max=i'=>\$expand_max, 'uni_direction'=>\$uni_direction, 'unirev_direction'=>\$unirev_direction, 'biotype=s'=>\$biotype,
	'expand3'=>\$expand3, 'expand5'=>\$expand5, 'forward_strand'=>\$forward_strand, 'reverse_strand'=>\$reverse_strand) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);

#checking the validity of arguments
$file or @ARGV or pod2usage ("Error in argument: please either specify the --file argument or provide identifiers in command line\n");
$file and @ARGV and pod2usage ("Error in argument: please do not specify the --file argument and provide identifiers in command line simultaneously\n");

#setting up default arguments
$intype ||= 'geneid';
$outtype ||= 'chr_region';
$species ||= 'human';
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
$xreftype ||= 'all';
$verbose and print STDERR "NOTICE: Loading registry from database\n";
Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host, -user => $user, -port => $port, -pass=>$password);
Bio::EnsEMBL::Registry->alias_exists ($species) or pod2usage ("Error in argument: the --species <$species> is not recognized by ensembl registry as valid species");
$| = 1;				#flush after each write so that we can monitor the progress of the program in real time
$expand3 and $expand5 and pod2usage ("Error in argument: you can specify either --expand3 or --expand5 but not both");

#preparing database adaptors
#although it is not recommended to use Bio::EnsEMBL::Registry->get_adaptor ($species, $group, $type) subroutine to get database adaptors, I found that the following method is about 10-30 times faster (generally the overhead is 1 second)
my $db = Bio::EnsEMBL::Registry->get_DBAdaptor ($species, 'core');
$db or confess "ERROR: Unable to get the core database adaptor for species $species\n";

my $slice_adaptor = $db->get_SliceAdaptor ();
my $gene_adaptor = $db->get_GeneAdaptor ();
my $transcript_adaptor = $db->get_TranscriptAdaptor ();
my $exon_adaptor = $db->get_ExonAdaptor ();
my $translation_adaptor = $db->get_TranslationAdaptor ();
my $marker_adaptor = $db->get_MarkerAdaptor ();						#MarkerAdaptor is stored in Bio::EnsEMBL::Map::DBSQL
#my $oligo_probe_adaptor = $db->get_OligoProbeAdaptor();
#my $oligo_array_adaptor = $db->get_OligoArrayAdaptor ();
#my $oligo_feature_adaptor = $db->get_OligoFeatureAdaptor ();				#not necessary to use them (they are not even supported in ensembl37)
my $reg_factor_adaptor = $db->get_RegulatoryFactorAdaptor ();
my $reg_feat_adaptor = $db->get_RegulatoryFeatureAdaptor ();

#preparing input/output operation to read identifiers
if (defined $file) {
	if ($file eq 'stdin') {
		$fh = \*STDIN;
	} else {
		open ($fh, $file) or confess "Error: cannot read from file $file: $!";
	}
}

#process identifiers one-by-one
while (defined (my $identifier = nextIdentifier ($fh))) {
	$identifier eq '' and writeResult ($identifier) and next;
	$comment and $identifier =~ m/^$comment/ and print "$query_line\n" and next;		#print the query line verbatim
	$verbose and print STDERR "NOTICE: Processing identifier <$identifier>\n";
	if ($intype eq 'geneid') {
		my $gene = $gene_adaptor->fetch_by_stable_id ($identifier);
		$gene or print STDERR "WARNING: No gene found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, gene2output ($gene, $identifier));
	} elsif ($intype eq 'transcriptid') {
		my $transcript = $transcript_adaptor->fetch_by_stable_id ($identifier);
		$transcript or print STDERR "WARNING: No transcript found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, transcript2output ($transcript, $identifier));
	} elsif ($intype eq 'exonid') {
		my $exon = $exon_adaptor->fetch_by_stable_id ($identifier);
		$exon or print STDERR "WARNING: No exon found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, exon2output ($exon, $identifier));
	} elsif ($intype eq 'translationid') {
		my $translation = $translation_adaptor->fetch_by_stable_id ($identifier);
		$translation or print STDERR "WARNING: No translation found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, translation2output ($translation, $identifier));		
	} elsif ($intype eq 'markerid') {
		#if multiple markers with the same name found, it contradicts the definition of marker, so we do not process such "marker"
		#presumably only one feature is associated with each marker, so we only use the first marker_feature found
		my $markers = $marker_adaptor->fetch_all_by_synonym ($identifier);
		@$markers or print STDERR "WARNING: No marker found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		@$markers == 1 or print STDERR "WARNING: Multiple markers found for identifier <$identifier>\n";
		my @output = map {marker2output ($_, $identifier)} @$markers;
		writeResult ($identifier, @output);
	} elsif ($intype =~ m/^(chr|chr_band|chr_region|cytogenetic|clone|supercontig|contig|toplevel|seqlevel)$/) {
		my $slice;
		if ($intype eq 'chr' or $intype eq 'toplevel') {								#toplevel is alias of chromosome
			$identifier =~ m/^(\d+|X|Y)$/ or print STDERR "ERROR: the supplied identifier <$identifier> is not a valid chromosome so skipeed\n" and writeResult ($identifier) and next;
			$slice = $slice_adaptor->fetch_by_region ('chromosome', $identifier);
		} elsif ($intype eq 'chr_band') {
			$identifier =~ m/^(\d+|X|Y)([pq][\d\.]+)$/ or print STDERR "ERROR: the supplied identifier <$identifier> is not a valid chromosome band so skipped\n" and writeResult ($identifier) and next;
			$slice = $slice_adaptor->fetch_by_chr_band ($1, $2);
		} elsif ($intype eq 'chr_region') {
			$identifier =~ m/(?:chr)?(\d+|X|Y):(\d+)\-(\d+)/ or print STDERR "ERROR: the supplied identifier <$identifier> is not a valid chromosome region so skipped\n" and writeResult ($identifier) and next;
			$slice = $slice_adaptor->fetch_by_region ('chromosome', $1, $2, $3);
		} elsif ($intype eq 'cytogenetic') {
			$identifier =~ m/^(\d+)([p|q]\d+(\.\d+)?)$/ or print STDERR "ERROR: the supplied identifier <$identifier> is not a valid cytogenetic identifier so skipped\n" and writeResult ($identifier) and next;
			$slice = $slice_adaptor->fetch_by_chr_band ($1, $2);
		} elsif ($intype eq 'clone' or $intype eq 'supercontig' or $intype eq 'contig' or $intype eq 'seqlevel') {	#seqlevel if alias of contig
			$slice = $slice_adaptor->fetch_by_region ($intype, $identifier);
			$slice or print STDERR "ERROR: the supplied identifier <$identifier> is not a valid clone/supercontig/contig\n" and writeResult ($identifier) and next;
			$slice = $slice->project ('chromosome')->[0]->to_Slice;
		}
		$slice or print STDERR "WARNING: No genomic region found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		if ($extend_slice) {
			$slice = $slice->expand ($extend_slice, $extend_slice);
		}
		writeResult ($identifier, slice2output ($slice, $identifier));
	} elsif ($intype eq 'xrefid') {
		writeResult ($identifier, xrefid2output ($identifier));
	} elsif ($intype eq 'regulatoryfactorid') {
		my $regfactor = $reg_factor_adaptor->fetch_by_name ($identifier);
		$regfactor or print STDERR "WARNING: No regulatory factor found for identifier <$identifier>\n" and writeResult ($identifier) and next;
		writeResult ($identifier, regfactor2output ($regfactor, $identifier));	
	} else {
		pod2usage ("Error in argument: --intype can only be geneid, transcriptid, exonid, chr, chr_band, chr_region, cytogenetic, clone, supercontig, probesetid, xrefid");
	}
}

sub regfactor2output {
	my ($regfactor, $identifier) = @_;
	my @output;
	#as of December 2006, the Ensembl database has not provided annotations for gene-regfactor correspondence in the regulatory_factor_coding table
	#therefore, this subroutine does not work at this moment.
	#the users can however manually search the cisRED database (http://www.cisred.org) or the miBase database (http://microrna.sanger.ac.uk) for the more information on regulatory factors
	if ($outtype eq 'geneid') {
		if ($regfactor->coding_gene) {
			push @output, $regfactor->coding_gene->stable_id;		#there is no such annotation yet in Ensembl database as of December 2006
		} else {
			print STDERR "WARNING: No gene found that corresponds to the regulatory factor ", $regfactor->name, "\n";
		}
	} elsif ($outtype eq 'regulatedgeneid') {
		my @regfeat = @{$reg_feat_adaptor->fetch_all_by_factor ($regfactor)};
		if (@regfeat) {
			for my $regfeat (@regfeat) {
				my @regulated_genes = @{$regfeat->regulated_genes};
				push @output, map {$_->stable_id} @regulated_genes;
			}
		}
	} elsif ($outtype eq 'type') {
		$regfactor->type and push @output, $regfactor->type;		#type means one of "miRNA", "transcription factor", "transcription factor complex"; however, most regfactor does not have "type" annotations yet
	} else {
		pod2usage ("Error in argument: the --outtype of $outtype cannot be used for the --intype $intype");
	}
	return @output;
}

sub xrefid2output {
	my ($identifier) = @_;
	my (@output, @new_common, @common);
	#how this works: Bio::EnsEMBL will call DBEntryAdaptor, which will check dbprimary_acc, display_label column in xref table, and check synonym column in external_synonym table
	if ($outtype eq 'geneid') {
		@common = @{$gene_adaptor->fetch_all_by_external_name ($identifier)};
		@common or print STDERR "WARNING: No gene found for xref id <$identifier>\n";
	} elsif ($outtype eq 'transcriptid') {
		@common = @{$transcript_adaptor->fetch_all_by_external_name ($identifier)};
		@common or print STDERR "WARNING: No transcript found for xref id <$identifier>\n";
	} elsif ($outtype eq 'translationid') {
		@common = @{$translation_adaptor->fetch_all_by_external_name ($identifier)};
		@common or print STDERR "WARNING: No translation found for xref id <$identifier>\n";
	} else {
		pod2usage ("Error in argument: the --outtype of $outtype cannot be uesd for xrefid input");
	}
	
	if ($xrefdb) {						#make sure that the external ID is indeed from the xref database that user specified
		for my $common (@common) {
			my @xrefoutput = common2xrefoutput ($common, $identifier);
			my $found;
			for my $xrefoutput (@xrefoutput) {
				$xrefoutput =~ /^$xrefdb,$identifier\b/i and $found = 1 and last;	#use i to ignore case
			}
			$found and push @new_common, $common;
		}
		@common = @new_common;
	}
	@output = map {$_->stable_id} @common;
	return @output;
}

sub slice2output {
	my ($slice, $identifier) = @_;
	my @output;
	if ($outtype eq 'geneid') {
		for my $next_gene (@{$slice->get_all_Genes}) {
			$biotype and $next_gene->biotype eq $biotype || next;
			push @output, $next_gene->stable_id;
		}

		if (not @output and $auto_expand || $expand_max) {		#do expansion when either --auto_expand or --expand_max flag is set
			my @default_expansion_list = (10_000, 100_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000);
			if ($expand_max) {
				my $original_array_dimension = scalar (@default_expansion_list);
				for my $i (0 .. $original_array_dimension-1) {
					$default_expansion_list[$original_array_dimension-1-$i] >= $expand_max and pop @default_expansion_list;
				}
				push @default_expansion_list, $expand_max;
			}

			my (@newoutput, $newslice);
			my $centerpoint = $slice->centrepoint;
			for my $expandsize (@default_expansion_list) {
				$verbose and print STDERR "NOTICE: Automatically expanding slice ${expandsize}bp to find gene\n";
				if ($expand5) {
					$newslice = $slice->expand ($expandsize, 0);
				} elsif ($expand3) {
					$newslice = $slice->expand (0, $expandsize);
				} else {
					$newslice = $slice->expand ($expandsize, $expandsize);
				}
				@newoutput = @{$newslice->get_all_Genes};
				@newoutput or next;			#make the slice size larger to make sure one slice is found
				for my $next_gene (@newoutput) {
					$biotype and $next_gene->biotype eq $biotype || next;	#check to see whether the biotype is correct if --biotype argument is specified
					if ($uni_direction) {
						if ($next_gene->seq_region_strand == 1 and $next_gene->seq_region_start > $centerpoint) {
							$reverse_strand and next;
							push @output, [$next_gene->stable_id, $next_gene->seq_region_start - $centerpoint];
						} elsif ($next_gene->seq_region_strand == -1 and $next_gene->seq_region_end < $centerpoint) {
							$forward_strand and next;
							push @output, [$next_gene->stable_id, $centerpoint - $next_gene->seq_region_end];
						}
					} elsif ($unirev_direction) {
						if ($next_gene->seq_region_strand == -1 and $next_gene->seq_region_start > $centerpoint) {
							$forward_strand and next;
							push @output, [$next_gene->stable_id, $next_gene->seq_region_start - $centerpoint];
						} elsif ($next_gene->seq_region_strand == 1 and $next_gene->seq_region_end < $centerpoint) {
							$reverse_strand and next;
							push @output, [$next_gene->stable_id, $centerpoint - $next_gene->seq_region_end];
						}
					} else {
						if ($next_gene->seq_region_start > $centerpoint) {
							$forward_strand and $next_gene->seq_region_strand == -1 and next;
							$reverse_strand and $next_gene->seq_region_strand == 1 and next;
							push @output, [$next_gene->stable_id, $next_gene->seq_region_start - $centerpoint];
						} elsif ($next_gene->seq_region_end < $centerpoint) {
							$forward_strand and $next_gene->seq_region_strand == -1 and next;
							$reverse_strand and $next_gene->seq_region_strand == 1 and next;
							push @output, [$next_gene->stable_id, $centerpoint - $next_gene->seq_region_end];
						} else {
							confess "FATAL ERROR";
						}
					}
				}
				if (@output) {
					@output = sort {$a->[1] <=> $b->[1]} @output;
					@output = ($output[0]->[0]);		#take the closest gene and discard all others
					last;
				}
			}
		}
	} elsif ($outtype eq 'groi') {
		my ($start, $end) = ($slice->start, $slice->end);
		my ($found_five, $found_three);			#found the boundaries for the five prime and three prime
		my ($newslice, @newgene);
		my ($new_five_prime_distance, $new_three_prime_distance);
		my @default_expansion_list = (10_000, 100_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000);
		if ($expand_max) {
			my $original_array_dimension = scalar (@default_expansion_list);
			for my $i (0 .. $original_array_dimension-1) {
				$default_expansion_list[$original_array_dimension-1-$i] >= $expand_max and pop @default_expansion_list;
			}
			push @default_expansion_list, $expand_max;
		}
			
		
		for my $expandsize (@default_expansion_list) {
			$verbose and print STDERR "NOTICE: Automatically expanding slice ${expandsize}bp to find gene\n";
			
			$expand5 || $expand3 and pod2usage ("Error: for GROI calculation the --expand5 and --expand3 argument are not supported at this moment!");
			$newslice = $slice->expand ($expandsize, $expandsize);
						($new_five_prime_distance, $new_three_prime_distance) = ($expandsize, $expandsize);

			@newgene = @{$newslice->get_all_Genes};
			@newgene or next;				#this slice is not large enough to find genes
			for my $next_gene (@newgene) {
				$biotype and $next_gene->biotype eq $biotype || next;	#check to see whether the biotype is correct if --biotype argument is specified
				$next_gene->seq_region_start == $start and $next_gene->seq_region_end == $end and next;
				$verbose and print STDERR "NOTICE: Next gene", $next_gene->stable_id, " start ", $next_gene->seq_region_start, " end ", $next_gene->seq_region_end, "\n";
				if ($next_gene->seq_region_start < $start) {
					$found_five = 1;
					if ($next_gene->seq_region_end < $start) {
						if ($start - $next_gene->seq_region_end < $new_five_prime_distance) {
							$new_five_prime_distance = $start - $next_gene->seq_region_end;
						}
					} else {
						$new_five_prime_distance = 1;	#there is nothing in five prime, because another gene overlap the five prime
						if ($next_gene->seq_region_end > $end) {
							$found_three = 1;
							$new_three_prime_distance = 1;	#this gene encompass the query gene
						}
					}
				} elsif ($next_gene->seq_region_start == $start) {
					$found_five = 1;
					$new_five_prime_distance = 1;
					if ($next_gene->seq_region_end >= $end) {
						$found_three = 1;
						$new_three_prime_distance = 1;
					}
				} else {
					if ($next_gene->seq_region_start <= $end) {
						if ($next_gene->seq_region_end > $end) {
							$found_three = 1;
							$new_three_prime_distance = 1;
						}
					} else {
						$found_three = 1;
						if ($next_gene->seq_region_start - $end < $new_three_prime_distance) {
							$new_three_prime_distance = $next_gene->seq_region_start - $end;
						}
					}
				}
				#$verbose and print STDERR "Next gene is ", $next_gene->stable_id, "\t", $next_gene->seq_region_start, "-", $next_gene->seq_region_end, "\n";
			}
			if ($found_five and $found_three) {
				push @output, $newslice->seq_region_name . ':' . ($start-$new_five_prime_distance+1) . '-' . ($end+$new_three_prime_distance-1);
				last;
			} elsif ($found_five) {
				$verbose and print STDERR "NOTICE: Five prime gene found, waiting to find three prime ones\n";
			} elsif ($found_three) {
				$verbose and print STDERR "NOTICE: Three prime gene found, waiting to find five prime ones\n";
			}
		}
		if (not @output) {		#if @output is still empty, check to see whether the gene is the first gene or last gene of a chromosome
			if ($start <= $default_expansion_list[$#default_expansion_list]) {
				print STDERR "WARNING: Cannot retrieve GROI for five primer upstream region for $identifier its slice contains the first gene in chromosome\n";	#example ENSMUSG00000050530
				push @output, $slice->seq_region_name . ':' . (1) . '-' . ($slice->end + $new_three_prime_distance-1);
			} else {
				my $chr_slice = $slice_adaptor->fetch_by_region ('chromosome', $slice->seq_region_name);
				if ($end + $default_expansion_list[$#default_expansion_list] >= $chr_slice->end) {
					print STDERR "WARNING: Cannot retrieve GROI for three prime downstream region for $identifier since its slice contains the last gene in chromosome\n";	#example ENSMUSG00000041936
					push @output, $slice->seq_region_name . ':' . ($start-$new_five_prime_distance+1) . '-' . ($chr_slice->end);
				} else {
					if ($found_five) {
						print STDERR "WARNING: Cannot retrieve GROI for three prime downstream region in $default_expansion_list[$#default_expansion_list] region\n";
						push @output, $slice->seq_region_name . ':' . ($start-$new_five_prime_distance+1) . '-' . ($end+$default_expansion_list[$#default_expansion_list]);
					} elsif ($found_three) {
						print STDERR "WARNING: Cannot retrieve GROI for five prime upstream region in $default_expansion_list[$#default_expansion_list] region\n";
						push @output, $slice->seq_region_name . ':' . ($start-$default_expansion_list[$#default_expansion_list]) . '-' . ($end+$new_three_prime_distance-1);
					} else {
						print STDERR "WARNING: Cannot retrieve GROI for both five prime and three prime $default_expansion_list[$#default_expansion_list] region\n";
						push @output, $slice->seq_region_name . ':' . ($start-$default_expansion_list[$#default_expansion_list]) . '-' . ($end+$default_expansion_list[$#default_expansion_list]);
					}
				}
			}
		}
	} elsif ($outtype eq 'transcriptid') {
		@output = map {$_->stable_id} @{$slice->get_all_Transcripts};
	} elsif ($outtype eq 'exonid') {
		@output = map {$_->stable_id} @{$slice->get_all_Exons};
	} elsif ($outtype eq 'translationid') {
		my @transcript = @{$slice->get_all_Transcripts};
		for my $transcript (@transcript) {
			if ($transcript->translation) {				#some transcripts do not have translation
				push @output, $transcript->translation->stable_id;
			}
		}
	} elsif ($outtype eq 'chr_region') {
		if ($nice_output) {
			push @output, $slice->name;
		} else {
			push @output, $slice->seq_region_name . ':' . $slice->start . '-' . $slice->end;
		}
	} elsif ($outtype eq 'cytogenetic') {
		my @band = sort sortCytogenetic map {$_->slice->seq_region_name . $_->display_id} @{$slice->get_all_KaryotypeBands};	#some genes such as ENSMUSG00000069049 does not have band
		if (@band > 1) {
			my ($band_chr1, $band_chr2);
			($band_chr1) = ($band[0] =~ m/^(\d+)/);
			($band_chr2) = ($band[$#band] =~ m/^(\d+)/);
			$band_chr1 == $band_chr2 and $band[$#band] =~ s/^\d+//;
			@output = ($band[0] . '-' . $band[$#band]);
		} elsif (@band) {
			@output = ($band[0]);
		} else {
			@output = ();
		}
	} elsif ($outtype eq 'size') {
		push @output, $slice->end - $slice->start + 1;
	} elsif ($outtype eq 'seq') {
		push @output, $slice->seq;
	} else {
		pod2usage ("Error in argument: when --intype is $intype, valid --outtype are geneid, transcriptid, exonid, translateid, location, externalname, cytogenetic, xref, allxref");
	}
	return @output;
}

sub gene2output {
	my ($gene, $identifier) = @_;
	my (@output);
	ref ($gene) and $gene->isa ("Bio::EnsEMBL::Gene") or confess "Error: a valid gene object must be supplied for this subroutine";
	if ($outtype eq 'geneid') {
		@output = ($gene->stable_id);
	} elsif ($outtype eq 'transcriptid') {
		my @transcript = @{$transcript_adaptor->fetch_all_by_Gene ($gene)};
		#my @transcript = @{$gene->get_all_Transcripts};					#this has the same effect as the above sentence, but both are extremely slow
		$verbose and @transcript > 1 and print STDERR "NOTICE: Sorting multiple transcripts for <$identifier> by length\n";
		@transcript = sort {$b->length <=> $a->length} @transcript;				#sort by the length of the sum of all exons
		@output = map {transcript2output ($_, $identifier)} @transcript;
	} elsif ($outtype eq 'exonid') {
		my @exon = @{$gene->get_all_Exons};
		@exon = sort {$a->seq_region_start <=> $b->seq_region_end} @exon;			#sort by start site of exons
		@output = map {exon2output ($_, $identifier)} @exon;
	} elsif ($outtype eq 'translationid' or $outtype eq 'domain' or $outtype eq 'protein_feature') {
		my @transcript = @{$gene->get_all_Transcripts};
		@transcript = sort {$b->translation->length <=> $a->translation->length} @transcript;	#sort by the length of the sum of all exons
		my @translation;
		for my $transcript (@transcript) {
			$transcript->translation and push @translation, $transcript->translation;	#some transcripts do not have translation
		}
		@translation = sort {$b->length <=> $a->length} @translation;				#sort by the length of the translations
		@output = map {translation2output ($_, $identifier)} @translation;
	} elsif ($outtype eq 'groi' or $outtype eq 'cytogenetic' or $outtype eq 'chr_region' or $outtype eq 'size') {
		my $slice = $slice_adaptor->fetch_by_gene_stable_id ($gene->stable_id);
		@output = slice2output ($slice, $identifier);
	} elsif ($outtype eq 'xref' or $outtype eq 'allxref' or $outtype eq 'displayxref') {
		@output = common2xrefoutput ($gene, $identifier);
	} elsif ($outtype eq 'regulatoryfeature') {
		my @regfeat = @{$gene->get_all_regulatory_features};
		for my $regfeat (@regfeat) {
			my $output;
			$output = 'regfeat:' . $regfeat->name;
			$output .= ',region:' . $regfeat->seq_region_name;
			$output .= ',starnd:' . $regfeat->seq_region_strand;
			$output .= ',range:' . $regfeat->seq_region_start . '-' . $regfeat->seq_region_end;
			push @output, $output;
		}
	} elsif ($outtype eq 'regulatoryfactorid') {
		my @regfeat = @{$gene->get_all_regulatory_features};
		for my $regfeat (@regfeat) {
			$regfeat->factor and $regfeat->factor->name and push @output, $regfeat->factor->name;
		}
	} elsif ($outtype eq 'biotype') {
		@output = ($gene->biotype);
	} else {
		pod2usage ("Error in argument: the --outtype $outtype cannot be used for gene object");
	}
	return @output;
}

sub transcript2output {
	my ($transcript, $identifier) = @_;
	my (@output);
	ref ($transcript) and $transcript->isa ("Bio::EnsEMBL::Transcript") or confess "Error: a valid transcript object (rather than $transcript) must be supplied for this subroutine";
	if ($outtype eq 'geneid') {
		my $gene = $gene_adaptor->fetch_by_transcript_stable_id ($transcript->stable_id);
		@output = $gene ? (gene2output ($gene, $identifier)) : ();
	} elsif ($outtype eq 'transcriptid') {
		@output = ($transcript->stable_id);
	} elsif ($outtype eq 'exonid') {
		my @exon = @{$transcript->get_all_Exons};
		@exon = sort {$a->seq_region_start <=> $b->seq_region_end} @exon;
		@output = map {exon2output ($_, $identifier)} @exon;
	} elsif ($outtype eq 'translationid' or $outtype eq 'domain' or $outtype eq 'protein_feature') {
		my $translation = $transcript->translation;
		$translation and @output = translation2output ($translation, $identifier);
	} elsif ($outtype eq 'chr_region' or $outtype eq 'cytogenetic') {
		my $slice = $slice_adaptor->fetch_by_transcript_stable_id ($transcript->stable_id);
		@output = slice2output ($slice, $identifier);
	} elsif ($outtype eq 'xref' or $outtype eq 'allxref' or $outtype eq 'displayxref') {
		@output = common2xrefoutput ($transcript, $identifier);
	} elsif ($outtype eq 'seq') {
		@output = ($transcript->spliced_seq);
	} elsif ($outtype eq 'length') {
		@output = ($transcript->length);
	} elsif ($outtype eq 'cdna_region') {
		@output = $transcript->cdna_coding_start . '-' . $transcript->cdna_coding_end;
	} elsif ($outtype eq 'coding_region') {		#coding_region is a collection of exonic regions excluding 5_UTR and 3_UTR
		my @exon = @{$transcript->get_all_Exons};
		@exon = sort {$a->seq_region_start <=> $b->seq_region_end} @exon;
		if ($transcript->seq_region_strand == 1) {
			while ($transcript->start_Exon->stable_id ne $exon[0]->stable_id) {
				print STDERR "WARNING: Exon that does not code sequence discarded\n" . $exon[0]->stable_id;
				shift @exon;
			}
			while ($transcript->end_Exon->stable_id ne $exon[$#exon]->stable_id) {
				print STDERR "WARNING: Exon that does not code sequence discarded\n" . $exon[$#exon]->stable_id;
				pop @exon;
			}
		} else {
			while ($transcript->start_Exon->stable_id ne $exon[$#exon]->stable_id) {
				print STDERR "WARNING: Exon that does not code sequence discarded\n" . $exon[$#exon]->stable_id;
				pop @exon;
			}
			while ($transcript->end_Exon->stable_id ne $exon[0]->stable_id) {
				print STDERR "WARNING: Exon that does not code sequence discarded\n" . $exon[0]->stable_id;
				shift @exon;
			}
		}
		for my $i (0 .. @exon-1) {
			push @output, $transcript->seq_region_name . ':' . ($i==0 ? $transcript->coding_region_start : $exon[$i]->seq_region_start) . '-' . ($i==$#exon ? $transcript->coding_region_end : $exon[$i]->seq_region_end);
		}
	} elsif ($outtype eq 'five_prime_utr_region') {
		my @exon = @{$transcript->get_all_Exons};
		@exon = sort {$a->seq_region_start <=> $b->seq_region_end} @exon;
		if ($transcript->seq_region_strand == 1) {
			while ($transcript->start_Exon->stable_id ne $exon[0]->stable_id) {
				push @output, $exon[0]->seq_region_name . ':' . $exon[0]->seq_region_start . '-' . $exon[0]->seq_region_end;
				shift @exon;
			}
			push @output, $exon[0]->seq_region_name . ':' . $exon[0]->seq_region_start . '-' . ($transcript->coding_region_start-1);
		} else {
			while ($transcript->start_Exon->stable_id ne $exon[$#exon]->stable_id) {
				push @output, $exon[$#exon]->seq_region_name . ':' . $exon[$#exon]->seq_region_start . '-' . $exon[$#exon]->seq_region_end;
				pop @exon;
			}
			push @output, $exon[$#exon]->seq_region_name . ':' . ($transcript->coding_region_end+1) . '-' . ($exon[$#exon]->seq_region_end);
		}
	} elsif ($outtype eq 'three_prime_utr_region') {
		my @exon = @{$transcript->get_all_Exons};
		@exon = sort {$a->seq_region_start <=> $b->seq_region_end} @exon;
		if ($transcript->seq_region_strand == 1) {
			while ($transcript->end_Exon->stable_id ne $exon[$#exon]->stable_id) {
				push @output, $exon[$#exon]->seq_region_name . ':' . $exon[$#exon]->seq_region_start . '-' . $exon[$#exon]->seq_region_end;
				pop @exon;
			}
			push @output, $exon[$#exon]->seq_region_name . ':' . ($transcript->coding_region_end+1) . '-' . ($exon[$#exon]->seq_region_end);
		} else {
			while ($transcript->end_Exon->stable_id ne $exon[0]->stable_id) {
				push @output, $exon[0]->seq_region_name . ':' . $exon[0]->seq_region_start . '-' . $exon[0]->seq_region_end;
				shift @exon;
			}
			push @output, $exon[0]->seq_region_name . ':' . $exon[0]->seq_region_start . '-' . ($transcript->coding_region_start-1);
		}
	} elsif ($outtype eq 'intron_region') {
		my @exon = @{$transcript->get_all_Exons};
		@exon = sort {$a->seq_region_start <=> $b->seq_region_end} @exon;
		for my $i (1 .. @exon-1) {
			push @output, $transcript->seq_region_name . ':' . ($exon[$i-1]->seq_region_end+1) . '-' . ($exon[$i]->seq_region_start-1);
		}
	} elsif ($outtype eq 'start_exonid') {
		@output = $transcript->start_Exon->stable_id;
	} elsif ($outtype eq 'end_exonid') {
		@output = $transcript->end_Exon->stable_id;
	} elsif ($outtype eq 'regulatoryfactorid') {
		my @regfeat = @{$reg_feat_adaptor->fetch_all_by_transcript ($transcript)};
		for my $regfeat (@regfeat) {
			$regfeat->factor and $regfeat->factor->name and push @output, $regfeat->factor->name;
		}
	} else {
		pod2usage ("Error: the --outtype $outtype cannot be used for transcript object");
	}
	return @output;
}

sub exon2output {
	my ($exon, $identifier) = @_;
	my (@output);
	ref ($exon) and $exon->isa ("Bio::EnsEMBL::Exon") or confess "Error: a valid exon object must be supplied for this subroutine";
	if ($outtype eq 'geneid') {
		my $gene = $gene_adaptor->fetch_by_exon_stable_id ($exon->stable_id);
		@output = $gene ? (gene2output ($gene, $identifier)) : ();
	} elsif ($outtype eq 'transcriptid') {
		my @transcript = @{$transcript_adaptor->fetch_all_by_exon_stable_id ($exon->stable_id)};
		@output = map {transcript2output ($_, $identifier)} @transcript;
	} elsif ($outtype eq 'exonid') {
		@output = ($exon->stable_id);
	} elsif ($outtype eq 'cytogenetic') {
		my $slice = $slice_adaptor->fetch_by_exon_stable_id ($exon->stable_id);
		@output = slice2output ($slice, $identifier);
	} elsif ($outtype eq 'chr_region') {
		@output = $exon->seq_region_name . ':' . $exon->seq_region_start . '-' . $exon->seq_region_end;
	} elsif ($outtype eq 'xref' or $outtype eq 'allxref' or $outtype eq 'displayxref') {
		@output = common2xrefoutput ($exon, $identifier);
	} else {
		pod2usage ("Error: the --outtype $outtype cannot be used for exon object");
	}
	return @output;
}

sub translation2output {
	my ($translation, $identifier) = @_;
	my (@output);
	ref ($translation) and $translation->isa ("Bio::EnsEMBL::Translation") or confess "Error: a valid translation object must be supplied for this subroutine";
	if ($outtype eq 'geneid') {
		my $gene = $gene_adaptor->fetch_by_translation_stable_id ($translation->stable_id);
		@output = $gene ? (gene2output ($gene, $identifier)) : ();
	} elsif ($outtype eq 'transcriptid' or $outtype eq 'exonid') {
		my $transcript = $transcript_adaptor->fetch_by_translation_stable_id ($translation->stable_id);
		@output = $transcript ? (transcript2output ($transcript, $identifier)) : ();
	} elsif ($outtype eq 'translationid') {
		@output = ($translation->stable_id);
	} elsif ($outtype eq 'domain') {		#domain include five categories: pfscan, scanprosite, superfamily, pfam, prints
		for my $df (@{$translation->get_all_DomainFeatures}) {
			push @output, join (',', $df->analysis->logic_name, $df->idesc, $df->start, $df->end);
		}
	} elsif ($outtype eq 'protein_feature') {	#protein is a superset of domain, but it also includes "analysis" categories such as ncoils
		for my $pf (@{$translation->get_all_ProteinFeatures}) {
			push @output, join (',', $pf->analysis->logic_name, $pf->idesc, $pf->start, $pf->end);
		}
	} else {
		pod2usage ("Error: the --outtype $outtype cannot be used for translation object");
	}
	return @output;
}

sub common2xrefoutput {
	my ($common, $identifier) = @_;
	my (@xref, @output);
	if ($outtype eq 'xref') {
		@xref = @{$common->get_all_DBEntries};
	} elsif ($outtype eq 'allxref') {
		@xref = @{$common->get_all_DBLinks};
	} elsif ($outtype eq 'displayxref') {
		@xref = ($common->display_xref);
	} else {
		@xref = @{$common->get_all_DBLinks};
	}
	if ($xreftype eq 'similarity') {
		for my $xref (@xref) {
			$xref->isa ("Bio::EnsEMBL::IdentityXref") or next;
			$xrefdb and $xrefdb eq $xref->dbname || next;			#the --xrefdb specify the xref database to be considered
			my $output = $nice_output ? join (',', $xref->db_display_name, $xref->display_id) : join (',', $xref->dbname, $xref->primary_id);
			$output .= ',' . join (',', $xref->query_identity, $xref->translation_start, $xref->translation_end, $xref->evalue||'');
			if ($xref->analysis) {
				$output .= ',' . $xref->analysis->logic_name;
			}
			push @output, $output;
		}
	} elsif ($xreftype eq 'go') {
		for my $xref (@xref) {
			$xref->isa ("Bio::EnsEMBL::GoXref") or next;
			$xrefdb and $xrefdb eq $xref->dbname || next;			#the --xrefdb specify the xref database to be considered
			my $output = $nice_output ? join (',', $xref->db_display_name, $xref->display_id) : join (',', $xref->dbname, $xref->primary_id);
			$output .= ',' . join (',', @{$xref->get_all_linkage_types});
			push @output, $output;
		}
	} elsif ($xreftype eq 'other') {
		for my $xref (@xref) {
			$xref->isa ("Bio::EnsEMBL::IdentityXref") || $xref->isa ("Bio::EnsEMBL::GoXref") and next;
			$xrefdb and $xrefdb eq $xref->dbname || next;			#the --xrefdb specify the xref database to be considered
			my $output = $nice_output ? join (',', $xref->db_display_name, $xref->display_id) : join (',', $xref->dbname, $xref->primary_id);
			push @output, $output;
		}
	} elsif ($xreftype eq 'all') {
		for my $xref (@xref) {
			defined $xref or print STDERR "WARNING: xref for $identifier not defined so skipped!\n" and next;
			if ($nice_output) {
				$xrefdb and $xrefdb eq $xref->db_display_name || next;
				push @output, join (',', $xref->db_display_name, $xref->display_id);
			} else {
				$xrefdb and $xrefdb eq $xref->dbname || next;			#the --xrefdb specify the xref database to be considered
				push @output, join (',', $xref->dbname, $xref->primary_id);
			}
		}
	} else {
		pod2usage ("Error in argument: the specified --xreftype $xreftype is not recognized. Valid values include 'similarity', 'go', 'other' and 'all' (default)");
	}
	my %output = map {$_, 1} @output;
	@output = sort keys %output;
	return @output;
}

sub marker2output {
	my ($marker, $identifier) = @_;
	my (@output);
	ref ($marker) and $marker->isa ("Bio::EnsEMBL::Map::Marker") or confess "Error: a valid marker object must be supplied for this subroutine";
	if ($outtype eq 'geneid' or $outtype eq 'transcriptid' or $outtype eq 'exonid' or $outtype eq 'translationid' or $outtype eq 'chr_region' or $outtype eq 'cytogenetic') {
		my $marker_features = $marker->get_all_MarkerFeatures;			#marker feature means that locations where this marker has been mapped to the genome via e-PCR
		$marker_features or print STDERR "WARNING: The marker <$identifier> is not found in genome by e-PCR\n" and $marker_features = [];
		for my $marker_feature (@$marker_features) {
			#if use $slice = $marker_feature->slice here, then the slice will be the entire chromosome
			my $slice = $slice_adaptor->fetch_by_Feature ($marker_feature);
			push @output, slice2output ($slice, $identifier);
		}			
	} elsif ($outtype eq 'synonym') {
		my @synonym = @{$marker->get_all_MarkerSynonyms};
		@output = map {join (',', $_->name, $_->source)} @synonym;
	} else {
		pod2usage ("Error: the --outtype $outtype cannot be used for marker object");
	}
	return @output;
}

#this subroutine is used to sort cytogenetic identifiers such as 3q14.2 and 5p7. The chromosome does not have to be identified, so a comparison of p23.2 and p12 is also allowed.
sub sortCytogenetic {
	#copy variable $a and $b to local variable for easy manipulation
	my ($element1, $element2) = ($a, $b);
	$verbose and print STDERR "NOTICE: Sorting chromosome $element1 and $element2\n";
	uc $element1 eq uc $element2 and return -1;		#when two cytogenetic locations are the same, we arbitrarily treat the first one as smaller

	$element1 =~ m/^(\d*|[A-Z]*)([a-z])(\d*)(\.(\d+))?$/ or print STDERR "WARNING: invalid cytogenetic format <$element1> found\n" and return -1;
	my ($a_chr, $a_pq, $a_secondary, $a_tertiary) = ($1, $2, $3, $5);
	$element2 =~ m/^(\d*|[A-Z]*)([a-z])(\d*)(\.(\d+))?$/ or print STDERR "WARNING: invalid cytogenetic format <$element2> found\n" and return 1;
	my ($b_chr, $b_pq, $b_secondary, $b_tertiary) = ($1, $2, $3, $5);
	
	if ($a_chr and not $b_chr or not $a_chr and $b_chr) {
		confess "Error: cytogenetic format for <$a> and <$b> does not compare: they should both have or not have the chromosome identifier";
	}
	defined $a_pq and defined $b_pq and defined $a_secondary and defined $b_secondary or confess "Error: invalid cytogenetic format for <$a> and <$b>";

	#chromosome can be numbers but can be letters, too. We arbitrarily sort letters after any numbers
	if ($a_chr and $b_chr) {
		if ($a_chr =~ m/^\d+$/ and $b_chr =~ m/^\d+$/) {
			$a_chr == $b_chr or return $a_chr <=> $b_chr;
		} else {
			$a_chr eq $b_chr or return $a_chr cmp $b_chr;
		}
		#$a_chr =~ m/[A-Z]/i and $a_chr = ord ($a_chr) + 1000;
		#$b_chr =~ m/[A-Z]/i and $b_chr = ord ($b_chr) + 1000;
		#$a_chr == $b_chr or return $a_chr <=> $b_chr;
	}
	$a_pq eq $b_pq or return $a_pq cmp $b_pq;
	$a_secondary == $b_secondary or return $a_secondary <=> $b_secondary;
	return $a_tertiary <=> $b_tertiary;

	#if two cytogenetic locations are the same
	confess "ERROR: cannot compare cytogenetic location <$element1> and <$element2>";
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

 retrieve_ensembl.pl [arguments] <identifiers | idfiles>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        -i, --intype <string>		input type
        -o, --outtype <string>		output type
        -f, --file <string>		file containing identifier
        -s, --species <string>		species for the variations

            --prefix_output		prefix output by identifier
            --suffix_output		suffix output by identifier
            --affix_query		prefix output with the query line
            --append_query		append output with the query line
        -s, --single_line_output	single line for each identifier
        -k, --keep_blank_line		output blank line for blank input line
        -n, --nice_output		use nice output for db name
            --keyword <string>		filter results by this keyword
            --comment_flag <letter>	a flag indicating a comment line in input
            --auto_expand		automatically expand slice until one result found
            --expand_max <int>		maximum size of slice expansion
            --uni_direction		expand to downstream of orginal point only
            --biotype <string>		filter results for gene by this biotype
            --extend_slice <int>	extend both ends when input is slice-type
            --expand3			expand 3 prime end
            --expand5			expand 5 prime end
            --forward_strand		consider forward strand only
            --reverse_strand		consider reverse strand only
        
            --host <string>		host URL for the database
            --user <string>		user name for the database
            --port <string>		port number for the database
            --password <string>		password for the database
            --ev			environment variable contains database info 
            
            --xreftype <string>		the type of cross-reference
            --xrefdb <string>		the name of cross-reference database

 Function: retrieve information from the ENSEMBL 
 database based on a given list of identifiers

 Example: 
 retrieve_ensembl.pl -i geneid -o exonid ENSG00000168952
 retrieve_ensembl.pl -i geneid -o displayxref -nice ENSG00000168952
 retrieve_ensembl.pl -i transcriptid -o length ENST00000367081
 retrieve_ensembl.pl -i chr_band 6q25.3 -o translationid
 retrieve_ensembl.pl -i chr_region 6:158992770-159105889
 retrieve_ensembl.pl -i markerid D5S1958
 retrieve_ensembl.pl -i geneid -o regulatoryfactorid ENSG00000168952 
 retrieve_ensembl.pl -i regulatoryfactorid -o regulatedgeneid crtHsap24718

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--intype>

the type of input identifier can be 'geneid', 'transcriptid', 
'exonid', 'chr', 'chr_band', 'chr_region', 'clone', 'supercontig', 
'externalname', 'probesetid', 'markerid'.

=item B<--outtype>

the output type depends on the input type.

=item B<--file>

a file that contains the identifiers to use. When --file is 'stdin', 
it will read from the standard input.

=item B<--species>

species database to retrieve data from. By default it will use 
'homo_sapeins' database. You can use any alias that can be recognized 
by ENSEMBL as valid species here.

=item B<--prefix_output>

add a prefix in every output line. This is useful when one input can 
generate several outputs. For example, one geneid may generate many 
transcriptid.

=item B<--suffix_output>

add a suffix in every output line. This is useful when one input can 
generate several outputs. For example, one geneid may generate many 
transcriptid.

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

=item B<--auto_expand>

automatically expand a slice until the desired --outtype is found. This is a 
highly experimental feature. I implement it only because I need to use this 
functionality in my own research, not because it is an important or commonly-
used feature.

=item B<--comment_flag>

A flag can be used to indicate that a line in the input is actually a comment 
line and should be skipped. This line will be printed verbatim in output. 
Usually you can use '#' as the comment character, but you can also use others 
such as '>' as comment character, too.


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

=item B<--nice_output>

use nice and human-readable output for database names and other display ids. For 
example, the nice output for "HUGO" is "HGNC Symbol", while the latter is much 
more intuitive for users to understand.

=item B<--keyword>

this argument is poorly defined and not used at this moment

=item B<--xreftype>

specify the cross-reference type, which can be 'similarity', 'go' and 'other'. 
When this argument is not specified, all types of cross-reference will be used. 
This argument is a conceptual superset of --xrefdb, which specifes the name of 
the individual database that provides the cross-refernece. Another side effect 
of using this argument is that type-specific information will be also included 
in the output. For example, when -xreftype is 'go', this program will also 
output the linkage type such as "IEA" and "TAS".

=item B<--xrefdb>

specify the cross-reference database name. This is a convenience argument that 
is used to filter out unnecessary results. You can always use the Unix/Linux 
command "grep" to do essentially similar things.

=item B<--ev>

specify that certain environment variables contains database information. The 
user can set up the ENSEMBLDB, ENSEMBLPORT, ENSEMBLUSER, ENSEMBLPASS environment 
variables so that a different database (other than the ensembldb.ensembl.org) 
can be queried for faster access.

=item B<--extend_slice>

when --intype indicate a slice-type 
(chr|chr_band|chr_region|cytogenetic|clone|supercontig|contig|toplevel|seqlevel)
, extend the slice by specified length at both ends.

=back

=head1 DESCRIPTION

This program is used to retrieve information on genomic regions from 
the ENSEMBL MySQL database. Users supply a list of identifiers (in 
command line or in files), and will be able to retrieve 
characteristics of corresponding genomic regions or new types of 
identifiers corresponding to the list of input identifiers. I have 
also written two other programs that can be used to retrieve 
information from the EMSEMBL-COMPARA and EMSEMBL-VARIATION database, 
respectively.

Some examples are explained below:

=over 8

=item 1. B<retrieve_ensembl.pl -i geneid -o transcriptid ENSG00000139618>

=item 2. B<retrieve_ensembl.pl -i geneid -o cytogenetic ENSG00000139618>

=item 3. B<retrieve_ensembl.pl -i geneid -o chr_region ENSG00000139618>

=item 4. B<retrieve_ensembl.pl -i geneid -o domain -file inputfile>

=item 5. B<retrieve_ensembl.pl -i geneid -o allxref ENSG00000139618 ENSG00000168952>

The above five commands all have the same --intype as 'geneid', which 
indicates the type of the identifier. The identifier can be supplied 
in command line, or supplied in a file (given by the --file argument, 
with one identifier per line in the file).

=item 6. B<retrieve_ensembl.pl -i transcriptid -o length ENST00000367081>

This command will take the transcriptid as --intype, then output the length of 
the transcript. You can try to use --outtype of 'seq' and see what will happen.

=item 7. B<retrieve_ensembl.pl -i chr X>

=item 8. B<retrieve_ensembl.pl -i chr_band 6q25.3>

=item 9. B<retrieve_ensembl.pl -i chr_region 6:158992770-159105889>

=item 10. B<retrieve_ensembl.pl -i clone AL359765.6 -o exonid>

=item 11. B<retrieve_ensembl.pl -i supercontig NT_011333 -o cytogenetic>

The above five commands all retrieve a slice from the genome, based on given 
identifier and --intype information. The default --outtype is a string 
representing the chromosome region, but you can also get the genes, transcripts 
or translations on this slice.

=item 12. retrieve_ensembl.pl -i xref -species mouse BRCA2

The above command retrieves the geneid for the given identifier in an external 
database for mouse.

=item 13. retrieve_ensembl.pl -i xrefid -o geneid 1810 -xrefdb EntrezGene

The above command retrieves genes that have an EntrezGene id '1810', then output 
the geneid.

=item 14. retrieve_ensembl.pl -i markerid D5S1958 -o synonym

The above command retrieves marker with name 'D5S1958', then print out all the 
synonyms for this marker, as well as the map name.

=item 14. retrieve_ensembl.pl -i markerid D5S1958 -o synonym -sing -pre

The above command do the same thing, but check the output format: all output are 
in single line, and prefixed with the input identifier.

=back

B<Some more technical details and caveats are described below>

1. Internally, there are two ways of outputting something: entity-based or 
slice-based. (Entity includes gene, transcript, translation, etc., while slice 
includes contig, chr, marker, etc.). These treatment are based on common sense of 
what is more appropriate and more desirable in a typical situation.

2. A summary of valid B<--outtype> for each B<--intype> is listed below:

=over 8

=item B<--intype geneid>

geneid, transcriptid, exonid, translationid, domain, protein_feature, chr_region, cytogenetic, xref, allxref, displayxref, biotype

=item B<--intype transcriptid>

geneid, transcriptid, exonid, translationid, domain, protein_feature, chr_region, cytogenetic, xref, allxref, displayxref, seq, length

=item B<--intype exonid>

geneid, transcriptid, exonid, chr_region, cytogenetic, xref, allxref, displayxref

=item B<--intype translationid>

geneid, transcriptid, exonid, translationid, domain, protein_feature

=item B<--intype markerid>

geneid, transcriptid, exonid, translationid, chr_region, cytogenetic, synonym

=item B<--intype chr|chr_band|chr_region|cytogenetic|clone|supercontig|contig|toplevel|seqlevel>

geneid, transcriptid, exonid, translationid, chr_region, cytogenetic

=item B<--intype xrefid>

geneid, transcriptid, translationid

=item B<--intype probesetid>

geneid (this intype is highly experimental and depreciated: use --intype xref instead to achieve higher coverage)

=back

=cut
                                                                                                                                                                                                                                                                                                                                                                                                                                                            
