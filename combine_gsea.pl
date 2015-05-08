#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Pod::Usage;

our $VERSION = 			'$Revision: bbb13c8a31de6a6e9a1e71ca347a7d02a855a27b $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2009-03-31 09:21:19 -0700 (Tue, 31 Mar 2009) $';

our ($verbose, $help, $man, $setmin, $setmax, $traditional, $variant, $skipfdr, $tolerate);
our @logfile;

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'setmin=i'=>\$setmin, 'setmax=i'=>\$setmax, 'traditional'=>\$traditional, 'variant'=>\$variant,
	'skipfdr'=>\$skipfdr, 'tolerate'=>\$tolerate) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);

(@logfile) = @ARGV;
my ($set_size, $set_nominalp, $set_es, $set_nes, $set_espi, $set_nespi) = ({}, {}, {}, {}, {}, {});
my ($cycle, @cycle_max, @cycle_min, @result);

for my $logfile (@logfile) {
	readLogFile ($logfile, $set_size, $set_es, $set_espi);
}

print STDERR "NOTICE: Total of ${\(scalar keys %$set_size)} gene sets will be analyzed\n";
for my $setid (keys %$set_size) {
	defined $cycle or $cycle = scalar (@{$set_espi->{$setid}}) and print STDERR "NOTICE: Combining data from $cycle permutation cycles to calculate statistical significance\n";
	$cycle == scalar (@{$set_espi->{$setid}}) or confess "Error: discordant occurance for set $setid (previous=$cycle current=" . scalar (@{$set_espi->{$setid}});

	my ($nominalp, $nes, $nespi);
	if ($traditional) {
		($nominalp, $nes, $nespi) = calculateNES_TRADITIONAL ($set_es->{$setid}, $set_espi->{$setid});
	} elsif ($variant) {
		($nominalp, $nes, $nespi) = calculateNES1 ($set_es->{$setid}, $set_espi->{$setid});
	} else {
		($nominalp, $nes, $nespi) = calculateNES1 ($set_es->{$setid}, $set_espi->{$setid});
	}
	$set_nominalp->{$setid} = $nominalp;
	$set_nes->{$setid} = $nes;
	$set_nespi->{$setid} = $nespi;
}

#calculate max and min for each cycle, for the FWER calculation
for my $i (0 .. $cycle-1) {
	my ($cycle_max, $cycle_min);				#maximum and minimum (over all gene sets) value of this cycle
	for my $setid (keys %$set_size) {
		$cycle_max ||= $set_nespi->{$setid}[$i];
		$cycle_min ||= $set_nespi->{$setid}[$i];
		$cycle_max < $set_nespi->{$setid}[$i] and $cycle_max = $set_nespi->{$setid}[$i];
		$cycle_min > $set_nespi->{$setid}->[$i] and $cycle_min = $set_nespi->{$setid}[$i];
	}
	push @cycle_max, $cycle_max;
	push @cycle_min, $cycle_min;
}

for my $setid (sort keys %$set_size) {
	if ($setid eq 'PATH303' or $setid eq 'PATH99') {	#debugging purposes
		$verbose = 1;
	} else {
		$verbose = 0;
	}
	$verbose and print STDERR "NOTICE: Processing set $setid ... ";

	my ($fdr, $fwer);
	if ($skipfdr) {
		($fdr, $fwer) = ('NA', 'NA');
	} elsif ($traditional) {
		$fdr = calculateFDR_TRADITIONAL ($set_nes->{$setid}, $set_size, $set_nes, $set_nespi);
		$fwer = calculateFWER_TRADITIONAL ($set_nes->{$setid}, \@cycle_max, \@cycle_min);
	} elsif ($variant) {
		$fdr = calculateFDR2 ($set_nes->{$setid}, $set_size, $set_nes, $set_nespi);
		$fwer = calculateFWER1 ($set_nes->{$setid}, \@cycle_max, \@cycle_min);
	} else {
		$fdr = calculateFDR1 ($set_nes->{$setid}, $set_size, $set_nes, $set_nespi);
		$fwer = calculateFWER1 ($set_nes->{$setid}, \@cycle_max, \@cycle_min);
	}
	#print STDERR "Done with FDR=$fdr FWER=$fwer\n";
	push @result, [$setid, $set_size->{$setid}, $set_es->{$setid}, $set_nes->{$setid}, $set_nominalp->{$setid}, $fdr, $fwer];
}

outputResult (@result);


sub readLogFile {
	my ($logfile, $set_size, $set_es, $set_espi) = @_;
	if ($logfile eq 'stdin') {
		*LOG = *STDIN;
	} else {
		open (LOG, $logfile) or confess "Error: cannot read from logfile $logfile: $!";
	}
	while (<LOG>) {
		s/[\r\n]+$//;
		my ($setid, $size, $es, $nes, $nominalp, @espi) = split (/\t/, $_);
		@espi or confess "Error: invalid record found in logfile $logfile: at least 6 tab-delimited fields expected: <$_>";
		$size =~ s/size=//;
		$es =~ s/ES=//;
		$setmin and $size < $setmin and next;
		$setmax and $size > $setmax and next;
		
		if ($set_size->{$setid}) {
			$size eq $set_size->{$setid} or confess "Error: setsize discordance between previous=$set_size->{$setid} and current=$size";
			if ($es ne $set_es->{$setid}) {
				if ($tolerate) {
					print STDERR "WARNING: setES discordance between previous=$set_es->{$setid} and current=$es\n";
				} else {
					confess "Error: setES discordance between previous=$set_es->{$setid} and current=$es";
				}
			}
		} else {
			$set_size->{$setid} = $size;
			$set_es->{$setid} = $es;
		}
		push @{$set_espi->{$setid}}, @espi;
	}
	close (LOG);
}

sub calculateNES_TRADITIONAL {
	my ($es, $espi) = @_;
	my ($nominalp, $nes, @nespi) = (0);
	my (@pos_espi, @neg_espi, $mean_posespi, $mean_negespi);
	if ($es >= 0) {
		@pos_espi = grep {$_>=0} @$espi;
		@neg_espi = grep {$_<0} @$espi;
		@pos_espi or push (@pos_espi, 0.1) and print STDERR "FATAL ERROR: unable to find positive mean (increase --cycle argument might help) but 0.1 is used here to prevent program exit\n";
		$mean_posespi = mean (\@pos_espi);
		@neg_espi and $mean_negespi = mean (\@neg_espi);
		$nes = $es / $mean_posespi;

		for (@$espi) {
			if ($_ >= 0) {
				push @nespi, $_/$mean_posespi;
			} else {
				push @nespi, -$_/$mean_negespi;
			}
		}
		for (@pos_espi) {
			$_ >= $es and $nominalp++;
		}
		$nominalp /= @pos_espi;
	} else {
		@pos_espi = grep {$_>=0} @$espi;
		@neg_espi = grep {$_<0} @$espi;
		@neg_espi or push (@neg_espi, -0.1) and print STDERR "FATAL ERROR: unable to find negative mean (increase --cycle argument might help) but -0.1 is used here to prevent program exit\n";
		$mean_negespi = mean (\@neg_espi);
		@pos_espi and $mean_posespi = mean (\@pos_espi);
		$nes = -$es / $mean_negespi;

		for (@$espi) {
			if ($_ >= 0) {
				push @nespi, $_/$mean_posespi;
			} else {
				push @nespi, -$_/$mean_negespi;
			}
		}

		for (@neg_espi) {
			$_ <= $es and $nominalp++;
		}
		$nominalp /= @neg_espi;
	}
	return ($nominalp, $nes, \@nespi);
}

sub calculateNES1 {
	my ($es, $espi) = @_;
	my ($nominalp, $nes, @nespi) = (0);
	my $mean_espi = mean ($espi);
	my $sd_espi = sd ($espi);
	$nes = ($es-$mean_espi)/$sd_espi;
	for (@$espi) {
		push @nespi, ($_-$mean_espi)/$sd_espi;
		$_>=$es and $nominalp++;
	}
	$nominalp /= @$espi;
	return ($nominalp, $nes, \@nespi);
}

sub calculateFDR1 {
	my ($current_nes, $geneset, $set_nes, $set_nespi) = @_;
	my ($count_nespi, $count_nespi_flag, $count_nes_flag, $count_num_geneset, $fdr) = (0, 0, 0, 0);
	for my $nextsetid (keys %$geneset) {						#retrieve information for all these gene sets
		for (@{$set_nespi->{$nextsetid}}) {
			$_ >= $current_nes and $count_nespi_flag++;
			$count_nespi++;
		}
	}

	$count_nes_flag = scalar (grep {$_ >= $current_nes} values %$set_nes);
	$count_num_geneset = scalar (values %$set_nes);

	$fdr = $count_nespi_flag / $count_nespi / ($count_nes_flag / $count_num_geneset);
	$verbose and print STDERR "NOTICE: FDR calculation details: $count_nespi_flag / $count_nespi / $count_nes_flag / $count_num_geneset\n";
	if ($fdr >= 1) {
		#printf ("WARNING: HUGE FDR: %i %i %i %i $fdr\n", $count_nespi_flag, $count_nespi, $count_nes_flag, $count_num_geneset);
		$fdr = 1;
	}
	return $fdr;
}

sub calculateFDR2 {
	my ($current_nes, $geneset, $set_nes, $set_nespi) = @_;
	my ($count_nespi, $count_nespi_flag, $count_nes_flag, $count_num_geneset, $fdr) = (0, 0, 0, 0);
	for my $nextsetid (keys %$geneset) {						#retrieve information for all these gene sets
		my $temp = 0;
		
		#$set_nes->{$nextsetid} >= $current_nes or next;				#this is the difference!!!
		
		for (@{$set_nespi->{$nextsetid}}) {
			$_ >= $current_nes and $count_nespi_flag++;
			$count_nespi++;
			$_ >= $current_nes and $temp++;
		}
		$current_nes > 2 and $temp > 10 and print STDERR "NOTICE: Bad set=$nextsetid count=$temp\n";
	}

	$count_nes_flag = scalar (grep {$_ >= $current_nes} values %$set_nes);
	$count_num_geneset = scalar (values %$set_nes);

	$fdr = $count_nespi_flag / $count_nespi / ($count_nes_flag / $count_num_geneset);
	$current_nes > 2 and print STDERR "NOTICE: FDR calculation details: $count_nespi_flag / $count_nespi / $count_nes_flag / $count_num_geneset\n";
	if ($fdr >= 1) {
		#printf ("WARNING: HUGE FDR: %i %i %i %i $fdr\n", $count_nespi_flag, $count_nespi, $count_nes_flag, $count_num_geneset);
		$fdr = 1;
	}
	return $fdr;
}

sub calculateFDR3 {
	my ($current_nes, $geneset, $set_nes, $set_nespi) = @_;
	my ($count_nespi, $count_nespi_flag, $count_nes_flag, $count_num_geneset, $fdr) = (0, 0, 0, 0);
	for my $nextsetid (keys %$geneset) {						#retrieve information for all these gene sets
		for (@{$set_nespi->{$nextsetid}}) {
			if ($current_nes >= 0) {					#process positive side
				$_ >= $current_nes and $count_nespi_flag++;
				#$_ >= 0 and $count_nespi++;				#there is no need to dichotomize the count_nespi (2008Oct)
				$count_nespi++;
			} else {
				$_ <= $current_nes and $count_nespi_flag++;
				#$_ < 0 and $count_nespi++;				#there is no need to dichotomize the count_nespi (2008Oct)
				$count_nespi++;
			}
		}
	}
	if ($current_nes >= 0) {
		$count_nes_flag = scalar (grep {$_ >= $current_nes} values %$set_nes);
		$count_num_geneset = scalar (grep {$_>=0} values %$set_nes);
	} else {
		$count_nes_flag = scalar (grep {$_ <= $current_nes} values %$set_nes);
		$count_num_geneset = scalar (grep {$_<0} values %$set_nes);
	}
	if (!$count_nespi or !$count_nes_flag or !$count_num_geneset) {
		print STDERR "WARNING: FDR cannot be calculated so 1 is assumed\n";
		$fdr = 1;
	} else {
		$fdr = $count_nespi_flag / $count_nespi / ($count_nes_flag / $count_num_geneset);
		if ($fdr >= 1) {
			#printf ("WARNING: HUGE FDR: %i %i %i %i $fdr\n", $count_nespi_flag, $count_nespi, $count_nes_flag, $count_num_geneset);
			$fdr = 1;
		}
	}
	return $fdr;
}

#traditional ways to calculate FDR that dichotimize the positive and negative values
sub calculateFDR_TRADITIONAL {
	my ($current_nes, $geneset, $set_nes, $set_nespi) = @_;
	my ($count_nespi, $count_nespi_flag, $count_nes_flag, $count_num_geneset, $fdr) = (0, 0, 0, 0);
	for my $nextsetid (keys %$geneset) {						#retrieve information for all these gene sets
		for (@{$set_nespi->{$nextsetid}}) {
			if ($current_nes >= 0) {					#process positive side
				$_ >= $current_nes and $count_nespi_flag++;
				#$_ >= 0 and $count_nespi++;				#there is no need to dichotomize the count_nespi (2008Oct)
				$count_nespi++;
			} else {
				$_ <= $current_nes and $count_nespi_flag++;
				#$_ < 0 and $count_nespi++;				#there is no need to dichotomize the count_nespi (2008Oct)
				$count_nespi++;
			}
		}
	}
	if ($current_nes >= 0) {
		$count_nes_flag = scalar (grep {$_ >= $current_nes} values %$set_nes);
		$count_num_geneset = scalar (grep {$_>=0} values %$set_nes);
	} else {
		$count_nes_flag = scalar (grep {$_ <= $current_nes} values %$set_nes);
		$count_num_geneset = scalar (grep {$_<0} values %$set_nes);
	}
	if (!$count_nespi or !$count_nes_flag or !$count_num_geneset) {
		print STDERR "WARNING: FDR cannot be calculated so 1 is assumed\n";
		$fdr = 1;
	} else {
		$fdr = $count_nespi_flag / $count_nespi / ($count_nes_flag / $count_num_geneset);
		if ($fdr >= 1) {
			printf ("WARNING: HUGE FDR: %i %i %i %i $fdr\n", $count_nespi_flag, $count_nespi, $count_nes_flag, $count_num_geneset);
			$fdr = 1;
		}
	}
	return $fdr;
}

sub calculateFWER1 {
	my ($current_nes, $cycle_max, $cycle_min) = @_;
	my (@more_extreme, $fwer);
	@more_extreme = grep {$_ >= $current_nes} @$cycle_max;
	$fwer = @more_extreme/@$cycle_max;
	return $fwer;
}

#traditional FWER calculate dichotomize the NES values so positive and negative values have different method
sub calculateFWER_TRADITIONAL {
	my ($current_nes, $cycle_max, $cycle_min) = @_;
	my (@more_extreme, $fwer);
	my $cycle_pos_max = scalar (grep {$_ >= 0} @$cycle_max);
	my $cycle_neg_min = scalar (grep {$_ <0 } @$cycle_min);
	if ($current_nes >= 0) {
		@more_extreme = grep {$_ >= $current_nes} @$cycle_max;
		$fwer = @more_extreme/$cycle_pos_max;
	} else {
		@more_extreme = grep {$_ <= $current_nes} @$cycle_min;
		$fwer = @more_extreme/$cycle_neg_min;
	}
	return $fwer;
}

sub outputResult {
	my (@result) = @_;
	@result = sort {$b->[3] <=> $a->[3]} @result;
	if ($traditional) {
		print "<-----------------Over-represented in head of ranked list-------------------------->\n";
		for my $result (@result) {
			$result->[3] >= 0 or next;
			print "Geneset=$result->[0]\tSize=", sprintf ("%-4d", $result->[1]), "\tES=", sprintf ("%.3f", $result->[2]), "\tNES=", sprintf ("%.3f", $result->[3]), "\tNominalP=", sprintf ("%.5f", $result->[4]), "\tFDR=", ($result->[5] eq 'NA')?'NA':sprintf ("%.3f", $result->[5]), "\tFWER=", ($result->[5] eq 'NA')?'NA':sprintf ("%.5f", $result->[6]), "\n";
		}
		print "<-----------------Over-represented in tail of ranked list-------------------------->\n";
		@result = sort {$a->[3] <=> $b->[3]} @result;
		for my $result (@result) {
			$result->[3] < 0 or next;
			print "Geneset=$result->[0]\tSize=", sprintf ("%-4d", $result->[1]), "\tES=", sprintf ("%.3f", $result->[2]), "\tNES=", sprintf ("%.3f", $result->[3]), "\tNominalP=", sprintf ("%.5f", $result->[4]), "\tFDR=", ($result->[5] eq 'NA')?'NA':sprintf ("%.3f", $result->[5]), "\tFWER=", ($result->[5] eq 'NA')?'NA':sprintf ("%.5f", $result->[6]), "\n";
		}
	} else {
		print "<-----------------Ranked list of over-represented gene sets/pathways-------------------------->\n";
		for my $result (@result) {
			print "Geneset=$result->[0]\tSize=", sprintf ("%-4d", $result->[1]), "\tES=", sprintf ("%.3f", $result->[2]), "\tNES=", sprintf ("%.3f", $result->[3]), "\tNominalP=", sprintf ("%.5f", $result->[4]), "\tFDR=", ($result->[5] eq 'NA')?'NA':sprintf ("%.3f", $result->[5]), "\tFWER=", ($result->[5] eq 'NA')?'NA':sprintf ("%.5f", $result->[6]), "\n";
		}
	}
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

#this subroutine calculate the arithmatic mean of a list of numbers.
sub mean {
	my @score = @{$_[0]};
	my $sum = 0;
	$sum += $_ for (@score);
	return $sum / @score;
}

=head1 SYNOPSIS

 combine_gsea.pl [arguments] <logfile | ... >

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --setmin <int>		minimum gene in set
            --setmax <int>		maximum gene in set
            --traditional		use traditional GSEA formula
            --variant			reserved argument (not used now)
            --skipfdr			skip FDR calculation (program finishes faster)
            --tolerate			tolerate discordances between input files

 Function: combine GSEA Enrichment Score values calculated in multiple log files 
 to calculate more accurate nominal P, FDR q and FWER P values.

 Example: combine_gsea.pl file1.log file2.log file3.log

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--setmin>

a threshold (default=20) specifying the minimum number of genes in a gene set 
that are also observed in the gene-stat-file. Smaller gene sets are not 
considered by this program.

=item B<--setmax>

a threshold (default=200) specifying the maximum number of genes in a gene set 
that are also observed in the gene-stat-file. Larger gene sets are not 
considered by this program.


=item B<--traditional>

apply the same dichotimization techniques used in the original GSEA publication.

=item B<--variant>

reserved argument that is not currently used but will be implemented in the 
future.

=back

=head1 DESCRIPTION

This program is used to calculate combined statistics from multiple log files 
generated by the calculate_gsea.pl program.

Pathway-based approaches for GWA study that use phenotype permutations are 
generally computationally expensive; therefore, typically one should be running 
the GWA analysis program on multiple computers (such as doing 10 permutations in 
each of 100 CPUs in a computational cluster), then calculate GSEA score for 100 
permutation files separately and generate separte log files, then finally use 
combine_gsea.pl program to combine results together and re-calculate test 
statistics for pathways.