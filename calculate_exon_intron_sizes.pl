#!/usr/bin/env perl

use warnings;
use strict;

# File: calculate_exon_intron_sizes.pl
# Author: Sean D. Gallaher
# Version: v1.0
# Date: 2022-SEP-01
# Description: 
# This is a script to 
# go through a gff3 file
# and calculate the sizes
# of exons and introns.

# Help message
my $help = "\nThis is a script to go through a gff3 file and calculate the sizes of exons and introns.\nPlease supply a gff3 file and the number of characters in the gene IDs as two arguements.\n\n";

# Check for requested info

my $gff = shift @ARGV;

if (!defined $gff) {
	die $help;
}

my $charLen = shift @ARGV;

if (!defined $charLen) {
	die $help;
}

# open the gff3 file

open (IN, "<", $gff)
	or die "Cannot open $gff: $!\n";




# go through file and put info in a hash

my %master;

while (my $line = <IN>) {
	my $gid;
	chomp $line;
	unless ($line =~ m/^#/) {
		my @lineArray = split(/\t/, $line);
		if ($lineArray[2] eq "gene") {
			my $desc = $lineArray[8];
			$desc =~ m/ID=(.{$charLen})/g ;
			$gid = $1;
			# print OUT "$gid\n";
			${master}{$gid}{'exon_count'} = 0;
			${master}{$gid}{'exon_total'} = 0;
			${master}{$gid}{'intron_total'} = 0;
			${master}{$gid}{'strand'} = $lineArray[6];
		}
		elsif ($lineArray[2] eq "CDS") {
			my $desc = $lineArray[8];
			$desc =~ m/.*ID=([^;]+).*/g ;
			my $tid = $1;
			#print OUT "$tid\n";
			$tid =~ m/(.{$charLen}).+/g;
			my $parent = $1;
			#print OUT "$tid\t$parent\n";
			my $oldExonCount = ${master}{$parent}{'exon_count'};
			my $newExonCount = $oldExonCount + 1 ;
			${master}{$parent}{'exon_count'} = $newExonCount;
			my $start = $lineArray[3];
			my $end = $lineArray[4];
			my $exonSize = 	$end - $start + 1;
			${master}{$parent}{'exons'}{$newExonCount}{'length'} = $exonSize;
			${master}{$parent}{'exons'}{$newExonCount}{'start'} = $start;
			${master}{$parent}{'exons'}{$newExonCount}{'end'} = $end;
			${master}{$parent}{'exon_total'} += $exonSize;
			my $exonSizesString = ${master}{$parent}{'exon_sizes'};
			if (defined $exonSizesString) {
				$exonSizesString = $exonSizesString . "," . $exonSize;
			}
			else {
				$exonSizesString = $exonSize;
			}
			${master}{$parent}{'exon_sizes'} = $exonSizesString;
			
		
		}	
	} 
}

# create output file
my $outfile = "exon_intron_size_report.tsv";
open (OUT, ">", $outfile) 
	or die "Cannot create $outfile: $!\n"; 
# go through %master and output the results

print OUT "#GID\texon_count\tintron_count\ttotal_exon_size\ttotal_intron_size\texon_sizes\tintron_sizes\n";

foreach my $gid (sort keys %master) {
	my $exonCount = ${master}{$gid}{'exon_count'};
	my $intronCount = $exonCount - 1;
	my $exonTotal = ${master}{$gid}{'exon_total'};
	my $exonSizes = ${master}{$gid}{'exon_sizes'};
	my $strand = ${master}{$gid}{'strand'};
	my $intronTotal = 0;
	my $intronSizes ;
	if ( $intronCount > 0 ) {
		if ($strand eq "+") {
			for (my $i = 1; $i < $exonCount; $i++){
				my $j = $i + 1;
				my $intronStart = ${master}{$gid}{'exons'}{$i}{'end'} + 1;
				my $intronEnd = ${master}{$gid}{'exons'}{$j}{'start'} - 1;
				my $intronSize = $intronEnd - $intronStart + 1; 
				$intronTotal += $intronSize;
				if (defined $intronSizes) {
					$intronSizes = $intronSizes . "," . $intronSize;
				}
				else {
					$intronSizes = $intronSize;
				}
			}
		}
		elsif ($strand eq "-") {
			for (my $k = $exonCount; $k > 1; $k--) {
				my $l = $k - 1;
				my $intronStart = ${master}{$gid}{'exons'}{$k}{'end'} + 1;
				my $intronEnd = ${master}{$gid}{'exons'}{$l}{'start'} - 1;
				my $intronSize = $intronEnd - $intronStart + 1; 
				$intronTotal += $intronSize;
				if (defined $intronSizes) {
					$intronSizes = $intronSizes . "," . $intronSize;
				}
				else {
					$intronSizes = $intronSize;
				}
				
			}	
		}
	}
	if (!defined $intronSizes) {
		$intronSizes = "0";
	}
	print OUT "$gid\t$exonCount\t$intronCount\t$exonTotal\t$intronTotal\t$exonSizes\t$intronSizes\n";
}


# close gff3 and output files
close (IN);
close (OUT);


