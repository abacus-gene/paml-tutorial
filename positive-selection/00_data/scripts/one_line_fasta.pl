#!/usr/bin/perl

use strict;
use warnings;

open (INP_FASTA, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
## Open the output file to save the reference genome ID matching the 
## input strain
my $cut_name = $ARGV[0];
$cut_name =~ s/\..*//;
my $out_name = "$cut_name" . "_one_line.fa";
open(OUT, ">>$out_name") or die "Cannot create the output file: $!";

my @oneline_fasta;
	
my $sequence;
my $linecount = 1;

while (<INP_FASTA>){
	chomp;
	if (m/^>/){ 							# If the line is a header line
		if ($linecount == 1){			# If this is the first header line, 
			push (@oneline_fasta, $_);		# just add the header to the oneline_fasta array
			$linecount++;					# Increase the linecount
		}	
		else{
			push (@oneline_fasta, $sequence); # print the concatenated sequence to output
			push (@oneline_fasta, $_);	# print the new header to oneline_fasta
			$sequence = '';			# empty the sequence string
		}
	}	
	else{ 							# not a header line? - must be sequence
    	chomp; 						# remove newline at end
    	$sequence .= $_; 				# append additional sequence
	}
}

push (@oneline_fasta, $sequence);			# print the last sequence to the file

foreach my $element (@oneline_fasta){
	print OUT "$element\n";
}

close OUT;
close INP_FASTA;