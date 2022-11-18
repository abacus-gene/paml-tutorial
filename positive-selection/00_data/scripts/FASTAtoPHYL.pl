#!/usr/bin/perl

use strict;
use warnings;
 
## Script that reads a FASTA alignment and converts it into 
## PHYLIP format.
##
## Usage:
## <path_to_script>/Get_partitions.pl <path_to_aln1> <numtaxa> <lenseq>
##
## NOTE: If the user wants to automatically extract <numtaxa> and <lenseq>
##       they can run the following code: 
##
##       ```  
##         num=$( grep '>' *fasta | wc -l )
##         len=$( sed -n '2,2p' *fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
##         perl ../../FASTAtoPHYL.pl *fasta $num $len
##       ```
##
##       This code snippet assumes that you have a FASTA file as an input 
##       file, thus the second line should be the sequence of the first 
##       taxon.
##
## Contact information: <sandra.ac93@gmail.com>

#IMPORTANT! CREATE ".TXT" FILES IN UNIX, NOT WINDOWS!!!!! IT THEN ADDS THE "\r" CARRIAGE AND IT MAKES THE SCRIPT CRASH...
		
## Open the input files
open (INFILE1, "<$ARGV[0]") or die "Cannot open $ARGV[0] file: $!";
## Get second argument to set output file name
my $outname = $ARGV[0];
chomp($outname);
$outname =~ s/\.fasta//;
$outname =~ s/..*\///;
print "Parsing alignment ".$outname." ... ... \n";
my $outname2 = $outname.".phy";
## Open the output file to save the alignment file 
open(OUT, ">$outname2") or die "Cannot create the output file: $!";
open(OUT2, ">log_lenseq.txt") or die "Cannot create log file: $!";

## Get lines of input file
my @aln1 = <INFILE1>;

## Print PHYLIP header 
print OUT $ARGV[1]."   ".$ARGV[2]."\n\n";

## Create variables 
my $count = 0;
my $species = "";
my $lenseq = 0;

## Loop over all the lines in input file 
foreach my $line (@aln1){

	chomp($line);
	
	# Ommit header
	if( $line =~ /^>/ ){
		$species = $line;
		$species =~ s/>//;
		print OUT $species."      ";
		$count += 1;
	}
	else{
		# Get sequence 
		print OUT $line."\n";
		# Keep track of length
		$lenseq = length($line);
		print OUT2 "Taxa: ".$species."\t"."Length of sequence: ".$lenseq."\n";
	}
}
					
print "Total no. species parsed: ".$count."\nLast seq length: ". $lenseq."\n";

## Close files
close(OUT);
close(OUT2);
close(INFILE1);