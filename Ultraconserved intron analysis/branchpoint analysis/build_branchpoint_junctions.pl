#!/usr/bin/perl
#
# build_branchpoint_junctions_v2.pl 
# graveley 8/25/13
# lawson 12/2/22

#This program creates "splice junctions" that would be observed at the 5' splice site:lariat intron junctions. 

#[Exon 1]->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->A->->->->->->->[Exon 2]
#       |-----1-----|------2----|-----3-----|-----4-----|-----5-----|-----6-----|-----7-----|

#		|->->->->A->->->->->->->|->->->->->->->
#		|-----6-----|-----7-----||-----1------|

#Steps

#	1. Extract target intron sequence
#
#	2. For each 100 nt window, take 1ast 100 nt of intron and fuse to first 100 nt of the intron. Then redo moving window 1 nt upstream for the downstream portion of the intron. 


use strict;
use FileHandle;

my $FILE = $ARGV[0]; #Junction Half Sequences
my $JUNCTION_LENGTH = $ARGV[1]; #Length of the segment for each half of the junction. Should be Readlength minus 6
my $OUT = new FileHandle ">$FILE"."_ref.fa"; #Branchpoint junction sequences in fasta format

my $line; 
my $header;
my $intron;
my $first_half;
my $second_half;
my $interval_1;
my $interval_2;
my $new_header;

open(FASTA, "$FILE") or die("can't open $FILE"); #INPUT FILE: FASTA file containing introns to be processed

printf("Processing $FILE \n");

  while ($line = <FASTA>){        #read fasta header
        chop($line);
        chop($line);
        $header = $line;
        $line = <FASTA>;
        chop($line);
        chop($line);
        $intron = $line;
        $first_half = substr($intron,0,$JUNCTION_LENGTH);
        for(my $i = 0; $i < 500; $i++ ){
        $interval_1 = -$JUNCTION_LENGTH - $i;
        last if ($JUNCTION_LENGTH > length($intron) + $interval_1);
        $interval_2 = $JUNCTION_LENGTH;
        $second_half = substr($intron,$interval_1,$interval_2);
        $new_header = "$header"."$interval_1";      
        $OUT -> print("$new_header\n$second_half"."$first_half\n");
		}
}
