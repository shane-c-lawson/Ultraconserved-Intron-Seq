#!/bin/bash

#SBATCH --qos=general
#SBATCH --mem 100G
#SBATCH --nodes=2

bowtie2 --no-unal -p 16 -x upstream_introns -1 R1.fastq.gz -2 R2.fastq.gz -S upstream_introns_hits.sam
