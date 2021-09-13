#!/bin/bash
#SBATCH -J phyloFlash
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o phyloFlash.out
#SBATCH -e phyloFlash.err

module load bbmap
module load perl
module load bedtools
module load mafft

# Use phyloFlash to extract 16S sequence information from metagenomic datasets
for file in `cat filelist2.txt`
do
phyloFlash.pl -lib ${file} -read1 ${file}_merged_R1_clean.fastq -read2 ${file}_merged_R2_clean.fastq -readlength 150 -CPUs 24
metaspades.py -1 ${file}.${file}_merged_R1_clean.fastq.SSU.1.fq -2 ${file}.${file}_merged_R1_clean.fastq.SSU.2.fq -t 24 -o ${file}_SSU
done

