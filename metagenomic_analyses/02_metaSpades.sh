#!/bin/bash
#SBATCH -J metaSpades
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o metaSpades.out
#SBATCH -e metaSpades.err

file=("Acanthamunnopsis") # Acanthamunnopsis, Cyclothone, Krill, Munneurycope, Myctophid, Poeobius, Sergestes, Tomopteris, Vampyroteuthis

# Combine all reads across samples for each species/organism to produce a coassembly
cat ${file}*_R1_clean.fastq > ${file}_merged_R1_clean.fastq
cat ${file}*_R2_clean.fastq > ${file}_merged_R2_clean.fastq
cat ${file}*_singles_clean.fastq > ${file}_merged_singles_clean.fastq

metaspades.py -1 ${file}_merged_R1_clean.fastq -2 ${file}_merged_R2_clean.fastq -s ${file}_merged_singles_clean.fastq -k 21,31,41,51,61,71,81,91,101,111,121 -t 24 -o ${file}_metaSpades_coassembly




