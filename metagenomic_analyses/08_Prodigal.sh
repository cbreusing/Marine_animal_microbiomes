#!/bin/bash
#SBATCH -J Prodigal
#SBATCH -t 100:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -o Prodigal.out
#SBATCH -e Prodigal.err

module load prodigal

# Predict protein-coding genes in prokaryotic fraction of each metagenome
for file in `cat filelist2.txt`
do
prodigal -a ${file}_metaSpades_coassembly/${file}_prokarya.faa -d ${file}_metaSpades_coassembly/${file}_prokarya.fna -f gff -i ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta -p meta -o ${file}_metaSpades_coassembly/${file}_prokarya.genes -s ${file}_metaSpades_coassembly/${file}_prokarya.all_genes 
done


