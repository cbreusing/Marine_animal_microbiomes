#!/bin/bash
#SBATCH -J Tiara
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -o Tiara.out
#SBATCH -e Tiara.err

source activate tiara
module load seqtk

# Sort prokaryotic, eukaryotic and organelle sequences from metagenomic assemblies with Tiara
for file in `cat filelist2.txt`
do
cd ${file}_metaSpades_coassembly
tiara -i scaffolds.fasta -o tiara_results.txt --tf all -t 24 --min_len 1000 --probabilities
grep ">" eukarya_scaffolds.fasta mitochondrion_scaffolds.fasta plastid_scaffolds.fasta unknown_scaffolds.fasta > non_proks.list
sed -i "s/.*>//g" non_proks.list
grep ">" scaffolds.fasta > scaffolds.list
sed -i "s/>//g" scaffolds.list
perl /gpfs/data/rbeinart/cbreusing/Scripts/extract_headers.pl scaffolds.list non_proks.list proks.list
seqtk subseq scaffolds.fasta proks.list > prokarya_scaffolds.fasta
cd ..
done
  

