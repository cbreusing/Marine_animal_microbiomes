#!/bin/bash
#SBATCH -J metaBAT
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -o metaBAT.out
#SBATCH -e metaBAT.err

source activate metabat2
module load samtools/1.12
module load bowtie2/2.3.5.1
module load hmmer
module load idba

# Try to bin prokaryotic contigs with metaBAT2 and MaxBin2 - this did not work well for this particular dataset
for file in `cat filelist2.txt`
do
bowtie2-build ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta
bowtie2 -p 24 -x ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta -1 ${file}_merged_R1_clean.fastq -2 ${file}_merged_R2_clean.fastq -U ${file}_merged_singles_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}_metaSpades_coassembly/${file}.bowtie2.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth ${file}_metaSpades_coassembly/depth.txt ${file}_metaSpades_coassembly/${file}.bowtie2.sorted.bam
metabat2 -i ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta -a ${file}_metaSpades_coassembly/depth.txt -m 1500 -o ${file}_metaSpades_coassembly/metabat/${file} --unbinned -t 24
mkdir ${file}_metaSpades_coassembly/maxbin
run_MaxBin.pl -contig ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta -reads ${file}_merged_R1_clean.fastq -reads2 ${file}_merged_R2_clean.fastq -out ${file}_metaSpades_coassembly/maxbin/${file} -thread 24 -min_contig_length 1000
done


