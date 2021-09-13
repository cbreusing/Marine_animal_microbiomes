#!/bin/bash
#SBATCH -J SeqFilter
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o SeqFilter.out
#SBATCH -e SeqFilter.err

module load samtools/1.12
module load seqtk
module load bowtie2/2.3.5.1
module load fastqc

for file in `cat filelist.txt`
do
# Check read quality
fastqc -t 24 ${file}_R1.fastq.gz
fastqc -t 24 ${file}_R2.fastq.gz
# Clip adapters and quality trim reads
java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 ${file}_R1.fastq.gz ${file}_R2.fastq.gz ${file}_R1_paired.fq ${file}_R1_unpaired.fq ${file}_R2_paired.fq ${file}_R2_unpaired.fq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/Illumina.fa:2:30:10 HEADCROP:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:75
# Remove common contaminants, such as human and PhiX DNA sequences
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 ${file}_R1_paired.fq -2 ${file}_R2_paired.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.cont.sorted.bam
samtools view -@ 24 -f12 ${file}.bowtie2.cont.sorted.bam > ${file}.cont.unmapped.sam
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -U ${file}_R1_unpaired.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.cont.singles.sorted.bam
samtools view -@ 24 -f4 ${file}.bowtie2.cont.singles.sorted.bam > ${file}.cont.singles.unmapped.sam
cut -f1 ${file}.cont.unmapped.sam | sort | uniq > ${file}.cont.unmapped_ids.lst
cut -f1 ${file}.cont.singles.unmapped.sam | sort | uniq > ${file}.cont.singles.unmapped_ids.lst
seqtk subseq ${file}_R1_paired.fq ${file}.cont.unmapped_ids.lst > ${file}_R1_clean.fastq
seqtk subseq ${file}_R2_paired.fq ${file}.cont.unmapped_ids.lst > ${file}_R2_clean.fastq
seqtk subseq ${file}_R1_unpaired.fq ${file}.cont.singles.unmapped_ids.lst > ${file}_singles_clean.fastq
done

