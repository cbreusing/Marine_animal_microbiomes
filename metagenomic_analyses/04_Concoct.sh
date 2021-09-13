#!/bin/bash
#SBATCH -J CONCOCT
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -o CONCOCT.out
#SBATCH -e CONCOCT.err

source activate concoct
module load samtools/1.12
module load bowtie2/2.3.5.1

file=("Acanthamunnopsis") # Acanthamunnopsis, Cyclothone, Krill, Munneurycope, Myctophid, Poeobius, Sergestes, Tomopteris, Vampyroteuthis

# Try to bin prokaryotic contigs with CONCOCT
cut_up_fasta.py ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta -c 10000 -o 0 --merge_last -b ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.bed > ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta
bowtie2-build ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta
  for i in `cat ${file}_metaSpades_coassembly/samples.txt`
  do
      bowtie2 -p 24 -x ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta -1 ${i}_R1_clean.fastq -2 ${i}_R2_clean.fastq -U ${i}_singles_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}_metaSpades_coassembly/${i}_10K.bowtie2.sorted.bam
      samtools index ${file}_metaSpades_coassembly/${i}_10K.bowtie2.sorted.bam ${file}_metaSpades_coassembly/${i}_10K.bowtie2.sorted.bam.bai
  done   
concoct_coverage_table.py ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.bed ${file}_metaSpades_coassembly/*10K.bowtie2.sorted.bam > ${file}_metaSpades_coassembly/coverage_table.tsv
CONTIGS=`grep -c ">" ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta`
fasta_to_features.py ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta ${CONTIGS} 4 ${file}_metaSpades_coassembly/kmer_4.csv
mkdir ${file}_metaSpades_coassembly/concoct
concoct --composition_file ${file}_metaSpades_coassembly/prokarya_scaffolds_10K.fasta --coverage_file ${file}_metaSpades_coassembly/coverage_table.tsv -t 1 -b ${file}_metaSpades_coassembly/concoct/${file} -l 1000
merge_cutup_clustering.py ${file}_metaSpades_coassembly/concoct/${file}_clustering_gt1000.csv > ${file}_metaSpades_coassembly/concoct/${file}_clustering_merged.csv
extract_fasta_bins.py ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta ${file}_metaSpades_coassembly/concoct/${file}_clustering_merged.csv --output_path ${file}_metaSpades_coassembly/concoct
