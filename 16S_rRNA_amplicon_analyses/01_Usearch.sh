#!/bin/bash
#SBATCH -J Usearch
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o Usearch.out
#SBATCH -e Usearch.err

module load seqtk
module load fastqc

for file in `cat filelist.txt`
do
# Check read quality
  fastqc -t 24 ${file}.fastq.gz
# Clip 16S-V4 adapters from raw reads
  java -jar /gpfs/data/rbeinart/Software/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar SE -threads 24 -phred33 ${file}.fastq.gz ${file}_singles.fq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/16S_V4_adaptors.fa:2:30:10 MINLEN:75
# Rename read IDs based on sample information so it is possible to link back to origin when merging in Usearch
  seqtk rename ${file}_singles.fq ${file}. > ${file}_renamed.fq
  sed -i "s/ .*//g" ${file}_renamed.fq
done

cat *_renamed.fq > all_merged.fq

# Filter reads, find unique amplicons, denoise and create zOTU / ASV table
usearch11 -fastq_filter all_merged.fq -fastq_maxee_rate 0.001 -fastaout all_filtered.fa -fastq_minlen 150 -fastq_truncqual 20
usearch11 -fastx_uniques all_filtered.fa -fastaout all_uniques.fa -sizeout -relabel Uniq
usearch11 -unoise3 all_uniques.fa -zotus zotus.fa -tabbedout unoise3.txt
usearch11 -otutab all_merged.fq -zotus zotus.fa -sample_delim . -otutabout zotus.txt -mapout zmap.txt -maxrejects 1000 -strand both

# Import results to Qiime2 for taxonomic annotation based on a naive Bayesian classifier corresponding to the particular amplicon; to obtain optimal annotations for all sequences I further compared the zotu.fa file to the new SILVA v138 database online and used the updated taxonomy for downstream analyses
export LC_ALL=en_US.utf8
export LANG=en_US.utf8

source activate qiime2-2019.10

biom convert -i zotus.txt -o zotu-table.biom --to-hdf5 --table-type="OTU table"

qiime tools import --input-path zotu-table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path zotu-table.qza
qiime tools import --input-path zotus.fa --output-path zotu_seqs.qza --type 'FeatureData[Sequence]'

qiime tools import --type 'FeatureData[Sequence]' --input-path silva_132_99_16S.fna --output-path silva_132_99_otus.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path silva_132_99_16S_consensus_taxonomy.txt --output-path silva_132_99-taxonomy.qza
qiime feature-classifier extract-reads --i-sequences silva_132_99_otus.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 150 --p-min-length 100 --p-max-length 400 --o-reads silva_132_99-ref-seqs.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva_132_99-ref-seqs.qza --i-reference-taxonomy silva_132_99-taxonomy.qza --o-classifier silva-132-99-515F-806R-classifier.qza

qiime feature-classifier classify-sklearn --i-classifier silva-132-99-515F-806R-classifier.qza --i-reads zotu_seqs.qza --o-classification taxonomy.qza --p-n-jobs 1

# Export all results for further filtering
qiime tools export --input-path taxonomy.qza --output-path unfiltered
