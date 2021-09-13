#!/bin/bash
#SBATCH -J Qiime2-Filter
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o Qiime2.out
#SBATCH -e Qiime2.err

# After manual filtering of ASVs into target categories (microbiome, diet, contaminants/host) create new (filtered) ASV tables and tree files
export LC_ALL=en_US.utf8
export LANG=en_US.utf8

source activate qiime2-2019.10
module load seqtk

# Microbiome
qiime feature-table filter-features --i-table zotu-table.qza --m-metadata-file microbiome.txt --o-filtered-table zotu-table-microbiome.qza
qiime tools export --input-path zotu-table-microbiome.qza --output-path microbiome
biom convert -i microbiome/feature-table.biom -o microbiome/zotu-table-microbiome.txt --to-tsv

seqtk subseq zotus.fa microbiome.txt > zotus_microbiome.fa
qiime tools import --input-path zotus_microbiome.fa --output-path zotu_microbiome-seqs.qza --type 'FeatureData[Sequence]'

qiime alignment mafft --i-sequences zotu_microbiome-seqs.qza --o-alignment aligned-microbiome.qza
qiime alignment mask --i-alignment aligned-microbiome.qza --o-masked-alignment masked-aligned-microbiome.qza
qiime phylogeny fasttree --i-alignment masked-aligned-microbiome.qza --o-tree unrooted-tree-microbiome.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree-microbiome.qza --o-rooted-tree rooted-tree-microbiome.qza

qiime tools export --input-path zotu_microbiome-seqs.qza --output-path microbiome
qiime tools export --input-path rooted-tree-microbiome.qza --output-path microbiome

# Diet
qiime feature-table filter-features --i-table zotu-table.qza --m-metadata-file diet.txt --o-filtered-table zotu-table-diet.qza
qiime tools export --input-path zotu-table-diet.qza --output-path diet
biom convert -i diet/feature-table.biom -o diet/zotu-table-diet.txt --to-tsv

seqtk subseq zotus.fa diet.txt > zotus_diet.fa
qiime tools import --input-path zotus_diet.fa --output-path zotu_diet-seqs.qza --type 'FeatureData[Sequence]'

qiime alignment mafft --i-sequences zotu_diet-seqs.qza --o-alignment aligned-diet.qza
qiime alignment mask --i-alignment aligned-diet.qza --o-masked-alignment masked-aligned-diet.qza
qiime phylogeny fasttree --i-alignment masked-aligned-diet.qza --o-tree unrooted-tree-diet.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree-diet.qza --o-rooted-tree rooted-tree-diet.qza

qiime tools export --input-path zotu_diet-seqs.qza --output-path diet
qiime tools export --input-path rooted-tree-diet.qza --output-path diet

# Convert to biom format after adding taxonomy information manually
biom convert -i microbiome/zotu-table-microbiome_tax.txt -o microbiome/zotu-table-microbiome_tax.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
biom convert -i diet/zotu-table-diet_tax.txt -o diet/zotu-table-diet_tax.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

