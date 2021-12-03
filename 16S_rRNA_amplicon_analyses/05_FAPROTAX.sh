#!/bin/bash
#SBATCH -J FAPROTAX
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4g
#SBATCH -o FAPROTAX.out
#SBATCH -e FAPROTAX.err

module load python/3.7.4

/gpfs/data/rbeinart/Software/FAPROTAX_1.2.4/collapse_table.py -i zotu-table-microbiome-filtered_tax.txt -g /gpfs/data/rbeinart/Software/FAPROTAX_1.2.4/FAPROTAX.txt -o microbiome-functions.txt -r microbiome-report.txt -s microbiome_subtables --omit_columns 0 --out_groups2records_table microbiome-groups2records.txt --group_leftovers_as "Miscellaneous" -d "taxonomy" --out_group_overlaps microbiome-group-overlaps.txt -v --normalize_collapsed colums_before_collapsing