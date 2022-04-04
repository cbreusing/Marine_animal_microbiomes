#!/bin/bash
#SBATCH -J SHOGUN
#SBATCH -t 24:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -o SHOGUN.out
#SBATCH -e SHOGUN.err

shi7 -i shi7_in -o shi7_out --flash True --strip_delim _,2 -t 24 -filter_l 75 -filter_q 30 -trim_q 25

shogun --log debug pipeline -a utree -i shi7_out/combined_seqs.fna -d /gpfs/data/rbeinart/Databases/SHOGUN -o SHOGUN_utree -l all --function -t 24