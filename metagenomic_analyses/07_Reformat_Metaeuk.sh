#!/bin/bash
#SBATCH -J Reformat
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -o Reformat.out
#SBATCH -e Reformat.err

module load perl

file=("Acanthamunnopsis") # Acanthamunnopsis, Cyclothone, Krill, Munneurycope, Myctophid, Poeobius, Sergestes

# Reformat MetaEuk output and link UniRef IDs, CDS and taxonomic information
grep ">" ${file}_eukarya.fas > ${file}_eukarya.headers
sed -i "s/|/ /g" ${file}_eukarya.headers
sed -e "s/ .*//g" ${file}_eukarya.headers > ${file}_eukarya.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_eukarya.acc > ${file}_eukarya.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_eukarya.headers ${file}_eukarya.ann ${file}_eukarya.annotations.txt

grep ">" ${file}_mito2.fas > ${file}_mito2.headers
sed -i "s/|/ /g" ${file}_mito2.headers
sed -e "s/ .*//g" ${file}_mito2.headers > ${file}_mito2.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_mito2.acc > ${file}_mito2.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_mito2.headers ${file}_mito2.ann ${file}_mito2.annotations.txt

grep ">" ${file}_mito4.fas > ${file}_mito4.headers
sed -i "s/|/ /g" ${file}_mito4.headers
sed -e "s/ .*//g" ${file}_mito4.headers > ${file}_mito4.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_mito4.acc > ${file}_mito4.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_mito4.headers ${file}_mito4.ann ${file}_mito4.annotations.txt

grep ">" ${file}_mito5.fas > ${file}_mito5.headers
sed -i "s/|/ /g" ${file}_mito5.headers
sed -e "s/ .*//g" ${file}_mito5.headers > ${file}_mito5.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_mito5.acc > ${file}_mito5.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_mito5.headers ${file}_mito5.ann ${file}_mito5.annotations.txt

grep ">" ${file}_plastid.fas > ${file}_plastid.headers
sed -i "s/|/ /g" ${file}_plastid.headers
sed -e "s/ .*//g" ${file}_plastid.headers > ${file}_plastid.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_plastid.acc > ${file}_plastid.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_plastid.headers ${file}_plastid.ann ${file}_plastid.annotations.txt

grep ">" ${file}_unk2.fas > ${file}_unk2.headers
sed -i "s/|/ /g" ${file}_unk2.headers
sed -e "s/ .*//g" ${file}_unk2.headers > ${file}_unk2.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_unk2.acc > ${file}_unk2.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_unk2.headers ${file}_unk2.ann ${file}_unk2.annotations.txt

grep ">" ${file}_unk4.fas > ${file}_unk4.headers
sed -i "s/|/ /g" ${file}_unk4.headers
sed -e "s/ .*//g" ${file}_unk4.headers > ${file}_unk4.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_unk4.acc > ${file}_unk4.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_unk4.headers ${file}_unk4.ann ${file}_unk4.annotations.txt

grep ">" ${file}_unk5.fas > ${file}_unk5.headers
sed -i "s/|/ /g" ${file}_unk5.headers
sed -e "s/ .*//g" ${file}_unk5.headers > ${file}_unk5.acc
parallel --pipepart -a /gpfs/data/rbeinart/Databases/uniref90.headers -j 24 --block 100M grep -w -f ${file}_unk5.acc > ${file}_unk5.ann
perl /gpfs/data/rbeinart/cbreusing/Scripts/reformat_metaeuk.pl ${file}_unk5.headers ${file}_unk5.ann ${file}_unk5.annotations.txt

# Extract all TaxIDs
cat *annotations.txt | sed -e "s/.*TaxID=//g" > tmp
cat tmp | sed -e "s/ .*//g" | sort -u > ${file}_all.taxIDs.txt
rm tmp
