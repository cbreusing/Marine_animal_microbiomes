#!/bin/bash
#SBATCH -J BLAST
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=10g
#SBATCH -o BLAST.out
#SBATCH -e BLAST.err

module load blast/2.9.0+
module load seqtk

file=("Acanthamunnopsis") # Acanthamunnopsis, Cyclothone, Krill, Munneurycope, Myctophid, Poeobius, Sergestes, Tomopteris, Vampyroteuthis

# Annotate prokaryotic protein-coding genes with UniRef90
blastp -query ${file}_prokarya.faa -db /gpfs/data/rbeinart/Databases/uniref90.fasta -out blast.${file}.uniref90.txt -evalue 1e-10 -num_threads 24 -max_target_seqs 1 -outfmt "6 std stitle"
# Create a new FASTA file with unannotated sequences
perl -anle 'print $F[0]' blast.${file}.uniref90.txt | sort -u > blast.${file}.uniref90.acc
grep ">" ${file}_prokarya.faa | sed -e "s/>//g" | sed -e "s/ .*//g" > ${file}_prokarya.acc
perl /gpfs/data/rbeinart/cbreusing/Scripts/extract_headers.pl ${file}_prokarya.acc blast.${file}.uniref90.acc ${file}_prokarya_nohit.acc
seqtk subseq ${file}_prokarya.faa ${file}_prokarya_nohit.acc > ${file}_prokarya_nohit.faa
# Annotate genes that had no hit in UniRef90 with RefSeq
blastp -query ${file}_prokarya_nohit.faa -db /gpfs/data/shared/databases/refchef_refs/refseq_protein/blast_db/refseq_protein -out blast.${file}.refseq.txt -evalue 1e-10 -num_threads 24 -max_target_seqs 1 -outfmt "6 std stitle staxids"

cat blast.${file}.uniref90.txt blast.${file}.refseq.txt | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 > blast.${file}.topHit.txt