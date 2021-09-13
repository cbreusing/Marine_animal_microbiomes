#!/bin/bash
#SBATCH -J MetaEuk
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=150g
#SBATCH -o MetaEuk.out
#SBATCH -e MetaEuk_.err

source activate metaeuk

file=("Acanthamunnopsis") # Acanthamunnopsis, Cyclothone, Krill, Munneurycope, Myctophid, Poeobius, Sergestes

# Predict and annotate eukaryotic and organelle CDS and assign taxonomy (confidence in taxonomic assignments was not high enough, I therefore used a separate script to extract the lineage information from the annotation hit)
mmseqs createdb ${file}_metaSpades_coassembly/eukarya_scaffolds.fasta ${file}_metaSpades_coassembly/eukarya_scaffolds.DB
mmseqs createdb	${file}_metaSpades_coassembly/mitochondrion_scaffolds.fasta ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB
mmseqs createdb ${file}_metaSpades_coassembly/plastid_scaffolds.fasta ${file}_metaSpades_coassembly/plastid_scaffolds.DB
mmseqs createdb ${file}_metaSpades_coassembly/unknown_scaffolds.fasta ${file}_metaSpades_coassembly/unknown_scaffolds.DB

metaeuk predictexons ${file}_metaSpades_coassembly/eukarya_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/eukarya_results.DB ${file}_metaSpades_coassembly/tmp -s 7.5 --use-all-table-starts 1 --metaeuk-eval 0.0001 -e 100 --min-length 40
metaeuk predictexons ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/mitochondrion_results.DB ${file}_metaSpades_coassembly/tmp -s 7.5 --use-all-table-starts 1 --metaeuk-eval 0.0001 -e 100 --min-length 40
metaeuk predictexons ${file}_metaSpades_coassembly/plastid_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/plastid_results.DB ${file}_metaSpades_coassembly/tmp -s 7.5 --use-all-table-starts 1 --metaeuk-eval 0.0001 -e 100 --min-length 40
metaeuk predictexons ${file}_metaSpades_coassembly/unknown_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/unknown_results.DB ${file}_metaSpades_coassembly/tmp -s 7.5 --use-all-table-starts 1 --metaeuk-eval 0.0001 -e 100 --min-length 40

metaeuk reduceredundancy ${file}_metaSpades_coassembly/eukarya_results.DB ${file}_metaSpades_coassembly/eukarya_predResults.DB ${file}_metaSpades_coassembly/eukarya_predGroups.DB
metaeuk reduceredundancy ${file}_metaSpades_coassembly/mitochondrion_results.DB ${file}_metaSpades_coassembly/mitochondrion_predResults.DB ${file}_metaSpades_coassembly/mitochondrion_predGroups.DB
metaeuk reduceredundancy ${file}_metaSpades_coassembly/plastid_results.DB ${file}_metaSpades_coassembly/plastid_predResults.DB ${file}_metaSpades_coassembly/plastid_predGroups.DB
metaeuk reduceredundancy ${file}_metaSpades_coassembly/unknown_results.DB ${file}_metaSpades_coassembly/unknown_predResults.DB ${file}_metaSpades_coassembly/unknown_predGroups.DB

metaeuk unitesetstofasta ${file}_metaSpades_coassembly/eukarya_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/eukarya_predResults.DB ${file}_metaSpades_coassembly/${file}_eukarya --translation-table 1
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/mitochondrion_predResults.DB ${file}_metaSpades_coassembly/${file}_mito2 --translation-table 2
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/mitochondrion_predResults.DB ${file}_metaSpades_coassembly/${file}_mito4 --translation-table 4
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/mitochondrion_predResults.DB ${file}_metaSpades_coassembly/${file}_mito5 --translation-table 5
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/plastid_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/plastid_predResults.DB ${file}_metaSpades_coassembly/${file}_plastid --translation-table 11
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/unknown_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/unknown_predResults.DB ${file}_metaSpades_coassembly/${file}_unk2 --translation-table 2
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/unknown_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/unknown_predResults.DB ${file}_metaSpades_coassembly/${file}_unk4 --translation-table 4
metaeuk unitesetstofasta ${file}_metaSpades_coassembly/unknown_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/unknown_predResults.DB ${file}_metaSpades_coassembly/${file}_unk5 --translation-table 5

metaeuk groupstoacc ${file}_metaSpades_coassembly/eukarya_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/eukarya_predGroups.DB ${file}_metaSpades_coassembly/eukarya_predGroups.tsv
metaeuk taxtocontig ${file}_metaSpades_coassembly/eukarya_scaffolds.DB ${file}_metaSpades_coassembly/${file}_eukarya.fas ${file}_metaSpades_coassembly/${file}_eukarya.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_eukarya ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 1

metaeuk groupstoacc ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/mitochondrion_predGroups.DB ${file}_metaSpades_coassembly/mitochondrion_predGroups.tsv
metaeuk taxtocontig ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB ${file}_metaSpades_coassembly/${file}_mito2.fas ${file}_metaSpades_coassembly/${file}_mito2.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_mito2 ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 2
metaeuk taxtocontig ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB ${file}_metaSpades_coassembly/${file}_mito4.fas ${file}_metaSpades_coassembly/${file}_mito4.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_mito4 ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 4
metaeuk taxtocontig ${file}_metaSpades_coassembly/mitochondrion_scaffolds.DB ${file}_metaSpades_coassembly/${file}_mito5.fas ${file}_metaSpades_coassembly/${file}_mito5.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_mito5 ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 5

metaeuk groupstoacc ${file}_metaSpades_coassembly/plastid_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/plastid_predGroups.DB ${file}_metaSpades_coassembly/plastid_predGroups.tsv
metaeuk taxtocontig ${file}_metaSpades_coassembly/plastid_scaffolds.DB ${file}_metaSpades_coassembly/${file}_plastid.fas ${file}_metaSpades_coassembly/${file}_plastid.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_plastid ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 11

metaeuk groupstoacc ${file}_metaSpades_coassembly/unknown_scaffolds.DB /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/unknown_predGroups.DB ${file}_metaSpades_coassembly/unknown_predGroups.tsv
metaeuk taxtocontig ${file}_metaSpades_coassembly/unknown_scaffolds.DB ${file}_metaSpades_coassembly/${file}_unk2.fas ${file}_metaSpades_coassembly/${file}_unk2.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_unk2 ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 2
metaeuk taxtocontig ${file}_metaSpades_coassembly/unknown_scaffolds.DB ${file}_metaSpades_coassembly/${file}_unk4.fas ${file}_metaSpades_coassembly/${file}_unk4.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_unk4 ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 4
metaeuk taxtocontig ${file}_metaSpades_coassembly/unknown_scaffolds.DB ${file}_metaSpades_coassembly/${file}_unk5.fas ${file}_metaSpades_coassembly/${file}_unk5.headersMap.tsv /gpfs/data/rbeinart/Databases/uniref90.DB ${file}_metaSpades_coassembly/${file}_unk5 ${file}_metaSpades_coassembly/tmp -s 7.5 --majority 0.5 --tax-lineage 1 --lca-mode 4 --translation-table 5




