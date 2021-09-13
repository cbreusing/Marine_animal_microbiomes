#!/bin/bash
#SBATCH -J graphbin
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=20g
#SBATCH -o GraphBin.out
#SBATCH -e GraphBin.err

source activate graphbin

# After manually evaluating the resulting bins try to refine with GraphBin
for file in `cat filelist2.txt`
do
python /gpfs/data/rbeinart/cbreusing/miniconda3/envs/graphbin/bin/prepresult.py --binned ${file}_metaSpades_coassembly/concoct --output ${file}_metaSpades_coassembly/concoct --prefix ${file}
graphbin --graph ${file}_metaSpades_coassembly/assembly_graph_with_scaffolds.gfa --binned ${file}_metaSpades_coassembly/concoct/${file}_initial_contig_bins.csv --output ${file}_metaSpades_coassembly/GraphBin --prefix ${file} --assembler SPAdes --paths ${file}_metaSpades_coassembly/scaffolds.paths --contigs ${file}_metaSpades_coassembly/prokarya_scaffolds.fasta
done

# All bins were low quality and could not be taxonomically characterized; I therefore decided to not go forward with the binning results and to annotate the combined prokaryotic section
