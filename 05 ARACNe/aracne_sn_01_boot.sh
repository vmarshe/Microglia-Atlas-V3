#!/bin/bash
#$ -cwd
#$ -l h_vmem=100G,h_rt=24:00:00
#$ -q csg2.q 
#$ -l t_pri
#$ -N aracne
#$ -e aracne_boot.out
#$ -o aracne_boot.out

for i in {1..200}
do
java -Xmx75G -jar ~/tools/ARACNe-AP/dist/aracne.jar -e ~/projects/microglia_atlas/aracne_sn/metacell_matrix.txt -o ~/projects/microglia_atlas/aracne_sn/ \
--tfs ~/projects/microglia_atlas/aracne_sn/tfs.symbols.txt \
--pvalue 1E-8 \
--seed $i
done

