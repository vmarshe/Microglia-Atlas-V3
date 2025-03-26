#!/bin/bash
#$ -cwd
#$ -l h_vmem=100G,h_rt=24:00:00
#$ -q csg2.q 
#$ -l t_pri
#$ -N aracne_mi
#$ -e aracne_mi.out
#$ -o aracne_mi.out

java -version

java -Xmx75G -jar ~/tools/ARACNe-AP/dist/aracne.jar -e ~/projects/microglia_atlas/aracne_sn/metacell_matrix.txt -o ~/projects/microglia_atlas/aracne_sn/ \
--tfs ~/projects/microglia_atlas/aracne_sn/tfs.symbols.txt \
--pvalue 1E-8 \
--seed 1234 \
--calculateThreshold