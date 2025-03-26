#!/bin/bash
#$ -cwd
#$ -l h_vmem=100G,h_rt=24:00:00
#$ -q csg2.q 
#$ -l t_pri
#$ -N aracne
#$ -e aracne_con.out
#$ -o aracne_con.out

java -Xmx75G -jar ~/tools/ARACNe-AP/dist/aracne.jar -o ~/projects/microglia_atlas/aracne_sn/ \
 --consolidate