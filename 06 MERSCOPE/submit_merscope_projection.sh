#!/bin/bash
#SBATCH --job-name=schpf
#SBATCH --output=outs_schpf_merscope-%a.txt
#SBATCH --error=outs_schpf_merscope-%a.txt
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-12:1

DIRS=($(ls -d /path/to/projects/microglia_atlas/merscope/*/ | grep "CP"))

HOME_PATH=${DIRS[$SLURM_ARRAY_TASK_ID-1]}'projection'
MOD_PATH='/path/to/projects/microglia_atlas/scHPF/00_model_fit'
SCHPF_PATH='/path/to/tools/consensus_scHPF_wrapper-main'

source ~/miniconda3/bin/activate schpf_p37v2

cd $HOME_PATH

python $SCHPF_PATH$'/downsample_loom.py' \
-l $HOME_PATH'/data.loom' \
-a 'orig.ident' \
-o $HOME_PATH'/data_dsamp_fullgenes.loom' \
-n 854

scHPF prep-like \
-i $HOME_PATH'/data_dsamp_fullgenes.loom' \
-r $MOD_PATH'/V3_full.v3tech.ss.ds.train.genes.txt' \
-p data \
-o $HOME_PATH'/data_out'

scHPF project \
-m $MOD_PATH'/V3_full.v3tech.ss.ds.consensus.final.joblib' \
-i $HOME_PATH'/data_out/data.filtered.mtx' \
-o $HOME_PATH'/data_out'

scHPF score \
-m $HOME_PATH'/data_out/V3_full.v3tech.ss.ds.consensus.final.proj.joblib' \
-o $HOME_PATH'/data_out' \
-p proj \
-g $MOD_PATH'/V3_full.v3tech.ss.ds.train.genes.txt'