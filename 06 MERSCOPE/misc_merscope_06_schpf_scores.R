#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE: Project MERSCOPE data into the reference UMAP space.
#
#===============================================================================
# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(ggrepel)
library(ggsignif)
library(ggtext)
library(logger)
library(patchwork)
library(schex)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tidyverse)
library(uwot)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "/mnt/vast/hpc/homes/vsm2116/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/mglia")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_project_scHPF.R"))

#===============================================================================
#
# PROJECTION
#
#===============================================================================
run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))
umap_uwot <- load_uwot(file.path(home_path, "umap/umap_uwot.rds"))

for(folder in run_metadata$run_name){
  
  res_path <- file.path(home_path, "merscope", folder, "projection")
  
  log_appender(appender_tee(file = file.path(res_path, "session_info.log")))
  log_threshold(TRACE)
  
  loom <- Connect(file.path(res_path, "data_dsamp_fullgenes.loom"),  mode = "r",  force = T)
  cell_names <- Cells(as.Seurat(x = loom))
  loom$close_all()
  
  gene_names <- read_table(file.path(res_path, "gene_input.txt"), col_names = F)$X1
  
  project_scHPF(loom_file = file.path(res_path, "data.loom"),
                name = "data",
                factors = paste0("scHPF_", c(1:5, 7:11, 14:26)),
                cell_names = cell_names,
                gene_names = gene_names,
                key = "scHPF",
                assay = "RNA",
                umap = umap_uwot,
                cell_score_file = file.path(res_path, "data_out", "proj.cell_score.txt"),
                gene_score_file = file.path(res_path, "data_out", "proj.gene_score.txt"),
                res_path = res_path)
  
}
#===============================================================================
#
# COMBINE DATA
#
#===============================================================================
umap <- c()
schpf <- c()
for(i in 1:nrow(run_metadata)){
  data <- readRDS(file.path(home_path, "merscope", run_metadata$run_name[i], "projection/data/data.rds"))
  data@meta.data$orig.ident <- run_metadata$batch_code[i]
  umap <- rbind(umap, data@reductions$schpf.umap@cell.embeddings)
  schpf <- rbind(schpf, data@reductions$scHPF@cell.embeddings)
}
umap <- CreateDimReducObject(umap[Cells(mglia), ], key="umap")
schpf <- CreateDimReducObject(schpf[Cells(mglia), ], key="scHPF")

saveRDS(umap, file.path(home_path, "merscope", "mglia","mglia.umap.rds"))
saveRDS(schpf, file.path(home_path, "merscope", "mglia","mglia.schpf.umap.rds"))

#===============================================================================
#
# DEFINE HIGH vs. LOW EXPRESSION CLASSES
#
#===============================================================================
scores <- schpf@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column("cell_names") %>%
  dplyr::select(cell_names, scHPF_1:scHPF_26) %>%
  mutate(scHPF_20_class = ifelse(scHPF_20 > median(scHPF_20)+(2*mad(scHPF_20)), "Mic.20.high", "Mic.20.low"),
         scHPF_26_class = ifelse(scHPF_26 > median(scHPF_26)+(2*mad(scHPF_26)), "Mic.26.high", "Mic.26.low"))

saveRDS(scores, file.path(home_path, "merscope/mglia/schpf_scores.rds"))
