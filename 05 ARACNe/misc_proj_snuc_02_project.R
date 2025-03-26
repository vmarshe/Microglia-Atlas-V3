#!/usr/bin/Rscript

#===============================================================================
# 
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(Matrix)
library(SeuratObject)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(logger)
library(uwot)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
file_names <- c(snuc = "ROSMAP_singlenuc", multiome = "multiome", kellis = "kellis_mic")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_project_scHPF.R"))

# UMAP
#-------------------------------------------------------------------------------
umap_uwot <- load_uwot(file.path(home_path, "umap/umap_uwot.rds"))

# PROCESS DATA
#-------------------------------------------------------------------------------
for(i in c("kellis")){ #,"snuc",  "multiome"
  
  res_path <- file.path(home_path, paste0("proj_", i, "_validation"))
  
  loom <- Connect(file.path(res_path, "input", "data_dsamp_fullgenes.loom"),  mode = "r",  force = T)
  loom <- as.Seurat(x = loom)
  
  cell_names <-  Cells(loom)
  gene_names <-  rownames(loom)
 
  project_scHPF(loom_file = file.path(res_path, "input", "data.loom"), 
                name = i, 
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
