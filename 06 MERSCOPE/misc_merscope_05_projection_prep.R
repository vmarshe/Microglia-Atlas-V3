#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE: Prep object for scHPF scoring.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib")
library(matrixStats)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(SeuratDisk)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas/merscope"

# DATA
#-------------------------------------------------------------------------------
run_metadata <- read_csv(file.path(home_path, "batch_names.csv"))
cell_type_annotations <- read_rds(file.path(home_path, "merged/cell_type_annotations.rds")) %>%
  dplyr::select(cell_names, orig.ident, cell_type)

schpf <- read_rds("~/projects/microglia_atlas/data/scHPF_v3_reduc_FULL.rds")
genes <- rownames(schpf@feature.loadings)

# GET COUNT DATA
#-------------------------------------------------------------------------------
for(j in 1:nrow(run_metadata)){
  
  res_path <- file.path(home_path, run_metadata$run_name[j], "projection")
  if(!dir.exists(res_path)) { dir.create(res_path)}
  
  data <- readRDS(file.path(home_path, run_metadata$run_name[j], "data_proc.rds"))
  data@meta.data <- data@meta.data %>%
    left_join(cell_type_annotations, by = c("fov.y"="orig.ident", "cell_names"))
  rownames(data@meta.data) <- data@meta.data$cell_names
  
  data <- subset(data, subset = cell_type %in% "Microglia")
  
  # IMPUTE MISSING GENES
  #-----------------------------------------------------------------------------
  missing_genes <- genes[!genes %in% rownames(data)]
  missing_mat <- matrix(data = 0, nrow = length(missing_genes), ncol = ncol(data))
  rownames(missing_mat) <- missing_genes
  colnames(missing_mat) <- colnames(data)
  
  res <- data.frame(genes = missing_genes)
  for(i in colnames(schpf)){
    ind <- data.frame(genes = names(sort(schpf@feature.loadings[,i], decreasing = T))) %>% 
      mutate(rank = row_number()) %>% 
      filter(genes %in% missing_genes)
    
    res <- left_join(res, ind) %>%rename_at("rank", ~i)
  }
  
  # SAVE
  #-----------------------------------------------------------------------------
  data <- CreateSeuratObject(rbind(GetAssayData(data), as(missing_mat, "dgTMatrix"))[genes, ])
  
  SaveLoom(data, filename = file.path(res_path, "data.loom"))
  write_delim(data.frame(Cells(data)), file.path(res_path, "cells.txt"), col_names = F)
  write_delim(data.frame(genes, genes), file.path(res_path, "gene_input.txt"), col_names = F, delim = "\t")
  write_delim(res, file.path(res_path, "missing_gene_importance.txt"), col_names = T, delim = "\t")
  
}
