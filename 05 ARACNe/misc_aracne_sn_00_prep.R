#!/usr/bin/Rscript
#===============================================================================
# 
# VALIDATION ARACNE MODEL
# This script create metacells in ROS-MAP CUIMC1 dataset
# to run ARACNe for the validation model
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(data.table)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sn")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/aracne.R"))

# DATA
#-------------------------------------------------------------------------------
hfile <- Connect(filename = file.path(home_path, "scHPF_snRNAseq/microglia.h5Seurat"), mode = "r", force = T)
test <- as.Seurat(x = hfile)
mat <- as.matrix(GetAssayData(test, assay = "RNA", slot = "count"))

# MAKE METACELLS
#-------------------------------------------------------------------------------
mat <- CreateSeuratObject(mat)
mat$cell_names <- rownames(mat@meta.data)
snuc_scHPF_projection <- readRDS(file.path(home_path, "scHPF_snRNAseq/snuc_scHPF_projection.rds"))

keep_cells <- intersect(colnames(mat), colnames(snuc_scHPF_projection))

mat <- subset(mat, subset = cell_names %in% keep_cells)

mat[["scHPF"]] <- CreateDimReducObject(embeddings = snuc_scHPF_projection@reductions$scHPF@cell.embeddings[keep_cells,],
                                       loadings = snuc_scHPF_projection@reductions$scHPF@feature.loadings,
                                       assay = "RNA")

meta_mat <- make_metacells(data = mat, 
                           dims = 1:23, 
                           reduction = "scHPF", 
                           k.param = 50,
                           features = rownames(mat@reductions$scHPF@feature.loadings))

meta_mat <- log2(meta_mat + 1)
meta_mat <- cbind(gene = rownames(meta_mat), meta_mat)

fwrite(meta_mat, 
       file = file.path(res_path, "metacell_matrix.txt"),
       sep = "\t", 
       quote = F, 
       row.names = F, 
       col.names = T)