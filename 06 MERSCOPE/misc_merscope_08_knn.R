#!/usr/bin/Rscript

#===============================================================================
# 
# MERSCOPE: Identification k-nearest neighbors.
#
#===============================================================================
# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(Seurat)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/spatial_stats")

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

#===============================================================================
# 
# DEFINE KNN
#
#===============================================================================
for(i in 1:12){
  for(niche in c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB")){
    mglia_file <- file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_mglia_sfe.rds"))
    if(file.exists(mglia_file)){
      keep_cells <- rownames(data@meta.data %>% filter(orig.ident == run_metadata$batch_code[i] & layer_niches %in% niche))
      subset <- subset(data, subset = cell_names %in% keep_cells)
      subset@images[["merged"]] <- CreateFOV(as.data.frame(coords[subset$cell_names,1:2]), type = "centroids", name = "merged")
      
      tmp <- GetTissueCoordinates(subset[["merged"]], which = "centroids")
      rownames(tmp) <- tmp[["cell"]]
      tmp <- as.matrix(tmp[ , c("x", "y")])
      
      set.seed(1234)
      nn <- FindNeighbors(tmp, k.param = 30, compute.SNN = FALSE, return.neighbor = T)
      
      saveRDS(nn, file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_nn.rds")))
      rm(list = c("keep_cells", "subset", "tmp", "nn")); gc()
    }
  }
}
