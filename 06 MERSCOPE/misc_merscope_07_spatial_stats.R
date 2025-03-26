#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE: Calculate spatial stats (Moran & Getis-Ord).
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib_spat")

library(gstat)
library(Seurat)
library(sf)
library(tidyverse)
library(Voyager)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/spatial_stats")

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

#===============================================================================
# 
# CALCULATE SPATIAL STATISTICS
#
#===============================================================================
for(i in 1:12){
  
  mglia <- subset(data, subset = cell_type == "Microglia" & orig.ident == run_metadata$batch_code[i])
  mglia <- DietSeurat(mglia, assays = "RNA")
  
  scores <- t(mglia@meta.data %>% dplyr::select(scHPF_1:scHPF_26))
  
  mglia@assays[["scHPF_scores"]] <- CreateAssayObject(scores, key = "scHPF")
  mglia@assays[["scHPF_scores"]]$data <- scores
  DefaultAssay(mglia) <- "scHPF_scores"
  
  spatial <- coords[colnames(mglia),]
  colnames(spatial) <- paste0("spatial_", 1:2)
  mglia@reductions[["spatial"]] <- CreateDimReducObject(embeddings = spatial, key = "spatial")
  
  for(niche in c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB")){
    keep_cells <- rownames(mglia@meta.data %>% filter(layer_niches %in% niche))
    
    if(length(keep_cells) > 100){
      tmp <- subset(mglia, subset = cell_names %in% keep_cells)
      colData <- as.data.frame(coords[colnames(tmp),1:2])
      meta <- cbind(colData, tmp@meta.data)
      counts <- t(tmp@meta.data %>% dplyr::select(scHPF_1:scHPF_26))
      sfe <- SpatialFeatureExperiment(list(counts = counts, logcounts= counts), 
                                      colData = meta,
                                      spatialCoords = colData,
                                      spatialCoordsNames = c("x", "y"))
      
      colGraph(sfe, "graph") <- findSpatialNeighbors(sfe, k = 5, method = "knearneigh", type = "spatialCoords")
      
      sfe <- runUnivariate(sfe, type = "localG_perm", features = c("scHPF_20", "scHPF_26"), colGraphName = "graph", include_self = TRUE)
      sfe <- runUnivariate(sfe, type = "localmoran_perm", features = c("scHPF_20", "scHPF_26"), colGraphName = "graph")
      
      saveRDS(sfe, file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_mglia_sfe.rds")))
      saveRDS(tmp, file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_mglia_seurat.rds")))
      rm(list = c("tmp", "colData", "meta", "counts", "sfe"))
      
    }
  }
  
  rm(list = c("mglia", "scores", "spatial")); gc()
  
}