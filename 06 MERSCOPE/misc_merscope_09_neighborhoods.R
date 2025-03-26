#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE: Identification of cells adjacent to high-expression microglia.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib")

library(Seurat)
library(SeuratObject)
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
# NUMBER OF ADJACENT MICROGLIA AND DISTANCE
#
#===============================================================================
res <- data.frame()
for(i in 1:12){
  for(niche in c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB")){
    
    mglia_file <- file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_mglia_sfe.rds"))
    
    if(file.exists(mglia_file)){
      
      keep_cells <- rownames(data@meta.data %>% filter(layer_niches == niche & orig.ident == run_metadata$batch_code[i]))
      subset <- subset(data, subset = cell_names %in% keep_cells)
      subset@images[["merged"]] <- CreateFOV(as.data.frame(coords[subset$cell_names,1:2]), type = "centroids", name = "merged")
      subset@meta.data <- subset@meta.data %>%
        left_join(cbind.data.frame(subset@images$merged@boundaries$centroids@coords , cell_names = subset@images$merged@boundaries$centroids@cells))
      rownames(subset@meta.data) <- subset@meta.data$cell_names
      
      nn <- readRDS(file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_nn.rds")))
      nn.dist <- nn@nn.dist
      rownames(nn.dist) <- nn@cell.names
      
      nn.idx <-  nn@nn.idx
      rownames(nn.idx) <- nn@cell.names
      
      for(factor in c(20, 26)){
        
        subset@meta.data <- subset@meta.data %>%
          mutate(lab = ifelse(cell_type == "Microglia", as.character(!!sym(paste0("scHPF_", factor, "_class"))), as.character(cell_type)))
        
        mglia_high <- rownames(subset@meta.data %>% filter(lab == paste0("Mic.", factor, ".high")))
        dist <- lapply(mglia_high,
                       function(x) {
                         neighbors <- rownames(nn.idx)[nn.idx[x,]][-1]
                         dist <- nn.dist[x,-1]
                         names(dist) <- neighbors
                         neighbors <- subset@meta.data %>% filter(cell_names %in% neighbors & lab == paste0("Mic.", factor, ".high") & orig.ident == subset@meta.data[x, "orig.ident"])
                         neighbors <- rownames(neighbors)
                         counts <- length(neighbors)
                         dist <- if(length(neighbors) > 0){mean(dist[neighbors])} else {NA}
                         return(cbind(niche, 
                                      factor, 
                                      cell_name = x, 
                                      x = subset@meta.data[x,"x"], 
                                      y = subset@meta.data[x,"y"],
                                      stat = "High", 
                                      counts, dist,
                                      neighbors = paste0(neighbors, collapse=",")))
                       }
        )
        
        res <- rbind(res, do.call(rbind.data.frame, dist))
        

        mglia_low <- rownames(subset@meta.data %>% filter(lab == paste0("Mic.", factor, ".low")))
        dist <- lapply(mglia_low,
                       function(x) {
                         neighbors <- rownames(nn.idx)[nn.idx[x,]][-1]
                         dist <- nn.dist[x,-1]
                         names(dist) <- neighbors
                         neighbors <- subset@meta.data %>% filter(cell_names %in% neighbors & lab == paste0("Mic.", factor, ".high") & orig.ident == subset@meta.data[x, "orig.ident"])
                         neighbors <- rownames(neighbors)
                         counts <- length(neighbors)
                         dist <- if(length(neighbors) > 0){mean(dist[neighbors])} else {NA}
                         return(cbind(niche, 
                                      factor, 
                                      cell_name = x, 
                                      x = subset@meta.data[x,"x"], 
                                      y = subset@meta.data[x,"y"],
                                      stat = "Low", 
                                      counts, 
                                      dist,
                                      neighbors = paste0(neighbors, collapse=",")))
                       }
        )
        res <- rbind(res, do.call(rbind.data.frame, dist))
        
      }
    }
  }
}
saveRDS(res, file.path(res_path, "adjacent_counts.rds"))
#===============================================================================
#
# IDENTIFY ADJACENT CELLS
#
#===============================================================================
res <- data.frame()

for(i in 1:12){
  for(niche in c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB")){
    
    mglia_file <- file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_mglia_sfe.rds"))
    if(file.exists(mglia_file)){
      
      keep_cells <- rownames(data@meta.data %>% filter(layer_niches == niche & orig.ident == run_metadata$batch_code[i]))
      subset <- subset(data, subset = cell_names %in% keep_cells)
      
      nn <- readRDS(file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_nn.rds")))
      nn.dist <- nn@nn.dist
      rownames(nn.dist) <- nn@cell.names
      nn.idx <- nn@nn.idx
      rownames(nn.idx) <- nn@cell.names
      
      for(factor in c(20, 26)){
        
        subset@meta.data <- subset@meta.data %>%
          mutate(lab = ifelse(cell_type == "Microglia", as.character(!!sym(paste0("scHPF_", factor, "_class"))), as.character(cell_type)))
        
        mglia_high <- rownames(subset@meta.data %>% filter(lab == paste0("Mic.", factor, ".high")))
        for(j in c("Adaptive Immune", "Astrocytes", "Endothelial", "Excitatory", "Inhibitory", "Oligodendrocytes", "OPCs", "Pericytes")){
          neighbors <- lapply(mglia_high,
                              function(x) {
                                neighbors <- rownames(nn.idx)[nn.idx[x,]][-1]
                                dist <- nn.dist[x,-1]
                                names(dist) <- neighbors
                                neighbors <- subset@meta.data %>% filter(cell_names %in% neighbors & lab == j & orig.ident == subset@meta.data[x, "orig.ident"])
                                neighbors <- rownames(neighbors)
                                
                                if(length(neighbors) > 0){
                                  return(cbind(orig.ident = run_metadata$batch_code[i],
                                               niche, factor,
                                               cell_name = x,
                                               x = subset@meta.data[x, "x"],
                                               y = subset@meta.data[x, "y"],
                                               cell_type = j,
                                               neighbors = paste0(neighbors, collapse=",")))
                                }
                              })
          res <- rbind(res, do.call(rbind.data.frame, neighbors))
        }
      }
      rm(list = c("subset", "nn", "nn.dist", "nn.idx", "keep_cells")); gc()
    }
  }
}

saveRDS(res, file.path(res_path, "adjacent_cells.rds"))