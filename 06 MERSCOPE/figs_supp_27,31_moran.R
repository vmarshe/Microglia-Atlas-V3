#!/usr/bin/Rscript

#===============================================================================
# 
# MICROGLIA ATLAS V3
# Figure S27.	Global (tissue-wide) spatial autocorrelation of IFN-I Response (scHPF_20) expression.
# Figure S31.	Global (tissue-wide) spatial autocorrelation of GPNMBHigh (scHPF_26) expression.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib_spat")

library(ggrastr)
library(ggtext)
library(gstat)
library(Nebulosa)
library(patchwork)
library(Seurat); options(Seurat.object.assay.version = "v5")
library(SeuratObject)
library(sf)
library(tidyverse)
library(Voyager)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/2025-02-14_cellular_niches")
figs <- c(20, 26)
names(figs) <- c("figs_s27_", "figs_s31_")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

#===============================================================================
# 
# CALCULATE SPATIAL STATISTICS
#
#===============================================================================
stats <- data.frame()
for(factor in c(20, 26)){
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
        
        sfe <- colDataUnivariate(sfe, type = "moran.mc", features = paste0("scHPF_", factor), colGraphName = "graph", include_self = TRUE, nsim = 10000)
        
        stats <- rbind(stats, cbind(orig.ident = run_metadata$batch_code[i], 
                                    niche, factor,
                                    stat = colFeatureData(sfe)[paste0("scHPF_", factor), ]$moran.mc_statistic_sample01, 
                                    p = colFeatureData(sfe)[paste0("scHPF_", factor), ]$moran.mc_p.value_sample01))
        
        rm(list = c("tmp", "colData", "meta", "counts", "sfe"))
        
      }
    }
    rm(list = c("mglia", "scores", "spatial")); gc()
  }
}
saveRDS(stats, file.path(res_path, "moran_scores.rds"))
#===============================================================================
# 
# AVERAGED
#
#===============================================================================
stats %>% 
  group_by(factor, niche) %>% 
  summarise(mean = mean(as.numeric(as.character(stat))),
            sd = sd(as.numeric(as.character(stat))))
#===============================================================================
# 
# PANEL
#
#===============================================================================
layers <- c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")

for(factor_num in c(20, 26)){
  plots <- list()
  for(i in 1:12){
    
    keep <- stats %>% filter(factor == factor_num & orig.ident == run_metadata$batch_code[i])
    
    moran <- formatC(as.numeric(as.character(keep$stat)), digits = 2, format = "g")
    p <- ifelse(as.numeric(as.character(keep$p)) < 0.001, "***", 
                ifelse(as.numeric(as.character(keep$p)) < 0.01, "**", 
                       ifelse(as.numeric(as.character(keep$p)) < 0.05, "*", "")))
    
    lab <- paste0(paste0("Moran's <i>I</i> (WM) = ", moran[1], p[1]), "<br>", 
                  paste0("Moran's <i>I</i> (GM L2-6) = ", moran[2], p[2]), "<br>", 
                  if(!is.na( moran[3])){paste0("Moran's <i>I</i> (GM L1) = ", moran[3], p[3])} else {""})
    
    mglia <- subset(data, subset = cell_type == "Microglia" & orig.ident == run_metadata$batch_code[i])
    mglia@meta.data <- mglia@meta.data %>% 
      left_join(data.frame(cell_names = mglia@images$merged@boundaries$centroids@cells, mglia@images$merged@boundaries$centroids@coords))
    rownames(mglia@meta.data) <- mglia@meta.data$cell_names
    scores <- t(mglia@meta.data %>% dplyr::select(scHPF_1:scHPF_26))
    
    mglia@assays[["scHPF_scores"]] <- CreateAssayObject(scores, key = "scHPF")
    mglia@assays[["scHPF_scores"]]$data <- scores
    DefaultAssay(mglia) <- "scHPF_scores"
    
    spatial <- coords[colnames(mglia),]
    colnames(spatial) <- paste0("spatial_", 1:2)
    mglia@reductions[["spatial"]] <- CreateDimReducObject(embeddings = spatial, key = "spatial")
    
    p <- Nebulosa::plot_density(mglia, 
                                features = paste0("scHPF-", factor_num), 
                                slot = "data", 
                                reduction = "spatial", 
                                pal = "magma", 
                                size = 0.2)+
      labs(title = run_metadata$batch_name[i], x = lab)+
      theme_umap()+
      theme(legend.key.height = grid::unit(0.15, "in"),
            legend.key.width = grid::unit(0.1, "in"),
            legend.text = element_text(size = 8, margin = margin(l=0.05, unit='in'), family = "Helvetica"),
            legend.title = element_text(size = 8, margin = margin(b=0.05, unit='in')),
            legend.margin = margin(l=-0.2, unit='in'),
            plot.title = element_markdown(size = 10, margin = margin(b = -0.05, unit='in'), face = "bold"),
            axis.title.x = element_markdown(size = 8, hjust = 0), 
            plot.margin = margin(t=0.05, b=0.05, r=0.05, l=0.05, unit = "in"))
    
    plots <- c(plots, list(p))
  }
  
  p <- wrap_plots(plots) + plot_layout(nrow = 4)
  ggsave(file.path(res_path, paste0(names(figs)[figs %in% factor_num], "panelA.png")),
         height = 8, width = 6, dpi = 600)
}

