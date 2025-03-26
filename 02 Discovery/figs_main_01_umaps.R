#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(ggrastr)
library(ggtext)
library(scCustomize)
library(schex)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(viridis)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "metadata")
fig <- "fig_main1_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))

#===============================================================================
#
# PLOTS
#
#===============================================================================
for(i in factors){
  
  genes <- paste0(names(sort(schpf@feature.loadings[,i], decreasing = T)[1:5]), collapse = ", ")
  p <- Plot_Density_Custom(data, features = i, reduction = "schpf.umap")
  p <- ggplot()+
    geom_point_rast(data = p$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5) +
    scale_colour_gradientn(colors = magma(10, direction = 1)) +
    theme_umap()+
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.margin  = margin(0, 0, 0, 0, "in"))
  
  print(genes)
  
  cairo_pdf(file.path(res_path, paste0(fig, "umap_", i, ".pdf")), height = 2.5, width = 2)
  print(p)
  dev.off()

}
