#!/usr/bin/Rscript

#===============================================================================
# 
# MERSCOPE: 
# Figure S28.	In situ distribution of IRHigh local cellular neighborhoods. 
# Figure S32.	In situ distribution of GPNMBHigh local cellular neighborhoods. 
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib_spat")

library(circlize)
library(ComplexHeatmap)
library(ggrastr)
library(ggsignif)
library(ggtext)
library(gstat)
library(mclust)
library(patchwork)
library(Seurat); options(Seurat.object.assay.version = "v5")
library(SeuratObject)
library(sf)
library(tidyverse)
library(viridisLite)
library(Voyager)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/2025-02-14_cellular_niches")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/load_merscope.R"))
adjacent_cells <- readRDS(file.path(res_path, "adjacent_cells.rds"))

#===============================================================================
# 
# PLOTs
#
#===============================================================================
cols <- proj_cols("pink", 'yellow', "teal")
names(cols) <- c("WM", "GM\n(L1)", "GM\n(L2-6)")

layers <- c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")

for(factor in c(20, 26)){
  plots <- list()
  for(i in 1:12){
    
    tmp <- subset(data, subset = orig.ident == run_metadata$batch_code[i])
    tmp@meta.data <- tmp@meta.data %>% 
      left_join(data.frame(cell_names = tmp@images$merged@boundaries$centroids@cells, tmp@images$merged@boundaries$centroids@coords))
    
    p0 <- coords_tmp %>%
      mutate(layer_niches = factor(layer_niches, 
                                   levels = c("GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB", "WM"),
                                   labels = c("GM\n(L2-6)", "GM\n(L1)",  "WM"))) %>%
      arrange(desc(layer_niches)) %>%
      ggplot(aes(x = x, y = y, color = layer_niches))+
      geom_point_rast(size = 0.1, scale = 0.75)+
      labs(y = run_metadata$batch_name[i])+
      scale_color_manual(values = cols)+
      guides(color= guide_legend(override.aes = aes(size = 4)))+
      theme_umap()+
      theme(legend.text = element_text(size = 8, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
            legend.title = element_blank(),
            legend.margin = margin(l=-0.2, unit='in'),
            legend.key = element_blank(),
            plot.margin = margin(0, 0, 0, 0),
            axis.title.y = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), angle = 90, size = 12, face = "bold", family = "Helvetica"))
    
    for(j in layers){
      file <- file.path(res_path, paste0(run_metadata$batch_code[i], "_", j, "_mglia_sfe.rds"))
      if(file.exists(file)){
        tmp <- readRDS(file)
        assign(paste0(run_metadata$batch_code[i], "_", j), tmp)
      }
    }
    
    coords_tmp <- data.frame()
    for(j in layers){
      if(exists(paste0(run_metadata$batch_code[i], "_", j))){
        tmp <- eval(parse(text = paste0(run_metadata$batch_code[i], "_", j)))
        coords_tmp <- rbind(coords_tmp, cbind(tmp@int_colData$localResults$localG_perm[[paste0("scHPF_", factor)]], tmp@int_colData$spatialCoords,  tmp@colData))
        rm(tmp)
      }
    }
    
    p1 <- coords_tmp %>%
      arrange(!!sym(paste0("scHPF_", factor))) %>%
      ggplot(aes(x = x, y = y, color = !!sym(paste0("scHPF_", factor))))+
      geom_point_rast(size = 0.1, scale = 0.75)+
      scale_color_gradient(low = "lightgrey", high = "blue")+
      theme_umap()+
      theme(legend.key.height = grid::unit(0.1, "in"),
            legend.key.width = grid::unit(0.1, "in"),
            legend.text = element_text(size = 8, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
            legend.title = element_blank(),
            legend.margin = margin(l=-0.2, unit='in'),
            plot.margin = margin(0, 0, 0, 0))
    
    mglia <- subset(data, subset = cell_type == "Microglia" & orig.ident == run_metadata$batch_code[i])
    
    scores <- t(mglia@meta.data %>% dplyr::select(scHPF_1:scHPF_26))
    
    mglia@assays[["scHPF_scores"]] <- CreateAssayObject(scores, key = "scHPF")
    mglia@assays[["scHPF_scores"]]$data <- scores
    DefaultAssay(mglia) <- "scHPF_scores"
    
    spatial <- coords[colnames(mglia),]
    colnames(spatial) <- paste0("spatial_", 1:2)
    mglia@reductions[["spatial"]] <- CreateDimReducObject(embeddings = spatial, key = "spatial")
    
    p2 <- Nebulosa::plot_density(mglia, 
                                 features = paste0("scHPF-", factor), 
                                 slot = "data", 
                                 reduction = "spatial", 
                                 pal = "magma", 
                                 size = 0.2)+
      theme_umap()+
      theme(legend.key.height = grid::unit(0.1, "in"),
            legend.key.width = grid::unit(0.1, "in"),
            legend.text = element_text(size = 8, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
            legend.title = element_blank(),
            legend.margin = margin(l=-0.2, unit='in'),
            plot.title = element_blank(),
            plot.margin = margin(0, 0, 0, 0))
    
    plot_data <- coords_tmp %>%
      rename_at(all_of("-log10p_adj Sim"), ~"p") %>%
      rename_at(all_of("localG"), ~"G") %>%
      mutate(col = ifelse(G > 0, p, -p)) %>%
      arrange(abs(col))
    
    p3 <- ggplot(plot_data, aes(x = x, y = y, color = col))+
      geom_point_rast(size = 0.1, scale = 0.75)+
      scale_color_gradient2(name = "Gi*Adj p", low = proj_cols("blue"), high = proj_cols("red"), midpoint = 0, mid = "#F0F0F0",
                            breaks = c(-1.5, 0, 1.5))+
      theme_umap()+
      theme(legend.key.height = grid::unit(0.1, "in"),
            legend.key.width = grid::unit(0.1, "in"),
            legend.text = element_text(size = 8, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
            legend.title = element_blank(),
            legend.margin = margin(l=-0.2, unit='in'),
            plot.margin = margin(0, 0, 0, 0))
    
    tmp <- subset(data, subset = orig.ident == run_metadata$batch_code[i])
    
    neighbors <- unlist(str_split(adjacent_cells %>% filter(factor == 20) %>% pull(neighbors), "[,]"))
    neighbors <- names(table(neighbors))[which(table(neighbors) > 0)]
    
    tmp@meta.data <- tmp@meta.data %>% 
      mutate(lab = ifelse(cell_names %in% neighbors, 1, 0)) %>%
      left_join(data.frame(cell_names = tmp@images$merged@boundaries$centroids@cells, tmp@images$merged@boundaries$centroids@coords))
    
    p4 <- ggplot() +
      geom_point_rast(data = coords_tmp %>% filter(layer_niches == "WM"), aes(x, y),
                      color = proj_cols("pink"), size = 0.1, alpha = 0.1, shape = ".")+
      geom_point_rast(data = coords_tmp %>% filter(layer_niches == "GM.AHNAK.FOS.ACTB"), aes(x, y),
                      color = proj_cols("yellow"), size = 0.1, alpha = 0.1, shape = ".")+
      geom_point_rast(data = coords_tmp %>% filter(layer_niches == "GM.SLC17A7.GAS7.HSP90AB1"), aes(x, y),
                      color = proj_cols("teal"), size = 0.1, alpha = 0.1, shape = ".")+
      geom_point_rast(data = tmp@meta.data %>% filter(lab == 1 & cell_type!="Microglia"), aes(x, y), color = proj_cols("blue"), size = 0.1, scale = 0.5) +
      geom_point_rast(data = coords_tmp %>% filter(!!sym(paste0("scHPF_", factor, "_class")) == paste0("Mic.", factor, ".high")), aes(x, y),color = "black", size = 0.1, scale = 0.5) +
      theme_umap()+
      theme(legend.position = "none",
            plot.margin = margin(0, 0, 0, 0))
    
    plots <- c(plots, list(wrap_plots(p0, p1, p2, p4, p3) + plot_layout(nrow = 1)))
  }
  plots[[1]][[1]] <- plots[[1]][[1]] + labs(title = "Tissue Layers")+
    theme(plot.title = element_textbox_simple(family = "Helvetica", face = "bold", size = 10,halign = 0.5,
                                              box.color = "#F0F0F0", fill = "#F0F0F0", 
                                              padding = margin(2.5, 2.5, 2.5, 2.5)))
  
  plots[[1]][[2]] <- plots[[1]][[2]] + labs(title = "Raw Score")+
    theme(plot.title = element_textbox_simple(family = "Helvetica", face = "bold", size = 10,halign = 0.5,
                                              box.color = "#F0F0F0", fill = "#F0F0F0", 
                                              padding = margin(2.5, 2.5, 2.5, 2.5)))
  
  plots[[1]][[3]] <- plots[[1]][[3]] + labs(title = "Density")+
    theme(plot.title = element_textbox_simple(family = "Helvetica", face = "bold", size = 10, halign = 0.5,
                                              box.color = "#F0F0F0", fill = "#F0F0F0", 
                                              padding = margin(2.5, 2.5, 2.5, 2.5)))
  
  plots[[1]][[4]] <- plots[[1]][[4]] + labs(title = "Niches")+
    theme(plot.title = element_textbox_simple(family = "Helvetica", face = "bold", size = 10,halign = 0.5,
                                              box.color = "#F0F0F0", fill = "#F0F0F0", 
                                              padding = margin(2.5, 2.5, 2.5, 2.5)))
  
  plots[[1]][[5]] <- plots[[1]][[5]] + labs(title = "Getis-Ord G*")+
    theme(plot.title = element_textbox_simple(family = "Helvetica", face = "bold", size = 10,halign = 0.5,
                                              box.color = "#F0F0F0", fill = "#F0F0F0", 
                                              padding = margin(2.5, 2.5, 2.5, 2.5)))
  
  p <- wrap_plots(plots) + plot_layout(ncol = 1)
  ggsave(file.path(res_path, paste0("spatial_stats_", factor, ".png")), height = 14, width = 9, dpi = 600)
}
