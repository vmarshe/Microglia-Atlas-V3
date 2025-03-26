#!/usr/bin/Rscript

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(ggrastr)
library(ggridges)
library(ggtext)
library(patchwork)
library(scCustomize)
library(schex)
library(Seurat)
library(tidyverse)
library(viridis)

# PATHS
#-------------------------------------------------------------------------------
dataset <- "Dolan_FinMg"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp13_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))

# REFERENCE DATA
#-------------------------------------------------------------------------------
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))
cell_embeddings <- data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names")
data@meta.data <- data@meta.data %>% left_join(cell_embeddings)
rownames(data@meta.data) <- data@meta.data$cell_names

#===============================================================================
#
# PANEL A
#
#===============================================================================
hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  labs(title = "Idiopathic HC")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(l = -0.2, unit ='in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 3.5, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
plots <- list()
for(i in factors){
  
  tmp <- data@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i]))
  
  p <- ggplot(tmp, aes(x = factor)) + 
    geom_density(fill = proj_cols("light blue"))+
    scale_x_continuous(trans = "log2", labels = function(x) formatC(x, digits = 1,  format = "g")) + 
    geom_vline(xintercept = 0.01, lty = 2) + 
    facet_wrap(~facet)+
    labs(x = "log(Cell score)", y = "Density")+
    theme_proj() +
    theme(axis.text = element_text(size = 10, family = "Helvetica"),
          axis.title.y = element_text(size = 14, family = "Helvetica"),
          axis.title.x =  element_blank(),
          legend.position = "none",
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          panel.grid.major = element_blank())
  
  
  if(grepl("-", names(factors)[factors %in% i])){
    p <- p + theme(strip.text = element_markdown(size = 10, family = "Helvetica", margin = margin(b = 0.05, t = 0.05, unit = "in")))
  } else {
    p <- p + theme(strip.text = element_text(size = 10, family = "Helvetica", margin = margin(b = 0.05, t = 0.05, unit = "in")))
  }
  
  plots <- c(plots, list(p))
}

plots[[21]] <- plots[[21]] + theme(axis.title.x = element_text(size = 14, family = "Helvetica"))
p <- wrap_plots(plots) + plot_layout(guides = "collect", nrow = 4, axes = "collect_y", axis_titles = "collect_y")

cairo_pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 7, width = 14)
print(p)
dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
prop <- data.frame()
for(i in factors){
  
  tmp <- data@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

p <- prop %>%
  mutate(prop = as.numeric(as.character(prop)), 
         factor = factor(i, levels = rev(factors), labels = rev(names(factors))),
         factor = gsub('-high', paste0("<sup>high</sup>"), factor),
         lab = paste0(round(as.numeric(prop)*100, 1), "%")) %>%
  arrange(factor) %>%
  ggplot(aes(x = prop, y = factor)) +
  geom_col(fill = proj_cols("yellow")) +
  geom_text(aes(label = lab), hjust = 0, family = "Helvetica", size = 3)+
  geom_vline(xintercept = 0, linewidth = 1)+
  scale_x_continuous(expand = c(0,0),limits = c(0, 1.5), breaks = c(0, 1), labels = c(0, 1))+
  labs(x = "Proportion of\nexpressing cells", y = "Expression Programs")+
  theme_proj()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        plot.margin = margin(t=0.1, r=-0.05, b=0.1, l=0.1, unit = "in"),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 11, margin = margin(b=0.05, 0, 0, 0, unit = "in"), family = "Helvetica"),
        axis.title.y = element_text(size = 12, family = "Helvetica"),
        axis.text.y = element_markdown(size = 11, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5))

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 5, width = 4)
print(p)
dev.off()
#===============================================================================
#
# PANEL D
#
#===============================================================================
xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

keep_factors <- prop %>% filter(prop >= 0.25) %>% pull(i)
plots <- list()
for(i in keep_factors){
  
  p <- Plot_Density_Custom(data, features = i, reduction = "schpf.umap", pt.size = 0.1)
  p$data$facet <- gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i])
  p <- ggplot()+
    geom_point_rast(data = p$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5) +
    facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(limits = c(ymin, ymax))+
    scale_color_gradientn(colors = magma(12))+
    theme_umap()+
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.title = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          strip.placement = "inside", 
          strip.background = element_rect(fill = NA, color = NA, linewidth = NA),
          plot.margin  = margin(0, 0, 0, 0, "in"))
  
  if(grepl("-", names(factors)[factors %in% i])){
    p <- p + theme(strip.text = element_markdown(size = 10, family = "Helvetica", margin = margin(b = 0.05, t = 0.05, unit = "in")))
  } else {
    p <- p + theme(strip.text = element_text(size = 10, family = "Helvetica", margin = margin(b = 0.05, t = 0.05, unit = "in")))
  }
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 4)

cairo_pdf(file.path(res_path,  paste0(fig, "panel_D.pdf")), height = 7, width = 14)
print(p)
dev.off()
