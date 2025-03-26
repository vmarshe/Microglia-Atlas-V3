#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(circlize)
library(colorspace)
library(ComplexHeatmap)
library(ggrepel)
library(ggsignif)
library(ggtext)
library(logger)
library(patchwork)
library(readxl)
library(rms)
library(scCustomize)
library(schex)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
fig <- "fig_main2_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_helpers.R"))

# REFERENCE DATA
#-------------------------------------------------------------------------------
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

#===============================================================================
#
# PANELS B, C: DRUG TREATED GBM
#
#===============================================================================
dataset <- "GBM_slice"
res_path <- file.path(home_path, paste0("proj_", dataset))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names

# PANEL B
#-------------------------------------------------------------------------------
hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(0, 0, 0, -0.1, 'in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, paste0(fig, "panel_B1.pdf")), height = 3, width = 3.5)
print(p)
dev.off()

plots <- list()
for(i in c("Ctrl", "Topotecan")){
  
  keep_id <- data@meta.data %>% filter(Drug_Tx == "Topotecan") %>% group_by(TB_ID, Drug_Tx) %>% pull(TB_ID)
  keep_id <- unique(keep_id)
  
  tmp <- subset(data, subset = TB_ID %in% keep_id)
  hex <- make_hexbin(subset(tmp, subset = Drug_Tx == i), nbins = 50, dimension_reduction = "schpf.umap")
  
  p <- ggplot() + 
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    ggtitle(i)+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(l=0.05, unit = 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, 'in'),
          plot.title = element_text(margin = margin(0, 0, 0.01, 0, 'in'), size = 14, 
                                    color = "black",
                                    face = "plain"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0.05, 0.05, 0, 0.05))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(ncol = 2)

pdf(file.path(res_path,  paste0(fig, "panel_B2.pdf")), height = 2.5, width = 5)
print(p)
dev.off()

# PANEL C
#-------------------------------------------------------------------------------
plots <- list()
for(i in c("scHPF_11", "scHPF_14", "scHPF_15", "scHPF_21")){
  
  thresh <- median(data@meta.data[,i]) + 2*mad(data@meta.data[,i])
  
  tmp <- data@meta.data %>% filter(Drug_Tx %in% c("Topotecan", "Ctrl"))
  
  keep_id <- tmp %>% 
    group_by(TB_ID, Drug_Tx) %>% 
    summarise(n = n()) %>% 
    group_by(TB_ID) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>% 
    pull(TB_ID)
  
  tmp = tmp %>%
    filter(TB_ID %in% keep_id) %>% 
    rename_at(i, ~"factor") %>%
    mutate(class = ifelse(factor > thresh, "High", "Low"),
           class = factor(class, levels = c("Low", "High"))) %>%
    filter(TB_ID %in% keep_id)
  
  plot_data <- tmp %>% 
    group_by(Drug_Tx, class) %>%
    summarise(n = n()) %>%
    group_by(Drug_Tx) %>% 
    reframe(class, n, prop = n/sum(n)) %>% 
    mutate(lab = paste0(round(prop*100, digits = 0), "%"),
           facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i]))
  
  anno <- chisq.test(tmp$Drug_Tx, tmp$class)$p.value
  anno <- ifelse(anno < 0.001, "***", ifelse(anno < 0.01, "**", ifelse(anno < 0.05, "*", "ns")))
  
  p <- ggplot(plot_data, aes(x = Drug_Tx, y = prop, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label = lab, y = prop + 0.01),
              position = position_dodge(width = 0.9),
              color = 'black',
              vjust = 0,
              size = 3,
              family = "Helvetica") +
    geom_signif(comparisons = list(c("Ctrl", "Topotecan")),
                step_increase = 0.1, 
                tip_length = 0, 
                margin_top = 0.15,
                vjust = 0.5,
                annotation = anno,
                family = "Helvetica", textsize = 4)+
    facet_wrap(~facet, strip.position = "top") +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 1.2),
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) formatC(x, format = "g"))+
    scale_fill_manual(name = "Expression\nstatus",
                      values = unname(proj_cols("light blue", "purple"))) +
    labs(y = "Proportion of cells")+
    theme_proj()+
    theme(panel.grid.major = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in"), hjust = 0),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          plot.margin = margin(t=0.1, r=0, b=0.05, l=0.05, "in"),
          legend.margin = margin(l = -0.1, unit = "in"),
          strip.text = element_text(size = 14, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")),
          strip.background = element_blank(),
          strip.clip = "off",
          plot.background = element_rect(fill = "transparent"))
  
  if(i %in% c("scHPF_15", "scHPF_11")){
    p <- p + 
      theme(strip.text = element_markdown(size = 14, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")))
  }
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) +
  plot_layout(guides = "collect", axis_title = "collect", nrow = 1)

pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 3, width = 8, bg = "transparent")
print(p)
dev.off()

#===============================================================================
#
# PANELS D-G: MURINE XENOGRAFT
#
#===============================================================================
dataset <- "xeno"
res_path <- file.path(home_path, paste0("proj_", dataset))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  as.data.frame() %>% 
  mutate(treatment = ifelse(grepl("WT", cell_names), 0, 1),
         treatment_name = ifelse(grepl("WT", cell_names), "WT", "5X"),
         sex = ifelse(grepl("Female", cell_names), 1, 0)) %>%
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names


# PANEL D
#-------------------------------------------------------------------------------
hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(0, 0, 0, -0.1, 'in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, "analysis", paste0(fig, "panel_D.pdf")), height = 3, width = 3.5)
print(p)
dev.off()

# PANEL E
#-------------------------------------------------------------------------------
wt <- subset(data, subset = treatment == 0)
non_wt <- subset(data, subset = treatment == 1)

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

plots <- list()
for(i in c(20, 26)){
  
  ind <- which(factors %in% paste0("scHPF_", i))
  
  p <- Plot_Density_Custom(wt, features = paste0("scHPF_", i))
  p$data$facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[ind])# names(factors)[ind]
  p <- p + 
    facet_wrap(~facet, strip.position = "left")+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    theme_umap()+
    theme(strip.placement = "inside",
          strip.text = element_markdown(size = 14),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0, 0, 0, 0, unit = "in"),
          legend.key.height = unit(0.25, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(t=0, r=0, b=0, l=0.05, 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t=0, r=0, b=0, l=-0.2, 'in'))
  
  plots <- c(plots, list(p))
  
  p <- Plot_Density_Custom(non_wt, features = paste0("scHPF_", i)) + 
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    theme_umap()+
    theme(plot.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "in"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          legend.key.height = unit(0.25, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(t=0, r=0, b=0, l=0.05, 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t=0, r=0, b=0, l=-0.2, 'in'))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(nrow = 2)

pdf(file.path(res_path, "analysis", paste0(fig, "panel_E.pdf")),  height = 4.5, width = 5.5)
print(p)
dev.off()

# PANEL F
#-------------------------------------------------------------------------------
plots <- list()
for(i in c("scHPF_20", "scHPF_26")){
  
  thresh <- median(data@meta.data[,i]) + 2*mad(data@meta.data[,i])
  
  tmp <- data@meta.data %>%
    rename_at(i, ~"factor") %>%
    mutate(class = ifelse(factor > thresh, "High", "Low"),
           class = factor(class, levels = c("Low", "High")), 
           treatment_name = factor(treatment_name, levels = c("WT", "5X")))
  
  plot_data <- tmp %>% 
    group_by(treatment_name, class) %>%
    summarise(n = n()) %>%
    group_by(treatment_name) %>% 
    reframe(class, n, prop = n/sum(n)) %>% 
    mutate(lab = paste0(round(prop*100, digits = 0), "%"),
           facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i]))
  
  anno <- chisq.test(tmp$treatment_name, tmp$class)$p.value
  anno <- ifelse(anno < 0.001, "***", ifelse(anno < 0.01, "**", ifelse(anno < 0.05, "*", "ns")))
  
  p <- ggplot(plot_data, aes(x = treatment_name, y = prop, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label = lab, y = prop + 0.01),
              position = position_dodge(width = 0.9),
              color = 'black',
              vjust = 0,
              size = 3,
              family = "Helvetica") +
    geom_signif(comparisons = list(c("WT", "5X")),
                step_increase = 0.1, 
                tip_length = 0, 
                margin_top = 0.1,
                vjust = 0,
                annotation = anno,
                family = "Helvetica", textsize = 4)+
    facet_wrap(~facet, strip.position = "top") +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1.1), 
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) formatC(x, format = "g"))+
    scale_fill_manual(name = "Expression\nstatus",
                      values = unname(proj_cols("light blue", "purple"))) +
    labs(y = "Proportion of cells")+
    theme_proj()+
    theme(panel.grid.major = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in"), hjust = 0),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          plot.margin = margin(t=0.1, r=0, b=0.05, l=0.05, "in"),
          legend.margin = margin(l = -0.1, unit = "in"),
          strip.text = element_markdown(size = 14, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")),
          strip.background = element_blank(),
          strip.clip = "off",
          plot.background = element_rect(fill = "transparent"))
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) +
  plot_layout(guides = "collect", axis_title = "collect", nrow = 1)

pdf(file.path(res_path, "analysis", paste0(fig,  "panel_F.pdf")), height = 3, width = 5)
print(p)
dev.off()

# PANEL G
#-------------------------------------------------------------------------------
Idents(data) <- "treatment_name"
res <- FindMarkers(data, 
                   ident.1 = "5X", 
                   ident.2 = "WT", 
                   assay = "RNA",
                   test.use = "MAST", 
                   latent.vars = c("nCount_RNA", "sex"), 
                   logfc.threshold = 0.1,
                   max.cells.per.ident = 5000, 
                   random.seed = 1234)

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("1", "2", "3")

p <- res %>%
  rownames_to_column("gene") %>%
  mutate(col = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, "2", "3"),
         col = ifelse(avg_log2FC < 0 & p_val_adj < 0.05, "1", col),
         lab = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, gene, NA),
         log10.p_val_adj = -log10(ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj!=0]), p_val_adj))) %>%
  ggplot(aes(x = avg_log2FC, y =  log10.p_val_adj, color = col)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_point(alpha = 0.6)+
  geom_vline(xintercept = 0,col = "black")+
  scale_color_manual(values = cols,
                     labels = c("Downregulated", "Upregulated", "NS (Adj. p > 0.05)"))+
  guides(color = guide_legend(override.aes = aes(size = 4)))+
  geom_text_repel(aes(label = lab), family = "Helvetica", seed=1234, size = 3, show.legend = F) +
  labs(x = "logFC (5X vs. WT)",
       y = "-log<sub>10</sub>(Adj. <i>p</i>-value)")+
  theme_proj()+
  theme(legend.margin = margin(l = -0.1, unit = "in"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
        panel.border = element_rect(linewidth = 1, color = "black"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"))

pdf(file.path(res_path, "analysis", paste0(fig,  "panel_G.pdf")), height = 3, width = 6)
print(p)
dev.off()

#===============================================================================
#
# PANELS H-K: DOLAN iPSCs
#
#===============================================================================
dataset <- "Dolan_iPSC"
res_path <- file.path(home_path, paste0("proj_", dataset))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))
tmp <- read_xlsx(file.path(res_path, "iPSC_H1_Integration_metadata2024.xlsx")) %>% 
  mutate(cell_names = ifelse(line == "H1", paste0(dataset, "_", bc), paste0(line, "_", bc)))

data@meta.data <- data@meta.data %>% 
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(tmp, by = "cell_names")

rownames(data@meta.data) <- data@meta.data$cell_names

data <- subset(data, subset = line %in% c("iCW50036", "iCW50118", "iCW70347"))

# PANEL H
#-------------------------------------------------------------------------------
hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(0, 0, 0, -0.1, 'in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, "analysis", paste0(fig, "panel_H.pdf")), height = 3, width = 3.5)
print(p)
dev.off()

# PANEL I
#-------------------------------------------------------------------------------
xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

untreated <- subset(data, subset = treatment == "Ctrl")
treated <- subset(data, subset = treatment == "Apop")

plots <- list()
for(i in c(5, 20)){
  
  ind <- which(factors %in% paste0("scHPF_", i))
  
  p <- Plot_Density_Custom(untreated, features = paste0("scHPF_", i))
  p$data$facet = names(factors)[ind]
  p <- p + 
    facet_wrap(~facet, strip.position = "left")+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    theme_umap()+
    theme(strip.placement = "inside",
          strip.text = element_text(size = 14),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0, 0, 0, 0, unit = "in"),
          legend.key.height = unit(0.25, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(t=0, r=0, b=0, l=0.05, 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t=0, r=0, b=0, l=-0.2, 'in'))
  
  plots <- c(plots, list(p))
  
  p <- Plot_Density_Custom(treated, features = paste0("scHPF_", i)) + 
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    theme_umap()+
    theme(plot.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "in"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          legend.key.height = unit(0.25, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(t=0, r=0, b=0, l=0.05, 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t=0, r=0, b=0, l=-0.2, 'in'))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(nrow = 2)

pdf(file.path(res_path, "analysis", paste0(fig, "panel_I.pdf")),  height = 4.5, width = 5.5)
print(p)
dev.off()

# PANEL J
#-------------------------------------------------------------------------------
plots <- list()
for(i in paste0("scHPF_", c(5,20))){
  
  thresh <- median(data@meta.data[,i]) + 2*mad(data@meta.data[,i])
  
  tmp <- data@meta.data %>% 
    rename_at(i, ~"factor") %>%
    mutate(class = ifelse(factor > thresh, "High", "Low"),
           class = factor(class, levels = c("Low", "High")),
           treatment = factor(treatment, levels = c("Ctrl", "Apop")))
  
  plot_data <- tmp %>%
    group_by(treatment, class) %>% 
    summarise(n = n()) %>%
    complete(class, fill = list(n = 0)) %>%
    group_by(treatment) %>% 
    reframe(class, n, prop = n/sum(n)) %>% 
    mutate(lab = paste0(round(prop*100, digits = 0), "%"),
           facet = names(factors)[factors %in% i])
  
  anno <- chisq.test(tmp$treatment, tmp$class)$p.value
  anno <- ifelse(anno < 0.001, "***", ifelse(anno < 0.01, "**", ifelse(anno < 0.05, "*", "ns")))
  
  p <- ggplot(plot_data, aes(x = treatment, y = prop, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label = lab, y = prop + 0.01),
              position = position_dodge(width = 0.9),
              color = 'black',
              vjust = 0,
              size = 3,
              family = "Helvetica") +
    geom_signif(comparisons = list(c("Ctrl", "Apop")),
                step_increase = 0.1, 
                tip_length = 0, 
                margin_top = 0.15,
                vjust = 0.5,
                annotation = anno,
                family = "Helvetica", 
                textsize = 4)+
    facet_wrap(~facet, strip.position = "top") +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1.1), 
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) formatC(x, format = "g"))+
    scale_fill_manual(name = "Expression status",
                      values = unname(proj_cols("light blue", "purple"))) +
    labs(y = "Proportion of cells")+
    theme_proj()+
    theme(panel.grid.major = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
          legend.title = element_text(size = 12, hjust = 1),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          plot.margin = margin(t=0.1, r=0, b=0.05, l=0.05, "in"),
          legend.margin = margin(l = -0.1, unit = "in"),
          strip.text = element_text(size = 14, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")),
          strip.background = element_blank(),
          strip.clip = "off",
          plot.background = element_rect(fill = "transparent"))
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + 
  plot_layout(nrow = 1, guides = "collect", axis_titles = "collect") &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal")

pdf(file.path(res_path, "analysis", paste0(fig, "panel_J.pdf")), height = 3.5, width = 4)
print(p)
dev.off()

# PANEL K
#-------------------------------------------------------------------------------
data <- NormalizeData(data, normalization.method = "LogNormalize", assay = "RNA") %>%
  ScaleData(assay = "RNA")

prop <- data.frame()
for(i in factors){
  
  tmp <- data@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

diff_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

genes <- data.frame()
for(i in diff_factors){
  genes <- rbind(genes, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}

keep_factors <- paste0("scHPF_", c(5, 15, 20))

keep_highlights <- data.frame()
for(i in keep_factors){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

data@meta.data <- data@meta.data %>% mutate(cols = paste0(line, "_", treatment))

id_cols <- c(proj_cols("teal", "purple", "orange1"), 
             lighten(proj_cols("teal", "purple", "orange1"), space ="HLS", amount = 0.5))
id_cols <- id_cols[c(1, 4, 2, 5, 3, 6)]
names(id_cols) <- c("iCW50036_Ctrl", "iCW50036_Apop", 
                    "iCW50118_Ctrl" ,"iCW50118_Apop", 
                    "iCW70347_Ctrl", "iCW70347_Apop")
names(id_cols) <- gsub("apop_neu", "AN", names(id_cols))

averages <- AverageExpression(data,
                              features = genes$name,
                              group.by = "cols",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data

averages <- averages[genes$name, names(id_cols)]
averages <- t(averages)

ha <- rowAnnotation(mod = anno_simple(names(id_cols), 
                                      col = id_cols, 
                                      height = unit(0.15, "in"), 
                                      gp = gpar(col = "white", lwd = 1.5)),
                    annotation_name_gp = gpar(fontsize = 0),
                    show_legend = F,
                    border = F)

col_split <- genes %>% filter(name %in% colnames(averages)) %>% pull(factor_name)
col_split <- factor(col_split, levels = unique(genes$factor_name))

Idents(data) <- "treatment"
sig <- FindMarkers(data, 
                   ident.1 = "Apop", 
                   ident.2 ="Ctrl",
                   features = keep_highlights$name,
                   assay = "RNA",
                   test.use = "MAST", 
                   latent.vars = c("nCount_RNA", "line"),
                   logfc.threshold = 0,
                   min.pct = 0) %>%
  rownames_to_column("gene") %>%
  mutate(gene = factor(gene, levels = keep_highlights$name),
         cols = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, proj_cols("red"), 
                       ifelse(avg_log2FC < 0 & p_val_adj < 0.05, proj_cols("blue"), "black"))) %>%
  arrange(gene)

bottom_ha <- HeatmapAnnotation(labs = anno_mark(at = which(colnames(averages) %in% keep_highlights$name), 
                                                labels = keep_highlights$name,
                                                side ="right",
                                                labels_gp = gpar(fontface = "italic", 
                                                                 fontfamily = "Helvetica", 
                                                                 fontsize = 10,
                                                                 col = sig$cols)),
                               which = "column")

cairo_pdf(file.path(res_path, "analysis", paste0(fig, "panel_K.pdf")), height = 5, width = 9)
set.seed(1234)
p <- Heatmap(averages, name = "mat",
             col = colorRamp2(breaks = c(min(averages), 0, max(averages)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             left_annotation = ha,
             top_annotation = columnAnnotation(foo = anno_empty(border = FALSE,  height = unit(0.1, "in"))),
             cluster_rows = T,
             row_names_side = "left",
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_dend_side = "right",
             row_dend_width = unit(0.2, "in"),
             column_split = col_split,
             row_title_side = "right",
             row_title_rot = 270,
             row_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 0),
             column_title_rot = 90,
             column_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 8, hjust = 0),
             column_names_side = "bottom",
             column_names_gp = gpar(fontface = "italic", fontfamily = "Helvetica", fontsize = 0),
             column_gap = unit(0.02, "in"),
             bottom_annotation = bottom_ha,
             cluster_columns = T,
             show_column_dend = FALSE,
             heatmap_legend_param = list(
               title = expression('Average scaled expression level'), 
               legend_width = unit(1, "in"),
               direction = "horizontal",
               title_position = "topcenter",
               border = "black"
             ))

draw(p, 
     align_heatmap_legend = "heatmap_center", 
     heatmap_legend_side = "top", 
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"),
     legend_title_gp = gpar(fontsize = 12, hjust = 1))

for(i in 1:length(diff_factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, height = unit(0.15, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}

dev.off()

#===============================================================================
#
# PANELS L, M: COMPOUND-TREATED HMC3s
#
#===============================================================================
dataset <- "HMC3stim"
res_path <- file.path(home_path, paste0("proj_", dataset))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names

# PANEL L
#-------------------------------------------------------------------------------
hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(0, 0, 0, -0.1, 'in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, "analysis", paste0(fig, "panel_L1.pdf")), height = 3, width = 3.5)
print(p)
dev.off()

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1]) 
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1]) 
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2]) 
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

plots <- list()
for(i in c("ctrl", "Camptothecin","Narciclasine", "Torin2")){
  
  tmp <- subset(data, subset = identity %in% ifelse(i == "ctrl", c("DMSO", "Untreated"), i))
  hex <- make_hexbin(tmp, nbins = 50, dimension_reduction = "schpf.umap")
  
  p <- ggplot() + 
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    scale_x_continuous(limits = c(xmin - abs(xmin)*0.05, xmax + abs(xmax)*0.05))+
    scale_y_continuous(limits = c(ymin - abs(ymin)*0.05, ymax + abs(ymax)*0.05))+
    labs(title = ifelse(i == "ctrl", "DMSO\nUntreated", i))+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(l=0.05, unit = 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, 'in'),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.title = element_text(margin = margin(b = 0.05, unit = 'in'), size = 14, color = "black", face = "bold", vjust = 0),
          plot.margin = margin(0, 0, 0, 0))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(nrow = 1)

pdf(file.path(res_path,"analysis",  paste0(fig, "panel_L2.pdf")), height = 3, width = 7)
print(p)
dev.off()

# PANEL M
#-------------------------------------------------------------------------------
plots <- list()
for(i in paste0("scHPF_", c(3,20))){
  
  thresh <- median(data@meta.data[,i]) + 2*mad(data@meta.data[,i])
  
  tmp <- data@meta.data %>% 
    rename_at(i, ~"factor") %>%
    mutate(class = ifelse(factor > thresh, "High", "Low"),
           class = factor(class, levels = c("Low", "High")),
           identity = recode(identity, "Untreated"="Untreated/DMSO", "DMSO"="Untreated/DMSO"),
           identity = factor(identity, levels = c("Untreated/DMSO", "Camptothecin", "Narciclasine", "Torin2")))
  
  plot_data <- tmp %>%
    group_by(identity, class) %>% 
    summarise(n = n()) %>%
    complete(class, fill = list(n = 0)) %>%
    group_by(identity) %>% 
    reframe(class, n, prop = n/sum(n)) %>% 
    mutate(lab = paste0(round(prop*100, digits = 0), "%"),
           facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i]))
  
  campto <- chisq.test(tmp$identity[tmp$identity %in% c("Untreated/DMSO", "Camptothecin")], 
                       tmp$class[tmp$identity %in% c("Untreated/DMSO", "Camptothecin")])$p.value
  narci <- chisq.test(tmp$identity[tmp$identity %in% c("Untreated/DMSO", "Narciclasine")], 
                      tmp$class[tmp$identity %in% c("Untreated/DMSO", "Narciclasine")])$p.value
  tor <- chisq.test(tmp$identity[tmp$identity %in% c("Untreated/DMSO", "Torin2")], 
                    tmp$class[tmp$identity %in% c("Untreated/DMSO", "Torin2")])$p.value
  
  anno <- c(campto, narci, tor)
  anno <- ifelse(anno < 0.001, "***", ifelse(anno < 0.01, "**", ifelse(anno < 0.05, "*", "ns")))
  
  p <- ggplot(plot_data, aes(x = identity, y = prop, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label = lab, y = prop + 0.01),
              position = position_dodge(width = 0.9),
              color = 'black',
              vjust = 0,
              size = 3,
              family = "Helvetica") +
    geom_signif(comparisons = list(c("Untreated/DMSO", "Camptothecin"),
                                   c("Untreated/DMSO", "Narciclasine"),
                                   c("Untreated/DMSO", "Torin2")),
                step_increase = 0.05,
                tip_length = 0,
                margin_top = 0.1,
                anno = anno,
                vjust = 0.5,
                family = "Helvetica",
                textsize = 3)+
    facet_wrap(~facet, strip.position = "top") +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1.25), 
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) formatC(x, format = "g"))+
    scale_fill_manual(name = "Expression\nstatus",
                      values = unname(proj_cols("light blue", "purple"))) +
    labs(y = "Proportion of cells")+
    theme_proj()+
    theme(panel.grid.major = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
          legend.title = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          plot.margin = margin(t=0.1, r=0, b=0.05, l=0.05, "in"),
          legend.margin = margin(l = -0.1, unit = "in"),
          strip.text = element_markdown(size = 14, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")),
          strip.background = element_blank())
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")

pdf(file.path(res_path, "analysis", paste0(fig, "panel_M.pdf")), height = 4, width = 9)
print(p)
dev.off()
