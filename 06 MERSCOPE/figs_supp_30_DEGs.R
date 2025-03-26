#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE
# Figure S30.	DEGs across cell types in IR-High and IR-Low neighborhoods.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/path/to/rlib")

library(patchwork)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggtext)
library(writexl)
library(ggrepel)
library(ggridges)
library(readxl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/2025-02-14_cellular_niches")
figs <- "fig_s30"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/load_merscope.R"))
adjacent_cells <- readRDS(file.path(res_path, "adjacent_cells.rds"))

layer_labs <- c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB")
names(layer_labs) <- c("WM", "GM (L2-6)", "GM (L1)")

#===============================================================================
# 
# DEGS
#
#===============================================================================
for(factor_name in c(20, 26)){
  all_res <- data.frame()
  for(niche_name in c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    for(j in c("Astrocytes", "Endothelial", "Excitatory", "Inhibitory", "Oligodendrocytes", "OPCs", "Pericytes")){
      
      neighbors <- adjacent_cells %>% 
        filter(cell_type == j & niche == niche_name & factor == factor_name) %>% 
        pull(neighbors)
      neighbors <- unlist(str_split(neighbors, "[,]"))
      neighbors <- names(table(neighbors))[which(table(neighbors) > 1)]
      
      keep_cells <- rownames(data@meta.data %>% filter(layer_niches == niche_name & cell_type == j))
      tmp <- subset(data, subset = cell_names %in% keep_cells)
      tmp <- NormalizeData(tmp, assay = "RNA") %>% 
        ScaleData(assay = "RNA", features = rownames(tmp))
      
      tmp@meta.data <- tmp@meta.data %>% 
        mutate(lab = ifelse(cell_names %in% neighbors, 1, 0))
      Idents(tmp) <- "lab"
      
      if(sum(table(tmp@meta.data$lab) > 100) == 2){
        resDE <- FindMarkers(tmp, 
                             ident.1 = "1",
                             ident.2 = "0",
                             test.use = "MAST", 
                             assay = "RNA",
                             slot = "data",
                             latent.vars = "nCount_RNA",
                             max.cells.per.ident = 2500,
                             random.seed = 1234, 
                             logfc.threshold = 0.25, 
                             verbose = F) %>% 
          rownames_to_column("genes")
        
        all_res <- rbind(all_res, cbind(cell_type = j, niche = niche_name, factor = factor_name, resDE))
      }
    }
  }
  
  write_xlsx(all_res, file.path(res_path, paste0("tab_adjacent_DEGs", factor_name, ".xlsx")))
}

#===============================================================================
# 
# PLOT
#
#===============================================================================
cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("1", "2", "3")
col_labs <- c("Downregulated", "Upregulated", "NS (Adj. p > 0.05)")

layer_labs <- c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")
names(layer_labs) <- c("WM", "GM (L1)", "GM (L2-6)")

plot_data <- read_xlsx(file.path(res_path, paste0("tab_adjacent_DEGs20.xlsx"))) %>% 
  mutate(log10.p.adj = -log10(p_val_adj),
         col = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, "2", "3"),
         col = ifelse(avg_log2FC < 0 & p_val_adj < 0.05, "1", col),
         lab = ifelse(p_val_adj < 0.05, genes, NA),
         niche = recode(niche, !!!layer_labs)) %>%
  filter(niche == "GM.SLC17A7.GAS7.HSP90AB1")

p <- ggplot(plot_data, aes(x = avg_log2FC, y = log10.p.adj, color = col)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_point(size = 1, alpha = 0.6)+
  geom_vline(xintercept = 0,col = "black")+
  scale_color_manual(values = cols,
                     labels = col_labs)+
  guides(color = guide_legend(override.aes = aes(size = 4)))+
  scale_x_continuous(labels = function(x) formatC(x, format = "g"))+
  scale_y_continuous(labels = function(x) formatC(x, format = "g"))+
  geom_text_repel(data = plot_data %>% drop_na(), aes(label = lab), 
                  family = "Helvetica", seed = 1234, size = 3, show.legend = F) +
  labs(x = "Log FC (IR<sup>High</sup> vs. IR<sup>Low</sup>)", 
       y = "log<sub>10</sub>(Adj. <i>p</i>-value)")+
  facet_grid(~ cell_type) +
  theme_proj()+
  theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
        legend.spacing.x = grid::unit(0.05, "in"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(b = -0.2, unit = "in"),
        panel.border = element_rect(linewidth = 1, color = "black"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        strip.text = element_text(size = 10, margin = margin(t = 0.05, b = 0.05, unit = "in")))

pdf(file.path(res_path, paste0(figs, "scHPF_20_DEGs.pdf")), height = 3, width = 10)
print(p)
dev.off()
