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
library(ComplexHeatmap)
library(ggrastr)
library(ggtext)
library(Hmisc)
library(patchwork)
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
fig <- "fig_supp01_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))

#===============================================================================
#
# PANEL A - SEX
#
#===============================================================================
sex_labs <- c("M", "F")
names(sex_labs) <- c("Males", "Females")

plots <- list()
for(i in c("M", "F")){
  
  tmp <- subset(data, subset = SEX %in% i)
  hex <- make_hexbin(tmp, nbins = 80, dimension_reduction = "schpf.umap")
  
  p <- ggplot() +
    geom_point(data = ref_hex, aes(x = x, y = y), color = "#F0F0F0") +
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    labs(title = paste0("**", names(sex_labs)[sex_labs %in% i], "** (N=", format(length(Cells(tmp)), big.mark= ","), ")"))+
    theme_umap()+
    theme(plot.title = element_markdown(size = 14, family = "Helvetica", margin = margin(t = 0, r = 0, b = -0.05, l = 0, unit = "in")),
          panel.border = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.title = element_text(size = 12, margin = margin(b = 0.05, unit = "in")),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(ncol = 2)

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 3, width = 7)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
age <- data@meta.data %>% group_by(`Donor.ID`) %>% sample_n(1) %>% pull(AGE)

annot <- paste0("M (SD) = ", format(mean(age), digits = 3), " (", format(sd(age), digits = 3), ")")

p <- FeaturePlot(data, feature = "AGE", reduction = 'schpf.umap') 

p <- ggplot(p$data, aes(x = schpfumap_1, y = schpfumap_2, color = AGE))+
  geom_point_rast(size = 0.1, scale = 0.3)+
  scale_color_gradientn(colours = magma(10, direction = 1))+
  guides(color = guide_colourbar(ticks.linewidth = 0.5))+
  labs(title = "**Age** (years)", subtitle = annot)+
  theme_umap() +
  theme(plot.title = element_markdown(size = 14, family = "Helvetica", margin = margin(b = 0.05, unit = 'in')),
        plot.subtitle = element_text(size = 10, hjust = 0.5, family = "Helvetica", margin = margin(b = 0.05, unit = 'in')),
        legend.text = element_text(size = 10, family = "Helvetica", margin = margin(l = 0.05, unit = 'in')),
        legend.title = element_blank(),
        legend.key.height = unit(0.2, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.margin = margin(l = -0.2, unit = "in"),
        legend.ticks.length = unit(0.025, "in"))

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 3, width = 3)
print(p)
dev.off()
#===============================================================================
#
# PANEL C - DX
#
#===============================================================================
dx_names <- c("AD", "ALS/FTD", "Control", "DNET", "ET", "HD", "MCI", "MS", "PART", "PD", "Stroke", "TLE")
names(dx_names) <- sort(unique(data@meta.data$DX_CAT))

dx_levels <- data@meta.data %>% 
  group_by(DX_CAT) %>% 
  summarise(n = n()) %>% 
  filter(n > 10000) %>%
  pull(DX_CAT)

plots <- list()
for(i in dx_levels){
  
  tmp <- subset(data, subset = DX_CAT %in% i)
  hex <- make_hexbin(tmp, nbins = 80, dimension_reduction = "schpf.umap")
  
  p = ggplot() +
    geom_point(data = ref_hex, aes(x = x, y = y), color = "#F0F0F0") +
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    labs(title = paste0("**", i, "**<br>(N=", format(length(Cells(tmp)), big.mark= ","), ")"))+
    theme_umap()+
    theme(plot.title = element_markdown(size = 14, family = "Helvetica", margin = margin(t = 0, r = 0, b = -0.05, l = 0, unit = "in")),
          panel.border = element_blank(),
          axis.line = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.title = element_text(size = 12, margin = margin(b = 0.05, unit = "in")),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 1)

pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 2.5, width = 16)
print(p)
dev.off()

rm(list = c("ref_hex", "hex", "tmp")); gc()
#===============================================================================
#
# PANEL D
#
#===============================================================================
prop <- data.frame()
for(i in factors){
  
  tmp <- schpf@cell.embeddings %>% as.data.frame() %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

p <- prop %>%
  mutate(prop = as.numeric(as.character(prop)), 
         factor = factor(i, levels = factors, labels = names(factors)),
         lab = paste0(round(as.numeric(prop)*100, 1), "%")) %>%
  ggplot(aes(y = prop, x = factor)) +
  geom_col(fill = proj_cols("yellow")) +
  geom_text(aes(label = lab), hjust = 1, angle = 90, family = "Helvetica")+
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.25), 
                     labels = function(x) formatC(x, digits = 2, format = "g"))+
  labs(y = "Proportion of\nexpressing cells", x = "Expression programs")+
  theme_proj()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line.y = element_line(linewidth = 0.75, color = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12, hjust = 1, angle = 45))

cairo_pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 4, width = 9)
print(p)
dev.off()

#===============================================================================
#
# PANEL E
#
#===============================================================================

# GENES
#-------------------------------------------------------------------------------
mat <- as.data.frame(rcorr(as.matrix(schpf@feature.loadings), type = "spearman")$r)
mat <- as.matrix(mat[factors, factors])
colnames(mat) <- gsub("scHPF_", "", colnames(mat))
rownames(mat) <- gsub("scHPF_", "", rownames(mat))

cairo_pdf(file.path(res_path, paste0(fig, "panel_E1.pdf")), height = 4, width = 3.5)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(-1, 0, 1), 
                              colors = c(proj_cols("blue"), "white", proj_cols("red"))),
             cluster_rows = T,
             row_names_side = "left",
             row_dend_side = "right",
             row_title = "Expression programs",
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 10),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 9),
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 9),
             cluster_columns = T,
             column_names_centered = F,
             column_title = "Expression programs",
             column_title_side = "bottom",
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 10),
             heatmap_legend_param = list(
               title = expression('Spearman '*rho), 
               legend_width = unit(1, "in"),
               direction = "horizontal",
               title_position = "topcenter",
               border = "black"
             ))

draw(p, heatmap_legend_side = "top", background = "transparent", padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2)) })
dev.off()

# CELLS
#-------------------------------------------------------------------------------
mat <- as.data.frame(rcorr(as.matrix(schpf@cell.embeddings), type = "spearman")$r)
mat <- as.matrix(mat[factors, factors])
colnames(mat) <- gsub("scHPF_", "", colnames(mat))
rownames(mat) <- gsub("scHPF_", "", rownames(mat))

cairo_pdf(file.path(res_path, paste0(fig, "panel_E2.pdf")), height = 4, width = 3.5)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(-1, 0, 1), 
                              colors = c(proj_cols("blue"), "white", proj_cols("red"))),
             cluster_rows = T,
             row_names_side = "left",
             row_dend_side = "right",
             row_title = "Expression programs",
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 10),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 9),
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 9),
             cluster_columns = T,
             column_names_centered = F,
             column_title = "Expression programs",
             column_title_side = "bottom",
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 10),
             heatmap_legend_param = list(
               title = expression('Spearman '*rho), 
               legend_width = unit(1, "in"),
               direction = "horizontal",
               title_position = "topcenter",
               border = "black"
             ))

draw(p, heatmap_legend_side = "top", background = "transparent", padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2)) })
dev.off()
