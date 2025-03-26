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
library(ggridges)
library(ggsignif)
library(ggtext)
library(grid)
library(lmerTest)
library(logger)
library(patchwork)
library(readxl)
library(schex)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
dataset <- "HMC3stim"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp12_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_helpers.R"))

# REFERENCE DATA
#-------------------------------------------------------------------------------
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
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
  labs(title = "Compound-treated<br>HMC3s")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(l=-0.2, unit='in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 3.5, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL B 
#
#===============================================================================
cols <- c(proj_cols_grey("light grey", "med grey"),  "#DF482E", "#4445A0", "#6095CB")
names(cols) <- c("Untreated", "DMSO", "Camptothecin", "Narciclasine", "Torin2")

fill_cols <- c("#F5F5F5", "#E0E0E0", "#D1E0F0", "#F7D1Ca", "#D4D5ED")
names(fill_cols) <- c("Untreated", "DMSO", "Torin2", "Camptothecin", "Narciclasine")

text_cols <- c("#282828", "#282828","#6095CB", "#DF482E", "#4445A0")
names(text_cols) <- c("Untreated", "DMSO", "Torin2", "Camptothecin", "Narciclasine")

keep_drugs <- sort(data@meta.data %>% group_by(identity) %>% summarise(n=n()) %>% filter(n > 500) %>% pull(identity))
keep_drugs <- factor(keep_drugs, levels = c("Untreated","DMSO", "Camptothecin", "Narciclasine", "Torin2"))
keep_drugs <- sort(keep_drugs)

ref_hex <- make_hexbin(data, nbins = 80, dimension_reduction = "schpf.umap")

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1]) 
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1]) 
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2]) 
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

plots <- list()
for(i in keep_drugs){
  
  tmp <- subset(data, subset = identity %in% i)
  hex <- make_hexbin(tmp, nbins = 50, dimension_reduction = "schpf.umap")
  plot_data <- data.frame(hex@misc$hexbin$hexbin.matrix, facet = i)
  
  p1 <- ggplot() + 
    geom_hex(data = plot_data, aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    scale_x_continuous(limits = c(xmin - abs(xmin)*0.05, xmax + abs(xmax)*0.05))+
    scale_y_continuous(limits = c(ymin - abs(ymin)*0.05, ymax + abs(ymax)*0.05))+
    facet_wrap(~facet)+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = 'in'), family = "Helvetica"),
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, -0.2, 'in'),
          plot.title = element_text(margin = margin(0, 0, 0.01, 0, 'in'), family = "Helvetica"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          strip.text = element_text(size = 16, color = text_cols[i], family = "Helvetica"),
          strip.placement = "outside",
          strip.background = element_rect(fill = fill_cols[i], color = cols[i], linewidth = 1.5))
  
  plots <- c(plots, list(p1))
}

p <- wrap_plots(plots) + plot_layout(nrow = 1)

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 3, width = 12)
print(p)
dev.off()

#===============================================================================
#
# PANEL C - RIDGE PLOT
#
#===============================================================================
cols <- c(proj_cols_grey("light grey", "med grey"),  "#DF482E", "#4445A0", "#6095CB")
names(cols) <- c("Untreated", "DMSO", "Camptothecin", "Narciclasine", "Torin2")

plots <- list()
for(i in factors){
  
  tmp <- data@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i]), 
           identity = factor(identity, levels = rev(names(cols))))
  
  p <- ggplot(tmp, aes(x = factor, y = identity, fill = identity)) + 
    geom_density_ridges()+
    geom_vline(xintercept = 0.01, lty = 2) + 
    scale_x_continuous(trans = "log2", labels = function(x) formatC(x, digits = 2,  format = "g")) + 
    scale_fill_manual(values = cols)+
    facet_wrap(~facet)+
    labs(x = "log(Cell score)")+
    theme_proj() +
    theme(strip.text = element_text(size = 10, margin = margin(b = 0.05, t = 0.05, unit = "in"), family = "Helvetica"),
          axis.text = element_text(size = 10),
          axis.title.x =  element_text(size = 12),
          axis.title.y =  element_blank(),
          legend.position = "none")
  
  if(grepl("-", names(factors)[factors %in% i])){
    p <- p + theme(strip.text = element_markdown(size = 10, family = "Helvetica", margin = margin(b = 0.05, t = 0.05, unit = "in")))
  } else {
    p <- p + theme(strip.text = element_text(size = 10, family = "Helvetica", margin = margin(b = 0.05, t = 0.05, unit = "in")))
  }
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + 
  plot_layout(guides = "collect", axes = "collect_y", axis_titles = "collect", nrow = 4)

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 7, width = 14)
print(p)
dev.off()

#===============================================================================
#
# PANEL D
#
#===============================================================================
prop <- data.frame()
for(i in factors){
  
  tmp <- data@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

keep_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

p1 <- prop %>%
  filter(i %in% keep_factors) %>%
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
        axis.title.x = element_text(size = 11, margin = margin(-0.3, 0, 0, 0, unit = "in"), family = "Helvetica"),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(size = 11, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5))

res <- data.frame()
for(i in keep_factors){
  for(j in c("Camptothecin", "Narciclasine", "Torin2")){
    
    tmp <- data@meta.data %>% 
      filter(identity %in%  c("Untreated", "DMSO", j)) %>%
      mutate(treatment = ifelse(identity %in% c("Untreated", "DMSO"), 0, 1),
             scale_UMI = c(scale(nCount_RNA))) %>%
      rename_at(all_of(i), ~"factor") %>%
      dplyr::select(treatment, factor, scale_UMI)
    
    out <- lm(log2(factor) ~ treatment + scale_UMI, data = tmp)
    
    coef <- coef(summary(out))
    confint <- confint(out)
    
    log2FC <- logfc(tmp$factor[tmp$treatment==1], tmp$factor[tmp$treatment==0])
    
    res <- rbind(res, cbind(factor = i,
                            drug = j,
                            log2FC,
                            n = nrow(subset), 
                            coef = coef["treatment", "Estimate"], 
                            se = coef["treatment", "Std. Error"], 
                            lowerCI = confint["treatment", "2.5 %"],
                            upperCI = confint["treatment", "97.5 %"],
                            p = coef["treatment", "Pr(>|t|)"]))
    
    rm(list = ls()[ls() %in% c("out", "coef", "confint", "log2FC")])
  }
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI", "p")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "bonferroni"),
         col = ifelse(log2FC > 0 & p.adj < 0.05, "Upregulated", "NS (Adj. *p* > 0.05)"),
         col = ifelse(log2FC < 0 & p.adj < 0.05, "Downregulated", col),
         col = factor(col, levels = c("Downregulated", "Upregulated", "NS (Adj. *p* > 0.05)")))

write_xlsx(res %>% dplyr::select(factor:p.adj), file.path(res_path, "HMC3stim.xlsx"))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated", "Upregulated", "NS (Adj. *p* > 0.05)")

plots <- list()
for(i in c("Camptothecin", "Narciclasine", "Torin2")){
  
  p <- res %>% 
    filter(drug == i) %>% 
    mutate(facet = i,
           facet = gsub('-high', paste0("<sup>high</sup>"), facet)) %>%
    ggplot(aes(x = log2FC, y = factor, fill = col)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5)+
    geom_col()+
    geom_vline(xintercept = -0.5, linewidth = 0.35, lty = 2)+
    geom_vline(xintercept = 0.5, linewidth = 0.35, lty = 2)+
    scale_fill_manual(values = cols)+
    labs(x = "log<sub>2</sub>FC(Treated vs. Untreated/DMSO)")+
    guides(size = guide_legend(title = "Log<sub>2</sub>FC", nrow = 1, override.aes = aes(shape = 19)),
           fill = guide_legend(title = "Significance", nrow = 1, override.aes = aes(shape = 21, size = 4)))+
    scale_x_continuous(labels = function(x) formatC(x, format = "g"),
                       limits = c(-4, 4), breaks = seq(-4, 4, 2))+
    facet_wrap(~facet)+
    theme_proj() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_markdown(size = 11, family = "Helvetica"),
          axis.text.x = element_text(size = 11, family = "Helvetica"),
          plot.margin = margin(0, r = 0.1, b=0.2, 0, "in"),
          panel.grid.major.x = element_line(linewidth = 0.5, color = "#F0F0F0"),
          panel.grid.major.y = element_blank(),
          axis.ticks.x = element_line(linewidth = 0.5), 
          axis.ticks.y = element_line(linewidth = 0.5), 
          panel.border = element_rect(linewidth = 1, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          legend.title = element_markdown(size = 10, family = "Helvetica", margin = margin(r = 0.05, unit = "in")),
          legend.text = element_markdown(size = 10, family = "Helvetica", margin = margin(l = 0.05, unit = "in")),
          legend.title.position = "left",
          legend.position = "top",
          legend.direction = "horizontal", 
          legend.key.spacing.y = unit(0.02, "in"),
          legend.spacing.x = unit(0, "in"),
          axis.text.y = element_blank(), 
          strip.text = element_text(size = 12, margin = margin(b = 0.05, t = 0.05, unit = "in"), color = text_cols[i]),
          strip.background = element_rect(fill = fill_cols[i]))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(nrow = 1, axis_titles = "collect", guides = "collect") & 
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.margin = margin(b = -0.05, unit = "in"))

p <- wrap_plots(p1, p) + plot_layout(nrow = 1, widths = c(0.25, 0.75)) + 
  theme(plot.margin = margin(l = 0.05, unit = "in"))

cairo_pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 5.5, width = 9)
print(p)
dev.off()

#===============================================================================
#
# PANEL E - GENE EXPRESSION HEATMAP
#
#===============================================================================
keep_factors <- res %>%
  filter(p.adj < 0.05 & abs(log2FC) > 0.5) %>%
  arrange(col, desc(abs(coef))) %>%
  mutate(factor = recode(factor, !!!factors)) %>%
  group_by(drug, col) %>%
  slice_head(n = 1) %>%
  pull(factor)

id_cols <- c(proj_cols_grey("light grey", "med grey"), "#DF482E", "#4445A0", "#6095CB")
names(id_cols) <- c("Untreated", "DMSO","Camptothecin", "Narciclasine", "Torin2")

data <- NormalizeData(data, normalization.method = "LogNormalize", assay = "RNA") %>%
  ScaleData(assay = "RNA")

diff_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

genes <- data.frame()
for(i in diff_factors){
  genes <- rbind(genes, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}
exp_genes <- names(which(rowSums(data@assays$RNA@data)[genes$name] > 0))
genes <- genes %>% filter(name %in% exp_genes)

keep_highlights <- data.frame()
for(i in keep_factors){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

averages <- AverageExpression(data,
                              features = unique(genes$name),
                              group.by = "identity",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data
averages <- averages[genes$name, names(id_cols)]

keep_highlights <- keep_highlights %>% filter(name %in% rownames(averages))

ha <- HeatmapAnnotation(mod = anno_simple(names(id_cols), 
                                          col = id_cols, 
                                          height = unit(0.15, "in"), 
                                          gp = gpar(col = "white", lwd = 1.5)),
                        annotation_name_gp = gpar(fontsize = 0),
                        show_legend = F,
                        border = F, 
                        which = "column")

row_split <- genes %>% filter(name %in% rownames(averages)) %>% pull(factor_name)
row_split <- factor(row_split, levels = unique(genes$factor_name))

right_ha <- rowAnnotation(labs = anno_mark(at = match(keep_highlights$name, rownames(averages)), 
                                           labels = keep_highlights$name,
                                           side ="right",
                                           labels_gp = gpar(fontface = "italic", fontfamily = "Helvetica", fontsize = 8)))

cairo_pdf(file.path(res_path, paste0(fig, "panel_E.pdf")), height = 8, width = 4)
set.seed(1234)
p <- Heatmap(averages, name = "mat",
             col = colorRamp2(breaks = c(min(averages), 0, max(averages)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             top_annotation = ha,
             left_annotation = rowAnnotation(foo = anno_empty(border = FALSE, width = unit(0.1, "in"))),
             cluster_rows = T,
             show_row_names = FALSE,
             row_split = row_split,
             row_title_side = "left",
             row_title_rot = 0,
             row_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 8),
             column_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 8),
             column_names_rot = 45,
             column_names_side = "bottom",
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_gap = unit(0.02, "in"),
             right_annotation = right_ha,
             column_dend_height = unit(0.3, "in"),
             cluster_columns = T,
             show_row_dend = FALSE,
             heatmap_legend_param = list(
               title = expression('Avg. scaled exp. level'), 
               legend_width = unit(1, "in"),
               direction = "horizontal",
               title_position = "topcenter",
               border = "black"
             ))

draw(p, 
     align_heatmap_legend = "heatmap_center", 
     heatmap_legend_side = "top", 
     background = "transparent", 
     padding = unit(c(b=0.1, l=0.2, b=0.3, r=0.1), "in"),
     legend_title_gp = gpar(fontsize = 12, hjust = 0.5, fontfamily = "Helvetica"))

for(i in 1:length(diff_factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}
dev.off()

