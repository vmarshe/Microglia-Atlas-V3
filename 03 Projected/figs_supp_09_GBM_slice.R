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
library(ggtext)
library(ggsignif)
library(grid)
library(lmerTest)
library(logger)
library(patchwork)
library(schex)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tidyverse)
library(writexl)

source("~/projects/microglia_atlas/src/00 Helpers/theme.R")
source("~/projects/microglia_atlas/src/00 Helpers/factors.R")
source("~/projects/microglia_atlas/src/00 Helpers/funcs_helpers.R")

# PATHS
#-------------------------------------------------------------------------------
dataset <- "GBM_slice"
home_path <- "/mnt/vast/hpc/homes/vsm2116/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp9_"

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

data[["pct_mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^MTRNR")
data@meta.data$Drug_Tx <- gsub(" ", "", data@meta.data$Drug_Tx)

#===============================================================================
#
# PANEL A - SCHEX PLOTS
#
#===============================================================================

# HEX
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

pdf(file.path(res_path, paste0(fig, "panel_A1.pdf")), height = 3, width = 3.5)
print(p)
dev.off()

keep_drugs <- sort(data@meta.data %>% group_by(Drug_Tx) %>% summarise(n = n()) %>% filter(n >= 1000) %>% pull(Drug_Tx))
keep_drugs <- keep_drugs[keep_drugs!="Ctrl"]

drugs_for_analysis <- data@meta.data %>% 
  filter(Drug_Tx!="Ctrl") %>% 
  group_by(Drug_Tx, TB_ID) %>% 
  summarise(n=n()) %>% 
  filter(n > 200) %>% 
  group_by(Drug_Tx) %>% 
  summarise(n=n()) %>% 
  filter(n>2) %>% 
  pull(Drug_Tx)

plots <- list()
for(i in keep_drugs){
  
  keep_id <- data@meta.data %>% filter(Drug_Tx == i) %>% group_by(TB_ID, Drug_Tx) %>% pull(TB_ID)
  keep_id <- unique(keep_id)
  
  tmp <- subset(data, subset = TB_ID %in% keep_id)
  hex_ctrl <- make_hexbin(subset(tmp, subset = Drug_Tx == "Ctrl"), nbins = 50, dimension_reduction = "schpf.umap")
  hex_trt <- make_hexbin(subset(tmp, subset = Drug_Tx == i), nbins = 50, dimension_reduction = "schpf.umap")
  
  p <- ggplot() + 
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(hex_ctrl@misc$hexbin$hexbin.matrix), aes(x = x, y = y), fill = proj_cols_grey("med grey"), stat = "identity") +
    geom_hex(data = data.frame(hex_trt@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    ggtitle(i)+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(l=0.05, unit = 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, 'in'),
          plot.title = element_text(margin = margin(0, 0, 0.01, 0, 'in'), size = 14, 
                                    color = ifelse(i %in% drugs_for_analysis, proj_cols("blue"), "black"),
                                    face = ifelse(i %in% drugs_for_analysis, "bold", "plain")),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0.05, 0.05, 0, 0.05))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(ncol = 5)

pdf(file.path(res_path, paste0(fig, "panel_A2.pdf")), height = 3.5, width = 9)
print(p)
dev.off()

#===============================================================================
#
# PANEL B - PROPORTION BARPLOT
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

p <- prop %>%
  filter(i %in% keep_factors) %>%
  mutate(prop = as.numeric(as.character(prop)), 
         factor = factor(i, levels = rev(factors), labels = rev(names(factors))),
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
        axis.text.y = element_text(size = 11, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5))

cairo_pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 5, width = 4.5)
print(p)
dev.off()
#===============================================================================
#
# PANEL C - VOLCANO PLOT
#
#===============================================================================
res <- data.frame()
for(j in drugs_for_analysis){
  
  tmp <- data@meta.data %>% filter(Drug_Tx %in% c(j, "Ctrl"))
  keep_id <- tmp %>% 
    group_by(TB_ID, Drug_Tx) %>% 
    summarise(n = n()) %>% 
    group_by(TB_ID) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>% 
    pull(TB_ID)
  
  for(i in keep_factors){
    
    subset <- tmp %>% 
      filter(TB_ID %in% keep_id & Drug_Tx %in% c("Ctrl", j)) %>% 
      rename_at(i, ~"factor") %>%
      mutate(treatment = ifelse(Drug_Tx == "Ctrl", 0, 1), 
             scale_UMI = scale(nCount_RNA), 
             scale_MT = scale(pct_mt)) %>%
      dplyr::select(treatment, scale_UMI, scale_MT, TB_ID, factor)
    
    out <- lmer(log2(factor) ~ scale_UMI + scale_MT + treatment + (1 | TB_ID),
                data = subset, 
                na.action = na.omit)
    
    coef <- coef(summary(out))
    confint <- confint(out)
    
    log2FC <- logfc(subset$factor[subset$treatment==1], subset$factor[subset$treatment==0])
    
    rand_eff <- VarCorr(out)$TB_ID[1]
    
    res <- rbind(res, cbind(factor = i,
                            treatment = j,
                            log2FC,
                            n = nrow(subset), 
                            coef = coef["treatment", "Estimate"], 
                            se = coef["treatment", "Std. Error"], 
                            lowerCI = confint["treatment", "2.5 %"],
                            upperCI = confint["treatment", "97.5 %"],
                            p = coef["treatment", "Pr(>|t|)"],
                            rand_eff))
    
    rm(list = ls()[ls() %in% c("out", "coef", "subset", "log2FC", "rand_eff","confint")]) 
  }
  
  rm(list = c("keep_id", "tmp"))
}

res <- res %>%
  mutate_at(all_of(c("n","log2FC", "coef", "se", "lowerCI", "upperCI","p", "rand_eff")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "bonferroni"),
         col = ifelse(log2FC > 0.5 & p.adj < 0.05, "Upregulated", "NS (Adj. *p* > 0.05)"),
         col = ifelse(log2FC < -0.5 & p.adj < 0.05, "Downregulated", col),
         col = ifelse(abs(log2FC) < 0.5 & p.adj < 0.05, "FC < 0.5", col),
         col = factor(col, levels = c("Downregulated", "Upregulated", "FC < 0.5", "NS (Adj. *p* > 0.05)")))

res %>% 
  dplyr::select(factor:p.adj) %>%
  write_xlsx(file.path(res_path, paste0(fig, dataset, "_res.xlsx")))

cols <- c(unname(proj_cols("blue", "pink")), "black",unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated", "Upregulated", "FC < 0.5", "NS (Adj. *p* > 0.05)")

p <- res %>%
  mutate(lab = ifelse(col %in% c("Downregulated", "Upregulated", "FC < 0.5"), 
                      paste0(treatment, " (", gsub(" \\([0-9]+\\)", "", factor), ")"), "")) %>%
  ggplot(aes(x = log2FC, y = -log10(p.adj), color =  col)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_point(alpha = 0.6)+
  geom_vline(xintercept = 0, col = "black")+
  scale_color_manual(values = cols)+
  guides(color = guide_legend(nrow = 1, title = "Significance", override.aes = aes(shape = 19, size = 4)))+
  scale_x_continuous(labels = function(x) formatC(x, format = "g"))+
  scale_y_continuous(labels = function(x) formatC(x, format = "g"))+
  geom_text_repel(aes(label = lab), family = "Helvetica", seed = 1234, size = 4, show.legend = F) +
  labs(x = "Log<sub>2</sub>FC (Treated vs. Control)",
       y = "-log<sub>10</sub>(Adj. *p*-value)")+
  theme_proj()+
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"),
        panel.grid.major = element_blank(),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        axis.title.y = element_markdown(size = 14, family = "Helvetica"),
        axis.title.x = element_markdown(size = 14, family = "Helvetica"),
        legend.title = element_markdown(size = 12, family = "Helvetica", face = "bold"),
        legend.text = element_markdown(size = 12, family = "Helvetica", margin = margin(l = 0.05, unit = "in")),
        legend.margin = margin(b = -0.2, unit = "in"),
        legend.key.spacing.x = unit(0.05, "in"), 
        legend.position = "top", 
        legend.direction = "horizontal",
        legend.justification = "left")

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 6, width = 7)
print(p)
dev.off()
#===============================================================================
#
# PANEL D - BARPLOTS
#
#===============================================================================
keep_factors <- res %>% 
  filter(treatment == "Topotecan" & p.adj < 0.001 & abs(log2FC) > 0.5) %>% 
  arrange(col, desc(abs(log2FC))) %>% 
  group_by(col) %>%
  slice_head(n = 3) %>%
  mutate(factor = recode(factor, !!!factors)) %>%
  pull(factor)

plots <- list()
for(i in keep_factors){
  
  thresh <- median(data@meta.data[,i]) + 2*mad(data@meta.data[,i])
  
  tmp <- data@meta.data %>% filter(Drug_Tx %in% c("Topotecan", "Ctrl"))
  
  keep_id <- tmp %>% 
    group_by(TB_ID, Drug_Tx) %>% 
    summarise(n = n()) %>% 
    group_by(TB_ID) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>% 
    pull(TB_ID)
  
  tmp <- tmp %>%
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
           facet = gsub(" \\(",  "\n\\(", names(factors)[factors %in% i]))
  
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
                margin_top = 0.1,
                vjust = 0,
                annotation = anno,
                family = "Helvetica", textsize = 4)+
    facet_wrap(~facet, strip.position = "top") +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1.2), 
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
          strip.text = element_text(size = 12, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")),
          strip.background = element_blank())
  
  plots <- c(plots, list(p))
}

plots[[1]] <- plots[[1]] + theme(axis.text.x = element_blank())
plots[[2]] <- plots[[2]] + theme(axis.text.x = element_blank())
plots[[3]] <- plots[[3]] + theme(axis.text.x = element_blank())

p <- wrap_plots(plots) + 
  plot_layout(guides = "collect", axis_title = "collect", nrow = 2) &
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        plot.margin = margin(t=0.05, b=0.05, l=0.05, r=0.05,'in'))

pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 5.5, width = 6.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL E - HEATMAP
#
#===============================================================================
keep_id <- data@meta.data %>% 
  filter(Drug_Tx %in% c("Topotecan", "Ctrl")) %>% 
  group_by(TB_ID, Drug_Tx) %>% 
  summarise(n = n()) %>% 
  group_by(TB_ID) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) %>% 
  pull(TB_ID)

sig_facs <- res %>% 
  filter(treatment == "Topotecan" & p.adj < 0.001 & abs(log2FC) > 0.5) %>% 
  mutate(factor = recode(factor, !!!factors)) %>%
  pull(factor)
sig_facs <- unique(sig_facs)

tmp <- subset(data, subset = TB_ID %in% keep_id & Drug_Tx %in% c("Ctrl", "Topotecan"))
tmp@meta.data$new_id <- paste0(gsub(" ", "", tmp@meta.data$TB_ID), "_", tmp@meta.data$Drug_Tx)
tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", assay = "RNA")

keep_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

genes <- data.frame()
for(i in keep_factors){
  genes <- rbind(genes, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}
exp_genes <- names(which(rowSums(data@assays$RNA@data)[genes$name] > 0))
genes <- genes %>% filter(name %in% exp_genes)

keep_highlights <- data.frame()
for(i in sig_facs){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

averages <- AverageExpression(tmp,
                              features = genes$name,
                              group.by = "new_id",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data

averages <- averages[genes$name, ]
averages <- t(averages)

id_cols <- c(proj_cols("teal", "purple", "yellow"))
id_cols <- c(id_cols, lighten(id_cols, space ="HLS", amount = 0.5))
id_cols <- id_cols[c(4, 1, 5, 2, 6, 3)]
names(id_cols) <- rownames(averages)

ha <- rowAnnotation(mod = anno_simple(names(id_cols), 
                                      col = id_cols, 
                                      height = unit(0.15, "in"), 
                                      gp = gpar(col = "white", lwd = 1.5)),
                    annotation_name_gp = gpar(fontsize = 0),
                    show_legend = F,
                    border = F)

row_split <- factor(gsub("_Ctrl|_Topotecan", "", rownames(averages)),
                    levels = unique(gsub("_Ctrl|_Topotecan", "", rownames(averages))))

col_split <- genes %>% filter(name %in% colnames(averages)) %>% pull(factor_name)
col_split <- factor(col_split, levels = unique(genes$factor_name))

Idents(tmp) <- "Drug_Tx"
sig <- FindMarkers(tmp, 
                   ident.1 = "Topotecan", 
                   ident.2 ="Ctrl",
                   features = unique(keep_highlights$name),
                   assay = "RNA",
                   test.use = "MAST", 
                   latent.vars = c("nCount_RNA"),
                   logfc.threshold = 0,
                   min.pct = 0) %>%
  rownames_to_column("gene") %>%
  filter(gene %in% colnames(averages)) %>%
  mutate(gene = factor(gene, levels = unique(colnames(averages))),
         cols = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, proj_cols("red"), 
                       ifelse(avg_log2FC < 0 & p_val_adj < 0.05, proj_cols("blue"), "black")))%>%
  rowwise() %>%
  mutate(at = which(colnames(averages) %in% gene)) %>%
  ungroup() %>%
  arrange(at)

bottom_ha <- HeatmapAnnotation(labs = anno_mark(at = sig$at, 
                                                labels = sig$gene,
                                                side ="right",
                                                labels_gp = gpar(fontface = "italic", 
                                                                 fontfamily = "Helvetica", 
                                                                 fontsize = 10,
                                                                 col = sig$cols)),
                               which = "column")

cairo_pdf(file.path(res_path, paste0(fig, "panel_E.pdf")), height = 4.5, width = 14)
set.seed(1234)
p <- Heatmap(averages, name = "mat",
             col = colorRamp2(breaks = c(min(averages), 0, max(averages)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             left_annotation = ha,
             top_annotation = columnAnnotation(foo = anno_empty(border = FALSE, height = unit(0.1, "in"))),
             cluster_rows = T,
             row_names_side = "left",
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_title = "sgRAN guide",
             row_title_side = "left",
             row_title_rot = 90,
             row_dend_side = "right",
             row_dend_width = unit(0.2, "in"),
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 0),
             column_split = col_split,
             column_gap = unit(0.02, "in"),
             column_title_rot = 45,
             column_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 10),
             show_column_names = F,
             cluster_columns = T,
             column_title_side = "top",
             show_column_dend = FALSE,
             bottom_annotation = bottom_ha,
             heatmap_legend_param = list(
               title = expression('Average\nscaled\nexpression\nlevel'), 
               legend_height = unit(1, "in"),
               border = "black"
             ))

draw(p, 
     heatmap_legend_side = "right", 
     background = "transparent", 
     padding = unit(c(0.2, 0.2, .2, 0.2), "in"),
     legend_title_gp = gpar(fontsize = 12, hjust = 1), align_heatmap_legend = "heatmap_center")

for(i in 1:length(keep_factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, height = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}

dev.off()