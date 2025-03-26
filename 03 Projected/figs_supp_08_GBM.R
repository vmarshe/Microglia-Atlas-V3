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

# PATHS
#-------------------------------------------------------------------------------
dataset <- "GBM"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp8_"

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

data@meta.data <- data@meta.data %>% 
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names

# below meta data obtained from records under https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103224.
meta <- data.frame(orig.ident = c("PJ016", "PJ017", "PJ018", "PJ025", 'PJ030', "PJ032", "PJ035", "PJ048"),
                   SEX = c("F", "M", "M", "M", "F", "F", "M", "M"),
                   AGE = c(49, 62, 65, 74, 56, 63, 50, 59), 
                   TYPE = c("proneural", "mesenchymal", "proneural", "classical", "classical", "mesenchymal", "classical", "proneural"),
                   ORIGIN = c("R frontal", "L temporal", "L temporal","R frontal","L temporal","L temporal","L temporal","R parietal"))

meta <- meta %>% arrange(TYPE, orig.ident)

data@meta.data$orig.ident <- recode(data@meta.data$orig.ident, "PJ030_701"="PJ030")
data@meta.data <- data@meta.data %>% left_join(meta, by="orig.ident")
rownames(data@meta.data) <- data@meta.data$cell_names

data <- subset(data, subset = orig.ident %in% meta$orig.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-|^MTRN|^MTRF|^MTND|^MTATP|^MTCYB|^MTCO")

# COLORS and LABELS
#-------------------------------------------------------------------------------
id_labs <- factor(meta$orig.ident)
names(id_labs) <- meta$orig.ident

meta$col <- NULL
meta$col[meta$TYPE == "classical"] <- lighten(proj_cols("orange1"), c(0, 0.3, 0.7), space="HCL")
meta$col[meta$TYPE == "mesenchymal"] <- lighten(proj_cols("teal"), c(0, 0.3), space="HCL")
meta$col[meta$TYPE == "proneural"] <- lighten(proj_cols("purple"), c(0, 0.3, 0.6), space="HCL")

id_cols <- meta$col
names(id_cols) <- meta$orig.ident

type_cols <- proj_cols("orange1","teal", "purple")
names(type_cols) <- c("classical","mesenchymal", "proneural")
type_labs <- c("classical","mesenchymal", "proneural")
names(type_labs) <- c("Classical","Mesenchymal", "Proneural")

sex_cols <- proj_cols("blue","yellow")
names(sex_cols) <- c("F","M")
sex_labs <- c("F", "M")
names(sex_labs) <- c("Female", "Male")

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
  labs(title = "GBM-derived myeloid<br>cells (scRNA-seq)")+
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
p1 <- DimPlot(data, group.by = "orig.ident") + 
  ggtitle("Donor ID") + 
  scale_color_manual(values = id_cols, breaks = id_labs, labels = names(id_labs))+
  guides(color = guide_legend(override.aes = aes(size = 3), nrow=2))+
  theme_umap()+
  theme(legend.key = element_rect(color = NA),
        legend.margin = margin(-0.2,0,0,0, "in"),
        legend.key.height = unit(0.05, "in"),
        legend.key.width = unit(0, "in"),
        legend.key.spacing.x = unit(0.1, 'in'),
        legend.key.spacing.y = unit(0.1, 'in'),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.text = element_text(size = 10, margin = margin(0, 0, 0, 0.02, "in")),
        plot.title = element_text(size =14, margin = margin(0, 0, 0, -0.2, "in"), family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"))

p2 <- DimPlot(data, group.by = "SEX") + 
  ggtitle("Sex") + 
  scale_color_manual(values = sex_cols, labels = names(sex_labs))+
  guides(color = guide_legend(override.aes = aes(size = 4), nrow=1))+
  theme_umap()+
  theme(legend.key = element_rect(color = NA),
        legend.margin = margin(-0.2,0,0,0, "in"),
        legend.key.height = unit(0.2, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(0.05, 'in'),
        plot.title = element_text(size =14, margin = margin(0, 0, 0, -0.2, "in"), family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"))

p3 <- DimPlot(data, group.by = "TYPE") + 
  ggtitle("Tumor Type") + 
  scale_color_manual(values = type_cols, breaks = type_labs, labels = names(type_labs))+
  guides(color = guide_legend(override.aes = aes(size = 4), nrow=2))+
  theme_umap()+
  theme(legend.key = element_rect(color = NA),
        legend.margin = margin(-0.2,0,0,0, "in"),
        legend.key.height = unit(0.05, "in"),
        legend.key.width = unit(0.15, "in"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.text = element_text(size = 10, margin = margin(0, 0, 0, 0.02, "in")),
        legend.key.spacing.x = unit(0.1, 'in'),
        legend.key.spacing.y = unit(0.1, 'in'),
        plot.title = element_text(size =14, margin = margin(0, 0, 0, -0.2, "in"), family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linewidth = 1.5, color = "black"),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"))

p <- wrap_plots(p1, p2, p3) + plot_layout(nrow = 1)

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 3.5, width = 8)
print(p)
dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
p <- data@meta.data %>%
  group_by(orig.ident) %>%
  summarise(n = n()) %>%
  mutate(orig.ident = factor(orig.ident, levels = meta$orig.ident)) %>%
  ggplot(aes(x = orig.ident, y = n, fill = orig.ident))+
  geom_col()+
  geom_text(aes(label = scales::comma(n)), vjust = 0, family = "Helvetica", size = 3, nudge_y = 50)+
  scale_y_continuous(label = label_number(suffix = "K", scale = 1e-3))+
  scale_fill_manual(values = meta$col)+
  labs(y = "Number of cells", x = "Tumor ID")+
  theme_proj()+
  theme(axis.line.x = element_line(linewidth = 0.25),
        axis.line.y = element_line(linewidth = 0.25),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
        panel.grid.major =element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 4, width = 4)
print(p)
dev.off()

#===============================================================================
#
# PANEL D - DIFFERENTIAL EXPRESSION
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
         lab = paste0(round(as.numeric(prop)*100, 1), "%")) %>%
  arrange(factor) %>%
  ggplot(aes(x = prop, y = factor)) +
  geom_col(fill = proj_cols("yellow")) +
  geom_text(aes(label = lab), hjust = 0, family = "Helvetica", size = 3)+
  geom_vline(xintercept = 0, linewidth = 1)+
  scale_x_continuous(expand = c(0,0),limits = c(0, 1.5), breaks = c(0, 1), labels = c(0, 1))+
  labs(x = "Proportion of\nexpressing cells")+
  theme_proj()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        plot.margin = margin(t=0.1, r=-0.05, b=0.1, l=0.1, unit = "in"),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 11, margin = margin(-0.3, 0, 0, 0, unit = "in"), family = "Helvetica"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5))

keep_ids <- data@meta.data %>% 
  group_by(orig.ident) %>% 
  dplyr::summarise(n=n()) %>% 
  filter(n >= 100) %>% 
  pull(orig.ident)

res <- data.frame()
for(i in keep_factors){
  
  tmp <- data@meta.data %>% 
    filter(TYPE %in% c("classical", "mesenchymal")) %>%
    filter(orig.ident %in% keep_ids) %>%
    rename_at(all_of(i), ~"factor") %>%
    mutate(scale_UMI = scale(nCount_RNA))
  
  out <- lmer(log2(factor) ~ scale_UMI + TYPE + (1 | orig.ident),
              data = tmp, 
              na.action = na.omit)
  
  coef <- coef(summary(out))
  confint <- confint(out)
  
  log2FC <- logfc(tmp$factor[tmp$TYPE=="mesenchymal"], tmp$factor[tmp$TYPE=="classical"])
  
  rand_eff <- VarCorr(out)$orig.ident[1]
  
  res <- rbind(res, cbind(factor = i,
                          n = nrow(tmp), 
                          log2FC,
                          coef = coef["TYPEmesenchymal", "Estimate"], 
                          se = coef["TYPEmesenchymal", "Std. Error"], 
                          lowerCI = confint["TYPEmesenchymal", "2.5 %"],
                          upperCI = confint["TYPEmesenchymal", "97.5 %"],
                          rand_eff,
                          p = coef["TYPEmesenchymal", "Pr(>|t|)"]))
  
  rm(list = ls()[ls() %in% c("out", "coef", "tmp", "p")]) 
  
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI","p", "rand_eff")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "BH"),
         col = ifelse(coef > 0 & p < 0.05, "Upregulated", "NS (Unadj. *p* > 0.05)"),
         col = ifelse(coef < 0 & p < 0.05, "Downregulated", col),
         col = factor(col, levels = c("Downregulated", "Upregulated", "NS (Unadj. *p* > 0.05)")))

write_xlsx(res, file.path(res_path, "GMB_res.xlsx"))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated", "Upregulated", "NS (Unadj. *p* > 0.05)")

p2 <- ggplot(res, aes(x = log2FC, y = factor, fill = col)) +
  geom_col()+
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5)+
  geom_vline(xintercept = 0.5, lty = 2)+
  geom_vline(xintercept = -0.5, lty = 2)+
  scale_fill_manual(values = cols)+
  labs(x = "Log<sub>2</sub>FC<br>(mesenchymal vs. classical)")+
  guides(fill = guide_legend(title = "Significance", ncol = 1, override.aes = aes(shape = 21, size = 4)))+
  theme_proj() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 11, family = "Helvetica"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, family = "Helvetica"),
        plot.margin = margin(0, r = 0.1, b=0.2, 0, "in"),
        panel.grid.major.x = element_line(linewidth = 0.25, color = "#F0F0F0"),
        panel.grid.major.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5), 
        axis.ticks.y = element_line(linewidth = 0.5), 
        panel.border = element_rect(linewidth = 1, color = "black"),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        legend.title = element_markdown(size = 10, family = "Helvetica", margin = margin(b = 0.05, unit = "in")),
        legend.text = element_markdown(size = 10, family = "Helvetica", margin = margin(l = 0.05, unit = "in")),
        legend.margin = margin(b = -0.2, unit = "in"),
        legend.title.position = "top",
        legend.position = "top",
        legend.direction = "horizontal", 
        legend.key.spacing.y = unit(0.02, "in"))

p <- wrap_plots(p1, p2) + plot_layout(nrow = 1, widths = c(0.4, 0.6))

cairo_pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 6, width = 6)
print(p)
dev.off()

#===============================================================================
#
# PANEL E - CLASS PROPORTIONS
#
#===============================================================================
plots <- list()
for(i in paste0("scHPF_", c(9, 23))){
  
  thresh <- median(data@meta.data[,i]) + 2*mad(data@meta.data[,i])
  
  tmp <- data@meta.data %>% 
    rename_at(i, ~"factor") %>%
    filter(TYPE %in% c("classical", "mesenchymal")) %>%
    filter(orig.ident %in% keep_ids) %>%
    mutate(class = ifelse(factor > thresh, "High", "Low"),
           class = factor(class, levels = c("Low", "High")),
           facet = names(factors)[factors %in% i]) %>%
    mutate(TYPE = str_to_sentence(TYPE)) %>%
    dplyr::select(TYPE, class)
  
  plot_data <- tmp %>% 
    group_by(TYPE, class) %>%
    summarise(n = n()) %>%
    group_by(TYPE) %>% 
    reframe(class, n, prop = n/sum(n)) %>% 
    mutate(lab = paste0(round(prop*100, digits = 0), "%"),
           facet = names(factors)[factors %in% i])
  
  anno <- chisq.test(tmp$TYPE, tmp$class)$p.value
  anno <- ifelse(anno < 0.001, "***", ifelse(anno < 0.01, "**", ifelse(anno < 0.05, "*", "ns")))
  
  p <- ggplot(plot_data, aes(x = TYPE, y = prop, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge())+
    geom_text(aes(label = lab, y = prop + 0.05),
              position = position_dodge(width = 0.9),
              color = 'black',
              vjust = 0,
              size = 3,
              family = "Helvetica") +
    geom_signif(comparisons = list(c("Classical", "Mesenchymal")),
                step_increase = 0.1, 
                tip_length = 0, 
                margin_top = 0.15,
                map_signif_level = T, 
                vjust = 0,
                annotation = anno) +
    facet_wrap(~facet, strip.position = "top") +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 1.1), 
                       breaks = seq(0, 1, 0.25),
                       labels = function(x) formatC(x, format = "g"))+
    scale_fill_manual(name = "Expression\nstatus",
                      values = unname(proj_cols("light blue", "purple"))) +
    guides(fill = guide_legend(nrow = 2))+
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
          plot.margin = margin(t=0.1, r=0, b=0.1, l=0.05, "in"),
          legend.margin = margin(l = -0.1, unit = "in"),
          legend.key.spacing.y = unit(0.01, "in"),
          strip.text = element_text(size = 8, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")),
          strip.background = element_blank())
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + 
  plot_layout(nrow = 2, guides = "collect", axis_titles = "collect", axes = "collect") &
  theme(legend.position = 'bottom',
        legend.direction = "horizontal")

pdf(file.path(res_path, paste0(fig, "panel_E.pdf")), height = 6.5, width = 2.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL F - HEATMAP
#
#===============================================================================
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
for(i in paste0("scHPF_", c(9, 23))){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

averages <- AverageExpression(data,
                              features = genes$name,
                              group.by = "orig.ident",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data
averages <- averages[genes$name, meta$orig.ident]

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

col_split <- factor(meta$TYPE, levels = c("classical", "mesenchymal", "proneural"))

Idents(data) <- "TYPE"
sig <- FindMarkers(data, 
                   ident.1 = "mesenchymal", 
                   ident.2 ="classical",
                   features = unique(keep_highlights$name),
                   assay = "RNA",
                   test.use = "MAST", 
                   latent.vars = c("nCount_RNA", "SEX"),
                   logfc.threshold = 0,
                   min.pct = 0) %>%
  rownames_to_column("gene") %>%
  filter(gene %in% rownames(averages)) %>%
  mutate(gene = factor(gene, levels = rownames(averages)),
         cols = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, proj_cols("red"), 
                       ifelse(avg_log2FC < 0 & p_val_adj < 0.05, proj_cols("blue"), "black")))%>%
  rowwise() %>%
  mutate(at = which(rownames(averages) %in% gene)) %>%
  ungroup() %>%
  arrange(at)

right_ha <- rowAnnotation(labs = anno_mark(at = sig$at, 
                                          labels = sig$gene,
                                          side ="right",
                                          labels_gp = gpar(fontface = "italic", 
                                                           fontfamily = "Helvetica", 
                                                           fontsize = 10,
                                                           col = sig$cols)))

cairo_pdf(file.path(res_path, paste0(fig, "panel_F.pdf")), height = 7, width = 7)
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
             column_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 0),
             column_names_rot = 45,
             column_names_side = "top",
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_gap = unit(0.05, "in"),
             row_gap = unit(0.02, "in"),
             right_annotation = right_ha,
             cluster_columns = T,
             column_dend_side = "bottom",
             column_split = col_split,
             show_row_dend = FALSE,
             column_dend_height = unit(0.3, "in"),
             heatmap_legend_param = list(
               title = expression('Average\nscaled\nexpression\nlevel'), 
               legend_height = unit(1, "in"),
               border = "black"
             ))

draw(p, 
     align_heatmap_legend = "heatmap_center", 
     background = "transparent", 
     padding = unit(c(b=0.1, l=0.2, b=0.3, r=0.1), "in"),
     legend_title_gp = gpar(fontsize = 8, hjust = 0.5, fontfamily = "Helvetica"))

for(i in 1:length(diff_factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}
dev.off()

