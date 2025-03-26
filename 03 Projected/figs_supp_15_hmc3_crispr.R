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
library(patchwork)
library(schex)
library(Seurat)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
dataset <- "HMC3crispr"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp15_"

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
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(0, 0, 0, -0.1, 'in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 3, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
target_genes <- data@meta.data %>% 
  group_by(gene_target) %>% 
  summarise(KO = sum(mixscape_class.global == "KO"), NP = sum(mixscape_class.global == "NP")) %>% 
  filter(KO >=30 & NP >= 30) %>% 
  pull(gene_target)

prop <- data.frame()
for(i in factors){
  
  tmp <- data@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

keep_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

res <- data.frame()
for(j in target_genes){
  for(i in keep_factors){
    
    tmp <- data@meta.data %>% 
      filter(gene_target == j) %>%
      rename_at(all_of(i), ~"factor") %>%
      mutate(condition = ifelse(mixscape_class.global == "KO", 1, 0), 
             scale_UMI = c(scale(nCount_RNA))) %>%
      dplyr::select(condition, factor, scale_UMI, gene_gRNA)
    
    out <- lmer(log2(factor) ~ scale_UMI + condition + (1 | gene_gRNA),
                data = tmp, 
                na.action = na.omit)
    
    rand_eff <- VarCorr(out)$gene_gRNA[1]
    sing <- any(grepl("singular", unlist(out@optinfo$conv)))
    
    if(rand_eff == 0 | sing == T){ out <- lm(log2(factor) ~ scale_UMI + condition, data = tmp) }
    
    coef <- coef(summary(out))
    confint <- confint(out)
    
    log2FC <- logfc(tmp$factor[tmp$condition==1], tmp$factor[tmp$condition==0])
    
    res <- rbind(res, cbind(factor = i,
                            gene = j,
                            n = nrow(tmp), 
                            log2FC,
                            coef = coef["condition", "Estimate"], 
                            se = coef["condition", "Std. Error"], 
                            lowerCI = confint["condition", "2.5 %"],
                            upperCI = confint["condition", "97.5 %"],
                            p = coef["condition", "Pr(>|t|)"] ))
    
  }
  
  rm(list = ls()[ls() %in% c("out", "d", "tmp")])
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI","p")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "bonferroni"),
         col = ifelse(log2FC > 0 & p.adj < 0.05, "Upregulated", "NS (Adj. *p* > 0.05)"),
         col = ifelse(log2FC < 0 & p.adj < 0.05, "Downregulated", col),
         col = factor(col, levels = c("Downregulated", "Upregulated", "NS (Adj. *p* > 0.05)")))

write_xlsx(res %>% dplyr::select(factor:p.adj), file.path(res_path, "HMC3crispr.xlsx"))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated", "Upregulated", "NS (Adj. *p* > 0.05)")

p <- res %>%
  mutate(lab = ifelse(col %in% c("Downregulated", "Upregulated"), 
                      paste0(gene, " (", gsub(" \\([0-9]+\\)", "", factor), ")"), "")) %>%
  ggplot(aes(x = log2FC, y = -log10(p.adj), color =  col)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_point(alpha = 0.6)+
  geom_vline(xintercept = 0, col = "black")+
  scale_color_manual(values = cols)+
  guides(color = guide_legend(title = "Significance", override.aes = aes(shape = 19, size = 4)))+
  scale_x_continuous(labels = function(x) formatC(x, format = "g"))+
  scale_y_continuous(labels = function(x) formatC(x, format = "g"))+
  geom_text_repel(aes(label = lab), family = "Helvetica", seed = 1234, size = 3, show.legend = F) +
  labs(x = "Log<sub>2</sub>FC (KO vs. NP)",
       y = "-log<sub>10</sub>(Adj. *p*-value)")+
  theme_proj()+
  theme(panel.border = element_rect(linewidth = 1.5, color = "black"),
        panel.grid.major = element_blank(),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        axis.title.y = element_markdown(size = 14, family = "Helvetica"),
        axis.title.x = element_markdown(size = 14, family = "Helvetica"),
        legend.title = element_markdown(size = 12, family = "Helvetica"),
        legend.text = element_markdown(size = 12, family = "Helvetica", margin = margin(l = 0.05, unit = "in")),
        legend.margin = margin(l = -0.2, unit = "in"),
        legend.key.spacing.y = unit(0.05, "in"))

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 4, width = 6.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
opts <- data.frame(factor = c("scHPF_8", "scHPF_10", "scHPF_21", "scHPF_24"),
                   gene = c("ABI3", "ABI3", "RIN3", "RIN3"), 
                   cols = proj_cols("teal", "teal", "orange2", "orange2"))
opts$factor_lab <- factor(opts$factor, levels = factors, labels = names(factors))

plots <- list()
for(i in 1:nrow(opts)){
  
  tmp <- data@meta.data %>% 
    filter(gene_target %in% c(opts$gene[i], "NTC")) %>%
    mutate(condition = factor(mixscape_class.global, 
                              levels = c("NTC","NP", "KO"), 
                              labels = c("NTC",paste0("*", opts$gene[i], "*-NP"), paste0("*", opts$gene[i], "*-KO")))) %>%
    rename_at(all_of(opts$factor[i]), ~"factor")
  
  tmp <- tmp %>%
    group_by(condition) %>% 
    summarise(n = n()) %>% 
    mutate(lab = paste0(condition, "<br>(*n*=", n, ")")) %>%
    right_join(tmp)
  
  labs <- c(unique(grep("NTC", tmp$lab, v=T)), unique(grep("NP", tmp$lab, v=T)), unique(tmp$lab[!grepl("NP|NTC", tmp$lab)]))
  cols <- unname(c(proj_cols_grey("med grey"), opts$cols[i], lighten(opts$cols[i],  space ="HLS", amount = 0.5)))
  names(cols) <- labs
  tmp$lab <- factor(tmp$lab, levels = labs)
  
  p <- tmp %>%
    mutate(facet = opts$factor_lab[i],
           facet = gsub('-high', paste0("<sup>high</sup>"), facet)) %>%
    ggplot(aes(x = lab, y = factor, fill = lab)) +
    geom_violin(alpha = 0.7) + 
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", color = "black", width = 0.7, lwd = 0.5)+
    geom_signif(comparisons = list(c(labs[2], labs[3])),
                tip_length = 0, margin_top = 0.05,
                map_signif_level = T, vjust = 0.5, family = "Helvetica")+
    scale_y_continuous(labels = function(x) formatC(x, digits = 2, format = "g"),
                       trans = "log2")+
    labs(y = "Cell score", x = "Condition")+
    scale_fill_manual(values = cols) +
    facet_wrap(~facet) +
    theme_proj() +
    theme(axis.line.x = element_line(linewidth = 0.2),
          axis.line.y = element_line(linewidth = 0.2),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
          panel.grid.major = element_blank(),
          axis.text.x = element_markdown(size = 12),
          axis.text.y = element_markdown(size = 12),
          legend.position = "none",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 11, margin = margin(t = 0.05, b = 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
  
}
p <- wrap_plots(plots) + plot_layout(nrow = 2, byrow = F, axis_titles = "collect", axes = "collect_x") &
  theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in"))

pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 5, width = 7)
print(p)
dev.off()

#===============================================================================
#
# PANEL D - HEATMAP
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

tmp <- subset(data, subset = gene_target %in% c("NTC","RIN3", "ABI3")) %>%
  NormalizeData(normalization.method = "LogNormalize", assay = "RNA") %>%
  ScaleData(assay = "RNA")

genes <- data.frame()
for(i in keep_factors){
  genes <- rbind(genes, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:10],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}
exp_genes <- names(which(rowSums(data@assays$RNA@data)[genes$name] > 0))
genes <- genes %>% filter(name %in% exp_genes)

sig_facs <-  res %>%
  filter(gene %in% c("ABI3", "RIN3") & p.adj < 0.05) %>%
  arrange(col, desc(abs(coef))) %>%
  mutate(factor = recode(factor, !!!factors)) %>%
  pull(factor)
sig_facs <- unique(sig_facs)

keep_highlights <- data.frame()
for(i in unique(sig_facs)){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:10],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

averages <- AverageExpression(tmp,
                              features = genes$name,
                              group.by = "mixscape_class",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data
averages <- averages[genes$name,]
averages <- t(averages)
averages <- averages[c("NTC","ABI3 NP", "ABI3 KO", "RIN3 NP", "RIN3 KO"), ] 

id_cols <- c(proj_cols_grey("med grey"), 
             proj_cols( "teal", "orange2"),
             lighten(proj_cols("teal", "orange2"), space ="HLS", amount = 0.5))
id_cols <- id_cols[c(1, 2, 4, 3, 5)]
names(id_cols) <- rownames(averages)

ha <- rowAnnotation(mod = anno_simple(names(id_cols), 
                                      col = id_cols, 
                                      height = unit(0.15, "in"), 
                                      gp = gpar(col = "white", lwd = 1.5)),
                    annotation_name_gp = gpar(fontsize = 0),
                    show_legend = F,
                    border = F)

row_split <- factor(c("NTC","ABI3 NP", "ABI3 KO", "RIN3 NP", "RIN3 KO"), 
                    levels = c("NTC","ABI3 NP", "ABI3 KO", "RIN3 NP", "RIN3 KO"))

averages <- averages[, which(colSums(averages) != 0)]
col_split <- genes %>% filter(name %in% colnames(averages)) %>% pull(factor_name)
col_split <- factor(col_split, levels = unique(genes$factor_name))

bottom_ha <- HeatmapAnnotation(labs = anno_mark(at = match(keep_highlights$name, colnames(averages)), 
                                                labels = keep_highlights$name,
                                                side ="right",
                                                labels_gp = gpar(fontface = "italic", 
                                                                 fontfamily = "Helvetica", 
                                                                 fontsize = 10)),
                               which = "column")

pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 4, width = 14)
set.seed(1234)
p <- Heatmap(averages, name = "mat",
             col = colorRamp2(breaks = c(min(averages), 0, max(averages)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             left_annotation = ha,
             top_annotation = columnAnnotation(foo = anno_empty(border = FALSE, height = unit(0.1, "in"))),
             cluster_rows = F,
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
               border = "black",
               title_gp = gpar(size = 12, fontfamily = "Helvetica"),
               labels_gp = gpar(size = 12, fontfamily = "Helvetica")
             ))

draw(p, 
     heatmap_legend_side = "right", 
     background = "transparent", 
     padding = unit(c(0.2, 0.2, .2, 0.2), "in"),
     align_heatmap_legend = "heatmap_center")

for(i in 1:length(keep_factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, height = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}

dev.off()