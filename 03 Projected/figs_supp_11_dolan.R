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
library(readxl)
library(ggridges)

# PATHS
#-------------------------------------------------------------------------------
dataset <- "Dolan_iPSC"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp11_"

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
cell_embeddings <- data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names")
data@meta.data <- data@meta.data %>% left_join(cell_embeddings)
rownames(data@meta.data) <- data@meta.data$cell_names

tmp <- read_xlsx(file.path(home_path, "proj_Dolan_iPSC/iPSC_H1_Integration_metadata2024.xlsx")) %>% 
  mutate(cell_names = ifelse(line == "H1", paste0(dataset, "_", bc), paste0(line, "_", bc)))

data@meta.data <- data@meta.data %>% left_join(tmp, by = "cell_names")
rownames(data@meta.data) <- data@meta.data$cell_names

h1 <- subset(data, subset = line %in% "H1")
ipsc <- subset(data, subset = line %in% c("iCW50036", "iCW50118", "iCW70347"))

tmp <- read_csv(file.path(home_path, "proj_Dolan_iPSC/iMGL_Figure1_metadata2024.csv"))

h1@meta.data <- h1@meta.data %>% left_join(tmp %>% dplyr::select(bc, final_clusters), by = 'bc')
rownames(h1@meta.data) <- h1@meta.data$cell_names

ipsc[["percent.mt"]] <- PercentageFeatureSet(ipsc, pattern ="^MT-|^MTRN")
h1[["percent.mt"]] <- PercentageFeatureSet(h1, pattern ="^MT-|^MTRN")

#===============================================================================
#
# PANEL A
#
#===============================================================================
hex <- make_hexbin(h1, nbins = 100, dimension_reduction = "schpf.umap")

p1 <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N cells\nper\nhexbin")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(l = -0.2, unit='in'))

treatments <- c("H1_Ctrl", "H1_Ab", "H1_Apop", "H1_Myln", "H1_Syn")
names(treatments) <- c("Control", "Amyloid", "AN", "Myelin", "Synaptosomes")

plots <- list()
for(i in treatments){
  
  tmp_hex <- make_hexbin(subset(h1, subset = line_treatment %in% i), nbins = 50, dimension_reduction = "schpf.umap")
  
  p <- ggplot() + 
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(tmp_hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    ggtitle(paste0(names(treatments)[treatments %in% i], " (", formatC(ncol(tmp_hex), big.mark = ","), ")"))+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(l=0.05, unit = 'in')),
          legend.title = element_blank(),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, 'in'),
          plot.title = element_text(margin = margin(0, 0, 0.01, 0, 'in'), size = 14, 
                                    face = "bold"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0.05, 0.05, 0, 0.05))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(c(list(p1), plots)) + 
  plot_layout(design = c("112
                         113
                         456"))

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 7, width = 7)
print(p)
dev.off()
#===============================================================================
#
# PANEL B - RIDGE PLOT
#
#===============================================================================
annot <- c(1:11)
names(annot) <- c("Other", "Disease-associated", "Antigen-presentation", "Antigen-presentation",
                  "Homeostatic", "Proliferative", "Antigen-presentation","Disease-associated",
                  "Proliferative", "Proliferative", "Interferon-responsive")

annot_cols <- c(`1`="grey", `2`="#fe611b", `3`="#009340", `4`="#009340",`5`= "#2f47a1",`6`= "#f2e300",`7`= "#009340",
                `8`= "#fe611b", `9`="#f2e300",`10`= "#f2e300",`11`= "#d2002f")
names(annot_cols) <- names(annot)

plots <- list()
for(i in factors){
  
  tmp <- h1@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(facet = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i]), 
           final_clusters = factor(final_clusters, levels = annot, labels = names(annot))) %>%
    drop_na(final_clusters)
  
  p <- ggplot(tmp, aes(x = factor, y = final_clusters, fill = final_clusters)) + 
    geom_density_ridges()+
    scale_x_continuous(trans = "log2", labels = function(x) formatC(x, digits = 1,  format = "g")) + 
    geom_vline(xintercept = 0.01, lty = 2) + 
    scale_fill_manual(values = annot_cols)+
    facet_wrap(~facet)+
    labs(x = "log(Cell score)")+
    theme_proj() +
    theme(axis.text = element_text(size = 10),
          axis.title.x =  element_blank(),
          axis.title.y =  element_blank(),
          legend.position = "none",
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"))
  
  if(grepl("-", names(factors)[factors %in% i])){
    p <- p + theme(strip.text = element_markdown(size = 10, margin = margin(b = 0.05, t = 0.05, unit = "in"), family = "Helvetica"))
  } else {
    p <- p + theme(strip.text = element_text(size = 10, margin = margin(b = 0.05, t = 0.05, unit = "in"), family = "Helvetica"))
    
  }
  
  plots <- c(plots, list(p))
}

plots[[21]] <- plots[[21]] + theme(axis.title.x = element_text(size = 14, family = "Helvetica"))
p <- wrap_plots(plots) + plot_layout(guides = "collect", nrow = 4, axes = "collect_y")

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
  
  tmp <- h1@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

keep_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

res <- data.frame()
for(i in keep_factors){
  for(j in 1:11){
    
    tmp <- h1@meta.data %>%
      rename_at(all_of(i), ~"factor") %>%
      mutate(scale_UMI = c(scale(nCount_RNA)), 
             scale_MT = c(scale(percent.mt)),
             cluster = ifelse(final_clusters == j, 1, 0)) %>%
      drop_na(cluster) %>%
      dplyr::select(factor, cluster, scale_UMI, scale_MT, line_replicate)
    
    log2FC <- logfc(tmp$factor[tmp$cluster==1], tmp$factor[tmp$cluster==0])
    
    out <- lmer(log2(factor) ~ cluster + scale_UMI + scale_MT + (1 | line_replicate), data = tmp)
    
    rand_eff <- VarCorr(out)$line_replicate[1]
    sing <- any(grepl("singular", unlist(out@optinfo$conv)))
    
    if(rand_eff == 0 | sing == T){ out <- lm(log2(factor) ~ cluster + scale_UMI + scale_MT, data = tmp) }
    
    coef <- coef(summary(out))
    confint <- confint(out, parm = "cluster")
    
    res <- rbind(res, cbind(factor = i, 
                            cluster = j, 
                            n = nrow(tmp), 
                            log2FC,
                            coef = coef["cluster", "Estimate"], 
                            se = coef["cluster", "Std. Error"], 
                            lowerCI = confint["cluster", "2.5 %"],
                            upperCI = confint["cluster", "97.5 %"],
                            rand_eff,
                            sing,
                            p = coef["cluster", "Pr(>|t|)"]))
    
  }
}

res <- res %>%
  mutate_at(all_of(c("n","log2FC", "coef", "se", "lowerCI", "upperCI", "rand_eff","p")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "bonferroni"))

write_xlsx(res, file.path(res_path, "H1_cluster_exp.xlsx"))

mat <- res %>%
  dplyr::select(factor, cluster, log2FC) %>%
  spread(cluster, log2FC) %>%
  column_to_rownames("factor") %>%
  as.matrix()
mat <- t(mat[,as.character(1:11)])

sig <- res %>%
  mutate(p.adj = ifelse(p.adj == 0, min(p.adj[p.adj!=0]), p.adj),
         p.adj = ifelse(log2FC > 0, -log10(p.adj), log10(p.adj))) %>% 
  dplyr::select(factor, cluster, p.adj) %>%
  spread(cluster, p.adj) %>%
  column_to_rownames("factor") %>%
  as.matrix()
sig <- t(sig[,as.character(1:11)])

col_fun <- colorRamp2(breaks = c(min(sig, na.rm=T), 0, max(sig, na.rm=T)), 
                      colors = c(proj_cols("blue"),"white", proj_cols("red")))

max_size <- 0.5

annot <- c(`1`="grey", `2`="#fe611b", `3`="#009340", `4`="#009340",`5`= "#2f47a1",`6`= "#f2e300",`7`= "#009340",
           `8`= "#fe611b", `9`="#f2e300",`10`= "#f2e300",`11`= "#d2002f")
names(annot) <- c("Other", "Disease-associated", "Antigen-presentation", "Antigen-presentation",
                  "Homeostatic", "Proliferative", "Antigen-presentation","Disease-associated",
                  "Proliferative", "Proliferative", "Interferon-responsive")

ha <- HeatmapAnnotation(bar = names(annot), 
                        col = list(bar = annot),  
                        gp = gpar(col = "white", lwd = 2),
                        show_annotation_name = c(bar = FALSE),
                        annotation_legend_param = list(bar = list(title = "Cluster annotations",
                                                                  labels_gp = gpar(fontsize = 10, 
                                                                                   fontfamily = "Helvetica"),
                                                                  title_gp = gpar(fontsize = 10, 
                                                                                  fontfamily = "Helvetica", 
                                                                                  fontface = "plain"))),
                        which = "row")

set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = col_fun,
             cell_fun = function(j, i, x, y, w, h, fill) {
               grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "white", fill = "white"))
               grid.circle(x=x,y=y,r = abs(mat[i, j]) * unit(max_size, "mm"),
                           gp = gpar(fill = col_fun(sig[i, j]), col = "black"))
             },
             row_names_side = "left",
             row_dend_side = "right",
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_centered = F,
             left_annotation = ha,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(Adj. p)'), 
               legend_width = unit(2, "in"),
               border = "black",
               labels_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
               title_gp = gpar(fontsize = 10, fontfamily = "Helvetica", fontface = "plain")
             )
)

lgd_list <- list(
  ComplexHeatmap::Legend(labels = c(0.5, 1, 2, 3, 4), 
                         title = expression('log'[2]*'FC'),
                         graphics = list(
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(1) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(2) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(3) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(4) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black"))),
                         labels_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
                         title_gp = gpar(fontsize = 10, fontfamily = "Helvetica", fontface = "plain")
  )
)

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 4.5, width = 9.5)
draw(p, 
     align_heatmap_legend = "heatmap_center",
     heatmap_legend_side = "right",
     background = "transparent",
     annotation_legend_list = lgd_list, 
     legend_grouping = "original",
     padding = unit(c(0.2, 0.6, 0.2, 0.2), "in"))

decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1)) })
dev.off()

#===============================================================================
#
# PANEL D
#
#===============================================================================
res <- data.frame()
for(i in keep_factors){
  for(j in treatments[-1]){
    
    tmp <- h1@meta.data %>%
      filter(line_treatment %in% c("H1_Ctrl", j)) %>%
      rename_at(all_of(i), ~"factor") %>%
      mutate(scale_UMI = c(scale(nCount_RNA)), 
             scale_MT = c(scale(percent.mt)),
             treatment = ifelse(line_treatment %in% j, 1, 0)) %>%
      drop_na(treatment) %>%
      dplyr::select(factor, treatment, scale_UMI, scale_MT, line_replicate)
    
    log2FC <- logfc(tmp$factor[tmp$treatment==1], tmp$factor[tmp$treatment==0])
    
    out <- lmer(log2(factor) ~ treatment + scale_UMI + scale_MT + (1 | line_replicate), data = tmp)
    
    rand_eff <- VarCorr(out)$line_replicate[1]
    sing <- any(grepl("singular", unlist(out@optinfo$conv)))
    assign("last.warning", NULL, envir = baseenv())
    
    if(rand_eff == 0 | sing == T){ out <- lm(log2(factor) ~ treatment + scale_UMI + scale_MT, data = tmp) }
    
    coef <- coef(summary(out))
    confint <- confint(out, parm = "treatment")
    
    res <- rbind(res, cbind(factor = i, 
                            treatment = j, 
                            n = nrow(tmp), 
                            log2FC,
                            coef = coef["treatment", "Estimate"], 
                            se = coef["treatment", "Std. Error"], 
                            lowerCI = confint["treatment", "2.5 %"],
                            upperCI = confint["treatment", "97.5 %"],
                            rand_eff,
                            sing,
                            p = coef["treatment", "Pr(>|t|)"]))
    
  }
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI", "rand_eff","p")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = factors, labels = names(factors)),
         p.adj = p.adjust(p, method = "bonferroni"))

write_xlsx(res, file.path(res_path, "H1_treatment.xlsx"))

mat <- res %>%
  dplyr::select(factor, treatment, log2FC) %>%
  spread(factor, log2FC) %>%
  column_to_rownames("treatment") %>%
  as.matrix()
rownames(mat) = names(treatments)[-1]

sig <- res %>%
  mutate(p.adj = ifelse(p.adj ==0, min(p.adj[p.adj!=0]), p.adj),
         p.adj = ifelse(log2FC > 0, -log10(p.adj), log10(p.adj))) %>%
  dplyr::select(factor, treatment, p.adj) %>%
  spread(factor, p.adj) %>%
  column_to_rownames("treatment") %>%
  as.matrix()
rownames(sig) <- names(treatments)[-1]

col_fun <- colorRamp2(breaks = c(min(sig, na.rm=T), 0, max(sig, na.rm=T)), 
                      colors = c(proj_cols("blue"),"white", proj_cols("red")))

max_size <- 1.25

set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = col_fun,
             cell_fun = function(j, i, x, y, w, h, fill) {
               grid.rect(x = x, y = y, width = w, height = h, gp = gpar(col = "white", fill = "white"))
               grid.circle(x=x,y=y,r = abs(mat[i, j]) * unit(max_size, "mm"),
                           gp = gpar(fill = col_fun(sig[i, j]), col = "black"))
             },
             row_names_side = "left",
             row_dend_side = "right",
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(Adj. p)'), 
               legend_width = unit(2, "in"),
               border = "black",
               labels_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
               title_gp = gpar(fontsize = 10, fontfamily = "Helvetica", fontface = "plain")
             )
)

lgd_list <- list(
  ComplexHeatmap::Legend(labels = c(0.1,0.5,1,1.5,2), title = expression('log'[2]*'FC'),
                         graphics = list(
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(1) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(1.5) * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black")),
                           function(x, y, w, h) grid.circle(x = x, y = y, r = 2 * unit(max_size, "mm"),
                                                            gp = gpar(fill = "black"))),
                         labels_gp = gpar(fontsize = 10, fontfamily = "Helvetica"),
                         title_gp = gpar(fontsize = 10, fontfamily = "Helvetica", fontface = "plain")
  )
)

cairo_pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 4, width = 8.5)
draw(p, 
     align_heatmap_legend = "heatmap_center",
     heatmap_legend_side = "right",
     background = "transparent",
     annotation_legend_list = lgd_list, 
     legend_grouping = "original",
     padding = unit(c(0.2, 0.2, 0.2, 0.2), "in"))

decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1)) })
dev.off()

#===============================================================================
#
# PANEL E
#
#===============================================================================
hex <- make_hexbin(ipsc, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N cells\nper\nhexbin")+
  labs(title = "iPSC-derived iMGLs")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(l = -0.2, unit='in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, paste0(fig, "panel_E.pdf")), height = 3.5, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL F
#
#===============================================================================
prop <- data.frame()
for(i in factors){
  
  tmp <- ipsc@meta.data %>% rename_at(all_of(i), ~"factor")
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
  labs(x = "Proportion of\nexpressing cells", y = "Expression Programs")+
  theme_proj()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        plot.margin = margin(t = 0.1, r = -0.05, b = 0.1, l = 0.1, unit = "in"),
        axis.line = element_blank(),
        axis.title.x = element_text(size = 11, margin = margin(b=0.05, 0, 0, 0, unit = "in"), family = "Helvetica"),
        axis.title.y = element_text(size = 12, family = "Helvetica"),
        axis.text.y = element_text(size = 11, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5))

res <- data.frame()
for(i in keep_factors){
  
  tmp <- ipsc@meta.data %>% 
    mutate(treatment = ifelse(treatment == "Ctrl", 0, 1)) %>%
    rename_at(all_of(i), ~"factor") %>%
    mutate(scale_UMI = c(scale(nCount_RNA)), 
           scale_MT = c(scale(percent.mt))) %>%
    dplyr::select(line, factor, treatment, scale_UMI, scale_MT)
  
  log2FC <- logfc(tmp$factor[tmp$treatment==1], tmp$factor[tmp$treatment==0])
  
  out <- lmer(log2(factor) ~ scale_UMI + scale_MT + treatment + (1 | line),
              data = tmp, 
              na.action = na.omit)
  
  coef <- coef(summary(out))
  confint <- confint(out, parm = "treatment")
  rand_eff <- VarCorr(out)$line[1]
  
  res <- rbind(res, cbind(factor = i, 
                          n = nrow(tmp), 
                          log2FC,
                          coef = coef["treatment", "Estimate"], 
                          se = coef["treatment", "Std. Error"], 
                          lowerCI = confint["treatment", "2.5 %"],
                          upperCI = confint["treatment", "97.5 %"],
                          rand_eff,
                          p = coef["treatment", "Pr(>|t|)"]))
  
  rm(list = ls()[ls() %in% c("out", "coef", "tmp", "p")]) 
  
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI","p", "rand_eff")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "BH"),
         col = ifelse(coef > 0 & p < 0.05, "Upregulated", "NS (Adj. *p* > 0.05)"),
         col = ifelse(coef < 0 & p < 0.05, "Downregulated", col),
         col = factor(col, levels = c("Downregulated", "Upregulated", "NS (Adj. *p* > 0.05)")))

write_xlsx(res, file.path(res_path, "iPSC_treatment.xlsx"))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated", "Upregulated", "NS (Adj. *p* > 0.05)")

p2 <- ggplot(res, aes(x = log2FC, y = factor, fill = col)) +
  geom_col()+
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5)+
  geom_vline(xintercept = 0.5, lty = 2, color = "black", linewidth = 0.3)+
  geom_vline(xintercept = -0.5, lty = 2, color = "black", linewidth = 0.3)+
  scale_fill_manual(values = cols)+
  labs(x = "Log<sub>2</sub>FC (Apop vs. Ctrl)")+
  guides(fill = guide_legend(title = "Significance", nrow = 1, override.aes = aes(shape = 21, size = 4)))+
  scale_x_continuous(limits = c(-9, 4), breaks = seq(-9, 4, 3), labels = function(x) formatC(x, format = "g"))+
  theme_proj() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_markdown(size = 11, family = "Helvetica"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11, family = "Helvetica"),
        plot.margin = margin(0, r = 0.1, b=0.2, 0, "in"),
        panel.grid.major = element_line(linewidth = 0.25, color = "#F0F0F0"),
        axis.ticks.x = element_line(linewidth = 0.5), 
        axis.ticks.y = element_line(linewidth = 0.5), 
        panel.border = element_rect(linewidth = 1, color = "black"),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        legend.title = element_markdown(size = 10, family = "Helvetica", margin = margin(b = 0.05, unit = "in")),
        legend.text = element_markdown(size = 10, family = "Helvetica", margin = margin(l = 0.05, unit = "in")),
        legend.title.position = "top",
        legend.position = "top",
        legend.direction = "horizontal", 
        legend.key.spacing.y = unit(0.02, "in"),
        legend.spacing.x = unit(0.1, "in"),
        legend.margin = margin(b = - 0.2, unit = "in"))

p <- wrap_plots(p1, p2) + plot_layout(nrow = 1, widths = c(0.4, 0.6))

cairo_pdf(file.path(res_path, paste0(fig, "panel_F.pdf")), height = 5.5, width = 6)
print(p)
dev.off()

#===============================================================================
#
# PANEL G
#
#===============================================================================
keep_factors <- res %>%
  filter(p.adj < 0.05 & abs(log2FC) > 0.5) %>%
  arrange(col, desc(abs(log2FC))) %>%
  mutate(factor = recode(factor, !!!factors)) %>%
  group_by(col) %>%
  slice_head(n = 3) %>%
  pull(factor)

plots <- list()
for(i in keep_factors){
  
  thresh <- median(ipsc@meta.data[,i]) + 2*mad(ipsc@meta.data[,i])
  
  tmp <- ipsc@meta.data %>% 
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
           facet = gsub(" \\(", "\n\\(", names(factors)[factors %in% i]), 
           facet = gsub('-high', paste0("<sup>high</sup>"), facet))
  
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
                margin_top = 0.1,
                vjust = 0,
                annotation = anno,
                family = "Helvetica", 
                textsize = 4)+
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
          strip.background = element_blank())
  
  if(grepl("-", names(factors)[factors %in% i])){
    p <- p + theme(strip.text = element_markdown(size = 8, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")))
  } else {
    p <- p + theme(strip.text = element_text(size = 8, family = "Helvetica", margin = margin(t=0.05, r=0, b=0.05, l=0, "in")))
  }
  
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

pdf(file.path(res_path, paste0(fig, "panel_G.pdf")), height = 5.5, width = 7)
print(p)
dev.off()
#===============================================================================
#
# PANEL H - iPSC UMAPs
#
#===============================================================================
xmin <- min(ipsc@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(ipsc@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(ipsc@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(ipsc@reductions$schpf.umap@cell.embeddings[,2])

ipsc_treated <- subset(ipsc, subset = treatment == "Apop")
ipsc_untreated <- subset(ipsc, subset = treatment == "Ctrl")

plots <- list()
for(i in paste0("scHPF_", c(21, 15, 5, 20))){
  
  p1 <- Plot_Density_Custom(ipsc_untreated, features = i, reduction = "schpf.umap", pt.size = 0.1)
  p1$data$facet <- names(factors)[factors %in% i]
  p1 <- p1 + facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(name = "Ctrl", limits = c(ymin, ymax))+
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
          strip.text = element_text(size = 10, margin = margin(b = 0.05, t = 0.05, unit = "in"), 
                                    color = ifelse(i %in% paste0("scHPF_", c(21, 15)), proj_cols("blue"), proj_cols("pink"))),
          strip.placement = "inside", 
          strip.background = element_rect(fill = NA, color = NA, linewidth = NA),
          plot.margin  = margin(0, 0, 0, 0, "in"),
          axis.title.y = element_textbox_simple(family = "Helvetica", 
                                                face = "bold", 
                                                size = 12, 
                                                orientation = "left-rotated",
                                                hjust = 1, 
                                                halign = 0.5, 
                                                fill = proj_cols_grey("light grey"), 
                                                padding = margin(2.5, 2.5, 2.5, 2.5),
                                                margin = margin(0, 0, -0.05, 0, unit = "in")))
  
  p2 <- Plot_Density_Custom(ipsc_treated, features = i, reduction = "schpf.umap", pt.size = 0.1)+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(name = "Apop", limits = c(ymin, ymax))+
    theme_umap()+
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.title = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          axis.title.y = element_textbox_simple(family = "Helvetica", 
                                                face = "bold", 
                                                size = 12, 
                                                orientation = "left-rotated",
                                                hjust = 1, 
                                                halign = 0.5, 
                                                color = "white",
                                                fill = proj_cols("purple"), 
                                                padding = margin(2.5, 2.5, 2.5, 2.5),
                                                margin = margin(0, 0, -0.05, 0, unit = "in")))
  
  p <- wrap_plots(p1, p2) +  theme(plot.margin  = margin(0, 0, 0, 0, "in")) + plot_layout(ncol = 1)
  plots <- c(plots, list(p))
  
}

plots[[2]] <- plots[[2]] & theme(axis.title.y = element_blank())
plots[[3]] <- plots[[3]] & theme(axis.title.y = element_blank())
plots[[4]] <- plots[[4]] & theme(axis.title.y = element_blank())

p <- wrap_plots(plots) + plot_layout(nrow = 1)

pdf(file.path(res_path,  paste0(fig, "panel_G.pdf")), height = 3.5, width = 8)
print(p)
dev.off()

#===============================================================================
#
# PANEL I
#
#===============================================================================
data <- NormalizeData(ipsc, normalization.method = "LogNormalize", assay = "RNA") %>%
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

keep_factors <- res %>%
  filter(p.adj < 0.05 & abs(log2FC) > 2) %>%
  arrange(col, desc(abs(coef))) %>%
  mutate(factor = recode(factor, !!!factors)) %>%
  pull(factor)

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
  filter(gene %in% colnames(averages)) %>%
  mutate(gene = factor(gene, levels = colnames(averages)),
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

cairo_pdf(file.path(res_path, paste0(fig, "panel_H.pdf")), height = 5, width = 9)
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
               title_position = "lefttop",
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
