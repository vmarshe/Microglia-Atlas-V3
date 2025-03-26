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
dataset <- "kamp"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp14_"

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

p1 <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  labs(title = "iTF-Microglia<br>(CROP-seq)")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(l=-0.2, unit='in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"))

pdf(file.path(res_path, paste0(fig, "panel_A1.pdf")), height = 3.5, width = 3.5)
print(p1)
dev.off()

plot_list <- list()
for(i in c("NTC", "ATP1A1", "CDK12", "MED1", "STK40")){
  hex <- make_hexbin(subset(data, subset = gene.name %in% i), nbins = 30, dimension_reduction = "schpf.umap")
  p <- ggplot() + 
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    ggtitle(paste0("<i>", i, "</i> (", nrow(hex@meta.data), ")"))+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
          legend.title = element_blank(),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, 'in'),
          plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 14),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0.05, 0.05, 0.05, 0.05))
  plot_list <- c(plot_list, list(p))
}

p <- wrap_plots(c(list(p1 + theme(plot.title = element_blank())), plot_list)) + 
  plot_layout(
    design = 
      "111122
       111133
       445566")

pdf(file.path(res_path, paste0(fig, "panel_A2.pdf")), height = 5, width = 6)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
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
        axis.title.x = element_text(size = 11, margin = margin(b = 0.1, 0, 0, 0, unit = "in"), family = "Helvetica"),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(size = 11, family = "Helvetica"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5))

cairo_pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 5, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
keep_genes <- data@meta.data %>% 
  group_by(gene.name) %>% 
  summarise(n = n()) %>% 
  filter(n > 100) %>% 
  pull(gene.name)
keep_genes <- keep_genes[keep_genes !="NTC"]

res <- data.frame()
for(i in keep_factors){
  for(j in keep_genes){
    
    tmp <- data@meta.data %>% filter(gene.name %in% c("NTC",j ))
    guides <- unique(tmp$zscoreAssignment[!grepl("NTC", tmp$zscoreAssignment)])
    
    tmp <- tmp %>% 
      filter(gene.name %in% c("NTC", j)) %>% 
      mutate(condition = ifelse(gene.name == "NTC", 0, 1),
             scale_UMI = scale(nCount_RNA)) %>%
      rename_at(all_of(i), ~"factor")
    
    out <- lmer(log2(factor) ~ scale_UMI + condition + (1 | sgRNA.name),
                data = tmp, 
                na.action = na.omit)
    
    rand_eff <- VarCorr(out)$sgRNA.name[1]
    sing <- any(grepl("singular", unlist(out@optinfo$conv)))
    
    log2FC <- logfc(tmp$factor[tmp$condition==1], tmp$factor[tmp$condition==0])
    
    if(rand_eff == 0 | sing == T){ out <- lm(log2(factor) ~ scale_UMI + condition, data = tmp) }
    
    coef <- coef(summary(out))
    confint <- confint(out)
    
    res <- rbind(res, cbind(factor = i,
                            gene = j,
                            n = nrow(tmp), 
                            log2FC,
                            coef = coef["condition", "Estimate"], 
                            se = coef["condition", "Std. Error"], 
                            lowerCI = confint["condition", "2.5 %"],
                            upperCI = confint["condition", "97.5 %"],
                            rand_eff, sing,
                            p = coef["condition", "Pr(>|t|)"]))
    
    rm(list = ls()[ls() %in% c("out", "coef", "tmp", "p")]) 
    
  }
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI", "p")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = rev(factors), labels = rev(names(factors))),
         p.adj = p.adjust(p, method = "bonferroni"),
         col = ifelse(coef > 0 & p.adj < 0.05, "Upregulated", "NS (Unadj. *p* > 0.05)"),
         col = ifelse(coef < 0 & p.adj < 0.05, "Downregulated", col),
         col = factor(col, levels = c("Downregulated", "Upregulated", "NS (Unadj. *p* > 0.05)")))

write_xlsx(res, file.path(res_path, "kamp_res.xlsx"))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated", "Upregulated", "NS (Unadj. *p* > 0.05)")

p <- res %>%
  mutate(lab = ifelse(col %in% c("Downregulated", "Upregulated", "FC < 0.5"), 
                      paste0(gene, " (", gsub(" \\([0-9]+\\)", "", factor), ")"), "")) %>%
  ggplot(aes(x = log2FC, y = -log10(p.adj), color =  col)) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_vline(xintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_point(alpha = 0.6)+
  geom_vline(xintercept = 0, col = "black")+
  scale_color_manual(values = cols)+
  scale_size(range = c(0.5, 4))+
  guides(size = guide_legend(title = "LMM<br>coefficient", override.aes = aes(shape = 19)),
         color = guide_legend(title = "Significance", override.aes = aes(shape = 19, size = 4)))+
  scale_x_continuous(labels = function(x) formatC(x, format = "g"))+
  scale_y_continuous(labels = function(x) formatC(x, format = "g"))+
  geom_label_repel(aes(label = lab), family = "Helvetica", seed=1234, size = 3, show.legend = F) +
  labs(x = "log<sub>2</sub> fold-change (sgRNA vs. NTC)",
       y = "log<sub>10</sub>(Adj. <i>p</i>-value)")+
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

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 5.5, width = 7)
print(p)
dev.off()

#===============================================================================
#
# PANEL D - RIDGE PLOT
#
#===============================================================================
id_cols <- c(proj_cols_grey("med grey"), 
             proj_cols( "pink","teal", "purple", "yellow"))
names(id_cols) <- c("NTC", "ATP1A1", "CDK12", "MED1", "STK40")

sig_facs <- res %>% 
  filter(gene %in% c("MED1", "CDK12", "STK40", "ATP1A1") & p.adj < 0.05) %>% 
  mutate(factor = recode(factor, !!!factors)) %>% 
  arrange(col, desc(abs(log2FC))) %>%
  group_by(gene, col) %>%
  slice_head(n = 1) %>%
  pull(factor)
sig_facs <- unique(sig_facs)

plots <- list()
for(i in factors[factors %in% sig_facs]){
  
  tmp <- data@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    filter(gene.name %in% c("MED1", "CDK12", "STK40", "ATP1A1", "NTC")) %>%
    mutate(facet = names(factors)[factors %in% i],
           facet = gsub('-high', paste0("<sup>high</sup>"), facet),
           gene.name = factor(gene.name, levels = names(id_cols)))
  
  p <- ggplot(tmp, aes(y = factor, x = gene.name, fill = gene.name)) + 
    geom_violin(alpha = 0.7) + 
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", color = "black", width = 0.7, lwd = 0.5)+
    geom_signif(comparisons = list(c("NTC", "ATP1A1"), c("NTC", "CDK12"), c("NTC", "MED1"), c("NTC", "STK40")),
                tip_length = 0, margin_top = 0.05, step_increase = 0.05, textsize = 3,
                map_signif_level =c("***"=0.001,"**"=0.01, "*"=0.05, " "=2), vjust = 0.5, family = "Helvetica")+
    scale_y_continuous(labels = function(x) formatC(x, digits = 2, format = "g"), trans = "log2")+
    scale_fill_manual(values = id_cols)+
    facet_wrap(~facet)+
    labs(y = "log(Cell score)", x = "Density")+
    theme_proj() +
    theme(strip.text = element_markdown(size = 10, margin = margin(b = 0.05, t = 0.05, unit = "in"), family = "Helvetica"),
          axis.text.x = element_text(size = 10, family = "Helvetica", hjust = 1, angle = 45),
          axis.text.y = element_text(size = 10, family = "Helvetica"),
          axis.title.y = element_text(size = 14, family = "Helvetica"),
          axis.title.x =  element_blank(),
          legend.position = "none",
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          panel.grid.major = element_blank())
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(guides = "collect", nrow = 1, axes = "collect", axis_titles = "collect")

cairo_pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 3, width = 14)
print(p)
dev.off()
#===============================================================================
#
# PANEL D - HEATMAP
#
#===============================================================================
tmp <- subset(data, subset = gene.name %in% c("NTC", "MED1", "CDK12", "STK40", "ATP1A1")) %>%
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

sig_facs <- res %>% 
  filter(gene %in% c("MED1", "CDK12", "STK40", "ATP1A1") & abs(log2FC) > 0.5 & p.adj < 0.05) %>% 
  mutate(factor = recode(factor, !!!factors)) %>% 
  group_by(gene, col) %>%
  arrange(gene, p.adj) %>%
  pull(factor)
sig_facs <- unique(sig_facs)

keep_highlights <- data.frame()
for(i in unique(sig_facs)){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

averages <- AverageExpression(tmp,
                              features = genes$name,
                              group.by = "zscoreAssignment",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data
averages <- averages[genes$name,]
averages <- t(averages)
averages <- averages[c("NTC_6300", "NTC_6301", "NTC_6302", "NTC_6303", 
                       sort(unique(grep("ATP1A1", tmp$zscoreAssignment, v=T))),
                       sort(unique(grep("CDK12", tmp$zscoreAssignment, v=T))),
                       sort(unique(grep("MED1", tmp$zscoreAssignment, v=T))),
                       sort(unique(grep("STK40", tmp$zscoreAssignment, v=T)))),]

id_cols <- c(proj_cols_grey(), 
             proj_cols( "pink","teal", "purple", "yellow"),
             lighten(proj_cols("pink", "teal", "purple", "yellow"), space ="HLS", amount = 0.5))
id_cols <- id_cols[c(1:4, 5, 9, 6, 10, 7, 11, 8, 12)]
names(id_cols) <- rownames(averages)

ha <- rowAnnotation(mod = anno_simple(names(id_cols), 
                                      col = id_cols, 
                                      height = unit(0.15, "in"), 
                                      gp = gpar(col = "white", lwd = 1.5)),
                    annotation_name_gp = gpar(fontsize = 0),
                    show_legend = F,
                    border = F)

row_split <- factor(c(rep("NTC", 4), rep("ATP1A1", 2), rep("CDK12", 2), rep("MED1", 2), rep("STK40", 2)), 
                    levels = c("NTC", "ATP1A1", "CDK12", "MED1", "STK40"))

col_split <- genes %>% filter(name %in% colnames(averages)) %>% pull(factor_name)
col_split <- factor(col_split, levels = unique(genes$factor_name))

bottom_ha <- HeatmapAnnotation(labs = anno_mark(at = match(keep_highlights$name, colnames(averages)), 
                                                labels = keep_highlights$name,
                                                side ="right",
                                                labels_gp = gpar(fontface = "italic", 
                                                                 fontfamily = "Helvetica", 
                                                                 fontsize = 10)),
                               which = "column")

pdf(file.path(res_path, paste0(fig, "panel_E.pdf")), height = 5, width = 14)
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
