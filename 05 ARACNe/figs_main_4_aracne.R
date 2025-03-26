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
library(fgsea)
library(clusterProfiler); options(enrichment_force_universe = FALSE)
library(org.Hs.eg.db)
library(ggtext)
library(grid)
library(metafor)
library(patchwork)
library(readxl)
library(rrvgo)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(writexl)
library(cowplot)
library(VennDiagram)
library(vegan)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sc")
fig <- "fig_4"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
aracne <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
regulators <- unique(aracne$Regulator)
schpf <- readRDS(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

#===============================================================================
#
# PANEL A
#
#===============================================================================
res <- read_xlsx(file.path(res_path, "tab_aracne_scHPF_enrichment.xlsx"))
top_regs <- res %>% 
  group_by(ID) %>% 
  top_n(-3, wt = p.adj.bonf) %>% 
  pull(signature)
top_regs <- unique(top_regs)

mat <- res %>% 
  filter(signature %in% top_regs) %>%
  dplyr::select(ID, p.adj.bonf, signature) %>%
  mutate(p.adj.bonf = -log10(p.adj.bonf),
         ID = factor(ID, levels = factors, labels = names(factors))) %>%
  spread(signature, p.adj.bonf) %>%
  column_to_rownames("ID") %>%
  as.matrix()

cairo_pdf(file.path(res_path, paste0(fig, "_panel_A.pdf")), height = 6, width = 16, bg = "transparent")
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(0, max(mat, na.rm = T)),  colors = c("white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & (mat[i, j] >-log10(0.05) | mat[i, j] <log10(0.05))) {
                 grid.points(pch = 8, x, y, size = unit(0.075, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             cluster_rows = F,
             row_title_side = "left",
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 16),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 14),
             row_names_side = "left",
             row_dend_side = "right",
             column_title = "sc-ARACNE TF Regulons",
             column_title_side = "bottom",
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 16),
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 14),
             column_names_rot = 45,
             cluster_columns = F,
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(Adj. '*italic(p)*'-value)'), 
               legend_width = unit(1.5, "in"),
               direction = "horizontal",
               title_position = "lefttop",
               border = "black",
               labels_gp = gpar(font = 14, fontfamily = "Helvetica"),
               title_gp = gpar(font = 14,  fontfamily = "Helvetica")
             ))

draw(p, 
     align_heatmap_legend = "heatmap_center", 
     heatmap_legend_side = "top",
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1.5)) })
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
mat_sn <- read_xlsx(file.path(home_path, "aracne_sn/tab_aracne_bulkRNA_fgsea.xlsx")) %>% 
  mutate(log10padj = -log10(p.adj.bonf),
         direction = ifelse(NES > 0, "up", "down"),
         signature = factor(signature, 
                            levels = c("ad_reagan_1_vs_0", "amyloid_sqrt","tangles", "cogng_random_slope_adjEduc"),
                            labels = c("Pathological AD", "Amyloid", "Tangles","Cognitive decline"))) %>%
  mutate(log10padj = ifelse(NES > 0, log10padj, -log10padj)) %>% 
  filter(pathway %in% c("ARID5B", "MITF", "PPARG", "CEBPA", "IRF7")) %>%
  dplyr::select(pathway, signature, log10padj) %>%
  arrange(pathway) %>%
  spread(signature, log10padj) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

mat_sc <- read_xlsx(file.path(home_path, "aracne_sc/tab_aracne_bulkRNA_fgsea.xlsx")) %>% 
  filter(pathway %in% c("ARID5B", "MITF", "PPARG", "CEBPA", "IRF7")) %>%
  mutate(p.adj.bonf = p.adjust(pval, method = "bonferroni"), 
         log10padj = -log10(p.adj.bonf),
         direction = ifelse(NES > 0, "up", "down"),
         signature = factor(signature, 
                            levels = c("ad_reagan_1_vs_0", "amyloid_sqrt","tangles", "cogng_random_slope_adjEduc"),
                            labels = c("Pathological AD", "Amyloid", "Tangles","Cognitive decline"))) %>%
  mutate(log10padj = ifelse(NES > 0, log10padj, -log10padj)) %>% 
  dplyr::select(pathway, signature, log10padj) %>%
  arrange(pathway) %>%
  spread(signature, log10padj) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

mat <- rbind(mat_sc, mat_sn)
row_split <- factor(c(rep("sc", 5), rep("sn", 5)), 
                    levels = c("sc", "sn"), 
                    labels = c("sc-ARACNE", "sn-ARACNE"))

pdf(file.path(res_path, paste0(fig, "_panel_B.pdf")), height = 4, width = 4.5)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(min(mat, na.rm=T), 0, max(mat, na.rm=T)), 
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & (mat[i, j] >-log10(0.05) | mat[i, j] < log10(0.05))) {
                 grid.points(pch = 8, x, y, size = unit(0.075, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             cluster_rows = T,
             row_title_side = "right",
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 14),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_names_side = "left",
             row_dend_side = "left",
             row_split = row_split,
             right_annotation =  rowAnnotation(foo = anno_empty(border = FALSE, width = unit(0.2, "in"))),
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 14),
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_rot = 45,
             cluster_columns = T,
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(Adj. '*italic(p)*'-value)'), 
               legend_height = unit(1, "in"),
               #direction = "horizontal",
               border = "black",
               labels_gp = gpar(font = 12),
               title_gp = gpar(font = 12)
             ))

draw(p, 
     heatmap_legend_side = "right", 
     align_heatmap_legend = "heatmap_center", 
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))

for(i in 1:2){
  decorate_heatmap_body("mat", slice = i, { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1.5)) })
  decorate_annotation("foo", slice = i, { grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left") })
}
dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
entrez_names <- read_rds(file.path(home_path, "data/entrez_mapping.rds"))

res <- data.frame()
for(regulator in c("PPARG", "ARID5B", "MITF", "CEBPA")){
  
  genes <- aracne %>% filter(Regulator == regulator & tfmode > 0) %>% pull(Target)
  genes <- entrez_names %>% filter(gene %in% genes) %>% pull(entrez)
  
  out <- enrichGO(gene = genes,
                  universe = as.character(na.omit(entrez_names$entrez)),
                  OrgDb  = org.Hs.eg.db,
                  ont  = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)@result %>%
    as.data.frame() %>% 
    filter(p.adjust < 0.05 & qvalue < 0.05)
  
  simMatrix <- calculateSimMatrix(out$ID, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
  scores <- setNames(-log10(out$qvalue), out$ID)
  reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
  
  pdf(file.path(res_path, paste0(fig, "_panel_C_", regulator, "_treemap.pdf")), height = 3, width = 3)
  treemapPlot(reducedTerms)
  dev.off()
  
}

#===============================================================================
#
# PANEL D
#
#===============================================================================
regulators <- c("IRF7", "ARID5B", "CEBPA", "MITF", "PPARG")

# ROSMAP METADATA
#-------------------------------------------------------------------------------
rosmap_meta <- read_delim(file.path(home_path, "data/dataset_707_basic_11-01-2023.txt"), delim="\t") %>% 
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  arrange(projid) %>%
  dplyr::select(projid, msex, age_death, educ, pmi, study) %>%
  distinct()

multiome_ids <- read_delim(file.path(home_path, "data/dataset_basic_n3638.txt")) %>%
  dplyr::select(projid, wgs_SampleID) %>%
  mutate(projid = str_pad(projid, 8, pad = "0")) %>%
  distinct()

non_overlapping_donors <- readRDS(file.path(home_path, "data/snuc_non_overlap_donors.rds"))

# CUIMC1
#-------------------------------------------------------------------------------
snuc_full <- readRDS(file.path(home_path, "proj_snuc/projection/snuc.rds"))
snuc_full@meta.data$cell_names <- rownames(snuc_full@meta.data)
snuc_full <- NormalizeData(snuc_full) %>% ScaleData(features = regulators)

snuc <- readRDS(file.path(home_path, "proj_snuc_validation/snuc/snuc.rds"))
snuc@meta.data$cell_names <- rownames(snuc@meta.data)
snuc_full@meta.data <- snuc_full@meta.data %>% 
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  left_join(snuc@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(snuc_full@assays$RNA@data[regulators, ] %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(rosmap_meta, by = "projid")
rownames(snuc_full@meta.data) <- snuc_full@meta.data$cell_names
snuc <- snuc_full
rm(snuc_full); gc()

# MIT
#-------------------------------------------------------------------------------
kellis_full <- readRDS(file.path(home_path, "proj_kellis/projection/kellis.rds"))
kellis_full@meta.data$cell_names <- rownames(kellis_full@meta.data)
kellis_full@meta.data$orig.ident <- as.character(kellis_full@meta.data$orig.ident)
kellis_full <- NormalizeData(kellis_full) %>% ScaleData(features = regulators)

kellis <- readRDS(file.path(home_path, "proj_kellis_validation/kellis/kellis.rds"))
kellis@meta.data$cell_names <- rownames(kellis@meta.data)
kellis_full@meta.data <- kellis_full@meta.data %>% 
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  left_join(kellis@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(kellis_full@assays$RNA@data[regulators, ] %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(rosmap_meta, by = "projid")
rownames(kellis_full@meta.data) <- kellis_full@meta.data$cell_names
kellis <- kellis_full
rm(kellis_full); gc()

# CUIMC2
#-------------------------------------------------------------------------------
multiome_full <- readRDS(file.path(home_path, "proj_multiome/projection/multiome.rds"))
multiome_full@meta.data$cell_names <- rownames(multiome_full@meta.data)
multiome_full@meta.data <- multiome_full@meta.data %>% 
  left_join(multiome_ids, by = c("orig.ident"="wgs_SampleID")) %>%
  mutate(orig.ident = as.character(str_pad(projid, 8, pad = "0")))
multiome_full <- NormalizeData(multiome_full) %>% ScaleData(features = regulators)

multiome <- readRDS(file.path(home_path, "proj_multiome_validation/multiome/multiome.rds"))
multiome@meta.data$cell_names <- rownames(multiome@meta.data)
multiome_full@meta.data <- multiome_full@meta.data %>% 
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  left_join(multiome@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(multiome_full@assays$RNA@data[regulators, ] %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("cell_names")) %>%
  left_join(rosmap_meta, by = "projid")

rownames(multiome_full@meta.data) <- multiome_full@meta.data$cell_names
multiome <- multiome_full
rm(multiome_full); gc()

# Keep unique IDs
#-------------------------------------------------------------------------------
kellis_ids <- non_overlapping_donors %>% filter(reference == "Kellis-Tsai") %>% pull(projid) #132
multiome_ids <- non_overlapping_donors %>% filter(reference == "Multiome") %>% pull(projid) #219

kellis <- subset(kellis, subset = orig.ident %in% kellis_ids) # 130
multiome <- subset(multiome, subset = orig.ident %in% multiome_ids) #212

# ASSOCIATIONS
#-------------------------------------------------------------------------------
combos <- data.frame(regulators, factor = c("scHPF_20", rep("scHPF_26", 4)))

res <- data.frame()
for(j in c("snuc", "kellis", "multiome")) { 
  
  data <- eval(parse(text = j))
  
  for(i in 1:5){
    
    counts <- data@meta.data %>% 
      mutate(reg = !!sym(combos$regulators[i])) %>%
      rename_at(all_of(combos$factor[i]), ~"factor") %>% 
      group_by(projid) %>% 
      summarise(factor = mean(factor), 
                reg = mean(reg), 
                umi = mean(nCount_RNA), 
                pmi = unique(pmi), 
                age = unique(age_death), 
                sex = unique(msex),
                study = unique(study))
    
    out <- lm(factor ~ reg + umi + pmi + age + sex + study, data = counts, na.action = na.omit)
    sum <- coefficients(summary(out))
    
    eff <- cbind.data.frame(data = j,
                            reg = combos$regulators[i], 
                            n = nrow(counts),
                            coef = sum["reg", "Estimate"], 
                            se = sum["reg", "Std. Error"], 
                            p = sum["reg", "Pr(>|t|)"])
    
    res <- rbind(res, eff)
    
  }
}

# META-ANALYSIS
#-------------------------------------------------------------------------------
meta_res <- data.frame()
for(i in regulators){
  
  eff <- res %>% filter(data %in% c("snuc", "kellis", "multiome") & reg == i)
  
  meta <- rma(yi = eff$coef, sei = eff$se, weights = eff$n, data = eff)
  
  meta_res <- rbind(meta_res, cbind(reg = i,
                                    meta.se = c(meta$se),
                                    meta.beta = c(meta$beta),
                                    meta.p = c(meta$pval),
                                    meta.ci.lb = c(meta$ci.lb),
                                    meta.ci.ub = c(meta$ci.ub),
                                    qep = meta$QEp))
  
}
meta_res <- meta_res %>%
  mutate_at(all_of(vars(contains("meta"))), ~as.numeric(as.character(.))) %>%
  mutate(meta.p.adj = p.adjust(meta.p, method = "bonferroni"))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("1", "2", "3")
col_labs <- c("Downregulating", "Upregulating", "NS (FDR p > 0.05)")

data_names <- c("snuc", "kellis", "multiome", "meta")

meta_row <- meta_res %>% 
  mutate(lab ="meta", 
         coef = ifelse(meta.beta > 0, -log10(meta.p.adj), log10(meta.p.adj))) %>%
  dplyr::select(lab, reg, coef)

res_row <- res %>% 
  filter(data %in% c("snuc", "kellis", "multiome")) %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
  mutate(lab = data, 
         coef = ifelse(coef > 0, -log10(p.adj), log10(p.adj))) %>%
  dplyr::select(lab, reg, coef)

mat <- rbind(res_row, meta_row) %>%
  mutate(lab = factor(lab, 
                      levels = c("snuc", "multiome", "kellis", "meta"), 
                      labels = c("CUIMC1 (n=424)", "CUIMC2 (n=212)", "MIT (n=130)", "Meta-analysis")),
         reg = factor(reg, levels = c("ARID5B", "CEBPA", "MITF", "PPARG", "IRF7"))) %>%
  spread(reg, coef) %>%
  column_to_rownames("lab") %>%
  as.matrix()

col_split <- factor(c(rep("GPNMB-high", 4), "IFN"))
row_split <- factor(c(rep("data", 3), "meta"))

p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(0, max(mat, na.rm=T)), 
                              colors = c("white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & abs(mat[i, j]) > -log10(0.05)) {
                 grid.points(pch = 8, x, y, size = unit(0.075, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1),
             cluster_rows = F,
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_names_side = "left",
             row_dend_side = "right",
             row_title_gp = gpar(fontsize = 0),
             row_gap = unit(0.1, "in"),
             row_split = row_split,
             column_split = col_split,
             column_title_side = "top",
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             top_annotation = HeatmapAnnotation(foo = anno_empty(border = FALSE,  height = unit(0.1, "in")), which = "column"),
             column_gap = unit(0.1, "in"),
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 10),
             column_names_rot = 45,
             cluster_columns = F,
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(Bonferroni p-value)'), 
               border = "black",
               labels_gp = gpar(font = 10),
               title_gp = gpar(font = 10),
               legend_width = unit(1.5, "in"),
               direction = "horizontal",
               title_position = "topcenter"
             ))

pdf(file.path(res_path, paste0(fig, "_panel_D.pdf")), height = 3, width = 4)
draw(p, 
     ht_gap = unit(0.2, "in"),
     heatmap_legend_side = "bottom", 
     align_heatmap_legend = "heatmap_center", 
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))

for(i in 1:2){
  decorate_annotation("foo", slice = i, { grid.rect(x = 0, height = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left") })
}

dev.off()
#===============================================================================
#
# PANEL E
#
#===============================================================================
aracne_sc <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
aracne_sn <- read_csv(file.path(home_path, "aracne_sn/aracne_w_tfmode.csv"))

targets_ARID5B <- aracne_sc %>% filter(Regulator == "ARID5B") %>% pull(Target)
targets_CEBPA <- aracne_sc %>% filter(Regulator == "CEBPA") %>% pull(Target)
targets_MITF <- aracne_sc %>% filter(Regulator == "MITF") %>% pull(Target)
targets_PPARG <- aracne_sc %>% filter(Regulator == "PPARG") %>% pull(Target)

p <- venn.diagram(list(targets_ARID5B, targets_CEBPA, targets_MITF, targets_PPARG),
                  category.names = c("ARID5B" , "CEBPA" , "MITF", "PPARG"),
                  filename = NULL,
                  lwd = 1.5,
                  col = proj_cols("blue", "yellow", "teal", "pink"),
                  fill = alpha(proj_cols("blue", "yellow", "teal", "pink"), 0.3),
                  cex = 0.75,
                  fontfamily = "Helvetica",
                  cat.pos = c(-10, 10, 0, 0),
                  cat.cex = 1,
                  cat.default.pos = "outer",
                  cat.fontfamily = "Helvetica",
                  cat.fontface = "bold",
                  cat.col = proj_cols("blue", "yellow", "teal", "pink"),
                  print.mode = c("raw", "percent"),
                  sigdigs = 2)

pdf(file.path(res_path, paste0(fig, "_panel_E.pdf")), height = 3.5, width = 3.5)
print(cowplot::plot_grid(p))
dev.off()
#===============================================================================
#
# PANEL F
#
#===============================================================================
regs <- c("ARID5B", "CEBPA", "MITF", "PPARG")

res <- data.frame()

for(j in c("snuc", "multiome", "kellis")) { 
  
  data <- eval(parse(text = j))
  
  counts <- data@meta.data %>% 
    group_by(projid) %>% 
    summarise(factor = mean(scHPF_26), 
              ARID5B = mean(ARID5B), 
              umi = mean(nCount_RNA), 
              MITF = mean(MITF), 
              CEBPA = mean(CEBPA), 
              PPARG = mean(PPARG),
              pmi = unique(pmi),
              age = unique(age_death), 
              sex = unique(msex),
              study = unique(study))%>%
    drop_na()
  
  other_vars <- subset(counts, select = c(umi, pmi, study, age, sex))
  tf_vars <- subset(counts, select = c(ARID5B, CEBPA, MITF, PPARG))
  
  part <- summary(varpart(Y = counts$factor, tf_vars, other_vars))
  sig_adj <- anova.cca(rda(Y = counts$factor, tf_vars, other_vars), permutations = 10000)
  res <- rbind(res, cbind(dataset = j, reg = "all.unique", part = part$uniqpart[1], sig =  sig_adj$`Pr(>F)`[1]))
  
  for(i in c("ARID5B", "CEBPA", "MITF", "PPARG")){
    other_vars <- counts %>% dplyr::select(-projid, -i, -factor)
    part <- varpart(Y = counts$factor, counts[,i], other_vars)
    sum_part <- summary(part)
    sig <- anova.cca(rda(Y = counts$factor, counts[,i], other_vars), permutations = 10000)
    res <- rbind(res, cbind(dataset = j, reg = i, part = sum_part$uniqpart[1], sig =  sig$`Pr(>F)`[1]))
  }
}

cols <- unname(c("#282828", proj_cols("blue", "yellow", "teal", "pink")))
names(cols) <- c("All 4 TFs", "ARID5B", "CEBPA", "MITF","PPARG")

values <- res %>% mutate(part = as.numeric(as.character(part)))

sig <- res %>%
  mutate(part = as.numeric(as.character(part))) %>%
  group_by(dataset) %>%
  mutate(sig = p.adjust(sig, method = "BH"),
         lab = ifelse(sig < 0.001, "***", ifelse(sig < 0.01, "**", ifelse(sig < 0.05, "*", "")))) %>%
  dplyr::select(dataset, reg, part, lab) %>% 
  ungroup()

plot_data <- left_join(values, sig) %>%
  mutate(part = as.numeric(as.character(part))) %>%
  mutate(lab = paste0(formatC(part*100, digits = 2, format = "g"), "%", lab)) %>%
  mutate(reg = factor(reg, 
                      levels = rev(c("all.unique", "ARID5B", "CEBPA", "MITF", "PPARG")), 
                      labels = rev(c("All 4 TFs", "ARID5B", "CEBPA", "MITF", "PPARG"))),
         dataset = factor(dataset, 
                          levels = c("snuc", "multiome", "kellis"),
                          labels = c("**CUIMC1** (*n*=424)", "**CUIMC2** (*n*=212)", "**MIT** (*n*=130)")))

p <- ggplot(plot_data, aes(y = reg, x = part, fill = reg)) + 
  geom_col() +
  geom_text(aes(label = lab,
                hjust = ifelse(part < 0.2, 0, 1),
                color = ifelse(part == "CEBPA" | part < 0.2, "1", "0")), 
            family = "Helvetica")+
  facet_wrap(~dataset, nrow = 1) + 
  scale_fill_manual(values = cols)+
  scale_color_manual(values = c("white", "black"))+
  labs(x = "Proportion of variance explained")+
  scale_x_continuous(limits = c(0, 0.6), breaks = c(0, 0.25, 0.5), labels = c(0, 0.25, 0.5))+
  theme_proj() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(size = 12),
        axis.text.x = element_markdown(size = 12),
        axis.title.x = element_markdown(size = 12),
        panel.grid.major = element_blank(), 
        strip.text = element_markdown(face = "plain", size = 10, margin = margin(t = 0.05, b = 0.05, unit = "in")))

pdf(file.path(res_path, paste0(fig, "_panel_F.pdf")), height = 3, width = 6)
print(p)
dev.off()