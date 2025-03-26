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
library(clusterProfiler)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(readxl)
library(rrvgo)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "metabolism")
fig <- "fig_supp03_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))

#===============================================================================
#
# PANEL A
#
#===============================================================================

# MAP GENES
#-------------------------------------------------------------------------------
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))
entrez_names <- read_rds(file.path(home_path, "data/entrez_mapping.rds"))

# GO TREE PLOTS
#-------------------------------------------------------------------------------
for(i in paste0("scHPF_", c(2, 8, 18, 21))){
  
  genes <- names(sort(schpf@feature.loadings[,i], decreasing = T))[1:100]
  
  genes <- entrez_names %>% filter(gene %in% genes) 
  n_missing_genes <- sum(is.na(genes$entrez))
  n_non_missing_genes <- sum(!is.na(genes$entrez))
  
  set.seed(1234)
  genes <- genes %>% filter(!is.na(entrez)) %>% group_by(gene) %>% sample_n(1) %>% pull(entrez)
  
  out <- enrichGO(gene = genes,
                  universe = as.character(na.omit(entrez_names$entrez)),
                  OrgDb  = org.Hs.eg.db,
                  ont  = "BP",
                  pAdjustMethod = "BH",
                  readable = TRUE)@result %>%
    as.data.frame() %>%
    filter(p.adjust < 0.05 & qvalue < 0.05)
  
  if(length(out$ID) > 0) {
    
    simMatrix <- calculateSimMatrix(out$ID, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
    scores <- setNames(-log10(out$qvalue), out$ID)
    reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
    
    pdf(file.path(res_path, paste0(fig, "panel_A_", i, ".pdf")), height = 3, width = 3)
    treemapPlot(reducedTerms)
    dev.off()
  }
}

#===============================================================================
#
# PANEL B
#
#===============================================================================
gene_lists <- list()
for(i in factors){
  tmp <- list(data.frame(gene = names(sort(schpf@feature.loadings[,i], decreasing = T)[1:100])) %>% pull(gene))
  gene_lists <- c(gene_lists, tmp)
}
names(gene_lists) <- names(factors)

mitocarta <- read_xls(file.path(home_path, "metabolism/Human.MitoCarta3.0.xls"), sheet = 4) %>% drop_na(MitoPathway)
all_genes <- read_delim(file.path(home_path, "data/scRNAseq_genes.txt"), delim = "\t")$genes

universe <- read_xls(file.path(res_path, "Human.MitoCarta3.0.xls"), sheet = 3) %>% pull(Symbol)
universe <- intersect(universe, all_genes)

terms <- data.frame()
for(i in 1:nrow(mitocarta)){
  terms <- rbind(terms, 
                data.frame(gsid = paste0(mitocarta$`MitoPathway`[i], " (L", length(unlist(str_extract_all(mitocarta$`MitoPathways Hierarchy`[i], "[>]"))), ")"),
                           Gene = gsub(" ", "", unlist(str_split(mitocarta$Genes[i], "[,]")))))
}
# this is a temporary solution to: https://github.com/YuLab-SMU/clusterProfiler/issues/283
# to avoid truncation of the universe
terms <- rbind(terms, cbind(gsid = "universe", Gene = universe))

# RUN ENRICHMENT
#-------------------------------------------------------------------------------
res <- data.frame()
for(i in 1:length(gene_lists)){
  
  genes <- gene_lists[[i]]
  
  out = enricher(genes, 
                 universe = universe, 
                 TERM2GENE = terms,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 minGSSize = 10,
                 maxGSSize = 500) %>%
    mutate(signature = names(gene_lists)[i]) 
  
  res <- rbind(res, out@result)
}
rownames(res) <- NULL

res %>% write_xlsx(file.path(res_path, "tab_MitoCarta_enrichment.xlsx"))

# keep pathways with at least one significant hit
keep <- res %>% filter(p.adjust < 0.05 & qvalue < 0.05) %>% pull(ID)
# keep factors with at least one significant hit
keep_factors <- res %>% 
  filter(p.adjust < 0.05 & qvalue < 0.05) %>% 
  group_by(signature) %>% 
  sample_n(1) %>% 
  pull(signature)
keep_factors <- c("Glycolysis (2)", keep_factors)

mat <- res %>% 
  filter(ID %in% keep & signature %in% keep_factors) %>%
  rowwise() %>%
  dplyr::select(ID, qvalue, signature) %>%
  mutate(qvalue = -log10(qvalue)) %>%
  mutate(signature = factor(signature, levels = names(factors))) %>%
  mutate(ID = factor(ID, levels = unique(terms$gsid))) %>%
  arrange(signature) %>%
  spread(ID, qvalue) %>%
  column_to_rownames("signature") %>%
  mutate_all(~as.numeric(as.character(.))) %>%
  as.matrix()

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 4.5, width = 12)
set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(0, max(mat, na.rm=T)), 
                              colors = c("white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & (mat[i, j] >-log10(0.05) | mat[i, j] <log10(0.05))) {
                 grid.points(pch = 8, x, y, size = unit(0.075, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             cluster_rows = F,
             row_names_side = "left",
             row_dend_side = "right",
             row_names_max_width = unit(6, "in"),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             cluster_columns = F,
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'\n(FDR q-value)'), 
               legend_width = unit(1.5, "in"),
               direction = "vertical",
               title_position = "topleft",
               border = "black"
             ))

draw(p, heatmap_legend_side = "right", background = "transparent", padding = unit(c(0.2, 1, 0.2, 0.2), "in"))
decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2)) })
dev.off()

