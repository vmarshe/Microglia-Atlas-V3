#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(SeuratObject)
library(Seurat)
library(tidyverse)
library(clusterProfiler); options(enrichment_force_universe = FALSE)
library(org.Hs.eg.db)
library(writexl)

source("~/projects/microglia_atlas/src/00 Helpers/factors.R")

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "GO")

#===============================================================================
#
# GO ENRICHMENT
#
#===============================================================================

# MAP GENES
#-------------------------------------------------------------------------------
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))
entrez_names <- read_rds(file.path(home_path, "data/entrez_mapping.rds"))

# ENRICHMENT
#-------------------------------------------------------------------------------
res <- data.frame()
for(i in factors){
  
  genes <- names(sort(schpf@feature.loadings[,i], decreasing = T))[1:100]
  genes <- entrez_names %>% filter(gene %in% genes) %>% pull(entrez)
  
  n_non_missing_genes <- length(genes)
  
  out <- enrichGO(gene = genes,
                  universe = as.character(na.omit(entrez_names$entrez)),
                  OrgDb  = org.Hs.eg.db,
                  ont  = "BP",
                  pAdjustMethod = "BH",
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = TRUE)@result %>%
    as.data.frame()
  
  res <- rbind(res, cbind(factor = i, out))
}
rownames(res) <- NULL

res <- res %>% filter(qvalue < 0.05 & p.adjust < 0.05)
write_xlsx(res, file.path(res_path, "tab_GO_enrichment.xlsx"))
