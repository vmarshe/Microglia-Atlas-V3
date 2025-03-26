#!/usr/bin/Rscript

#===============================================================================
#
# DISCOVERY ARACNE MODEL
# This script calculates the enrichment of regulons across scHPF factors.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sc")

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
# PANEL B - ENRICHMENT
#
#===============================================================================

# GENE MAPPING
#-------------------------------------------------------------------------------
universe <- rownames(schpf@feature.loadings)

terms <- data.frame()
for(i in 1:length(factors)){
  terms <- rbind(terms, 
                 data.frame(gsid = unname(factors[i]),
                            Gene = names(sort(schpf@feature.loadings[,factors[i]], decreasing = T)[1:100])))
}
# this is a temporary solution to: https://github.com/YuLab-SMU/clusterProfiler/issues/283 to avoid truncation of the universe
terms <- rbind(terms, cbind(gsid = "universe", Gene = universe))

# GENE LIST
#-------------------------------------------------------------------------------
gene_list <- list()
for(i in 1:length(regulators)){
  tmp <- data.frame(gene = aracne %>% filter(Regulator %in% regulators[i] & dir == "up") %>% pull(Target))
  gene_list[[i]] <- tmp$gene
  names(gene_list)[i] <- regulators[i]
}

# ENRICHER
#-------------------------------------------------------------------------------
res <- data.frame()
for(i in 1:length(gene_list)){
  
  genes <- gene_list[[i]]
  
  out <- enricher(genes, 
                  universe = universe, 
                  TERM2GENE = terms,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  minGSSize = 10,
                  maxGSSize = 2500)
  
  if(length(out) > 0){
    out@result$signature <- names(gene_list)[i]
    res <- rbind(res, out@result)
  } else {
    log_info("No enriched TFs: ", names(gene_list)[i])
  }
}

res$p.adj.bonf = p.adjust(res$pvalue, method = "bonferroni")
write_xlsx(res, file.path(res_path, "tab_aracne_scHPF_enrichment.xlsx"))