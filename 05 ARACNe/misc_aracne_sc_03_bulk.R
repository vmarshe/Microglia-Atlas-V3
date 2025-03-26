#!/usr/bin/Rscript

#===============================================================================
#
# DISCOVERY ARACNE MODEL
# This script calculates the enrichment of AD-related genes across ARACNe
# regulons.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(fgsea)
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
# PANEL E - BULK RNASEQ
#
#===============================================================================
data_path <- "~/data/Bulk RNAseq association [Annie Lee]"
files <- grep("dge_", list.files(data_path), value = T)
files <- grep("cogng_random_slope|ad_reagan|amyloid_sqrt|tangles", files, value = T)

# GENE MAPPING
#-------------------------------------------------------------------------------
ensembl_names <- read_rds(file.path(home_path, "data/ensembl_mapping.rds")) %>% drop_na()

sig_regs <- c("ARID5B", "MITF", "PPARG", "CEBPA", "IRF7")

fgsea_lists <- list()
for(i in 1:length(sig_regs)){
  tmp <- data.frame(gene = aracne %>% filter(Regulator %in% sig_regs[i] & tfmode > 0) %>% pull(Target))
  fgsea_lists[[i]] <- na.omit(tmp %>% left_join(ensembl_names, by = "gene") %>% pull(ensembl))
  names(fgsea_lists)[i] <- sig_regs[i]
}

#===============================================================================
#
# FGSEA
#
#===============================================================================
res <- data.frame()
for(i in 1:length(files)){
  
  file_name <- gsub("dge_|\\.txt", "", files[i])
  
  tmp <- read.table(file.path(data_path, files[i]), row.names = 1) %>% 
    rownames_to_column("gene") %>%
    filter(gene %in% ensembl_names$ensembl) %>%
    filter(padj < 0.05) %>% 
    arrange(desc(log2FoldChange))
  
  genes <- tmp$log2FoldChange
  names(genes) <- tmp$gene
  genes <- sort(genes, decreasing = T)
  
  set.seed(1234)
  fgsea_res <- fgsea(pathways = fgsea_lists, 
                     stats = genes,
                     minSize = 10,
                     scoreType = "std",
                     nPermSimple = 10000) %>%
    mutate(signature = file_name)
  
  res <- rbind(res, fgsea_res)
}
res$p.adj.bonf = p.adjust(res$pval, method = "bonferroni")
write_xlsx(res, file.path(res_path, "tab_aracne_bulkRNA_fgsea.xlsx"))
