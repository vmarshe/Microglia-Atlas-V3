#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(clusterProfiler)
library(ggtext)
library(org.Hs.eg.db)
library(patchwork)
library(rrvgo)
library(tidyverse)
library(viridis)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "DAM")
fig <- "fig_supp02_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))

#===============================================================================
#
# PANEL A
#
#===============================================================================

# GENES
#-------------------------------------------------------------------------------
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))
entrez_names <- read_rds(file.path(home_path, "data/entrez_mapping.rds"))

# GO TREE PLOTS
#-------------------------------------------------------------------------------
for(i in paste0("scHPF_", c(11,26))){
  
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
    
    pdf(file.path(res_path, paste0(fig, i, ".pdf")), height = 3, width = 3)
    treemapPlot(reducedTerms)
    dev.off()
  }
}

#===============================================================================
#
# PANEL B
#
#===============================================================================
dam <- sort(c("HLA-DRB1", "MAFB", "CHI3L1", "LGALS3", "MS4A7", "CD9", 
              "TREM2", "LPL", "ITGAX", "APOE", "TYROBP", "CTS7",
              "HEXA", "CLEC7A", "CSF1"))

homeostatic <- sort(c("LYVE1", "ABCG2", "OLFML3", "CST3", "ITM2C", "P2RY13", 
                      "CX3CR1", "SELPLG", "TMEM119", "CSF1R", 
                      "SPARC", "SERINC3", "RGS2","CCR5"))

features <- c(dam, homeostatic)[c(dam, homeostatic) %in% rownames(schpf@feature.loadings)]

new_labs <- gsub('-high', paste0("<sup>high</sup>"), rev(names(factors)))

tmp <- schpf@feature.loadings %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  gather("factor", "loading", scHPF_1:scHPF_26) %>%
  group_by(factor) %>%
  arrange(factor, desc(loading)) %>%
  mutate(rank = row_number(), 
         rank_pct = 1-(row_number()/n())) %>%
  filter(genes %in% features) %>%
  mutate(genes = factor(genes, levels = features), 
         factor = factor(factor, levels = rev(factors), labels = new_labs),
         top = ifelse(rank_pct >= 1-(100/nrow(schpf@feature.loadings)), "*", ""))


p <- ggplot(tmp, aes(x = factor, y = genes, fill = rank_pct, size = rank_pct))+
  geom_point(shape = 21, stroke = 0.5) +
  geom_text(aes(label = top), size = 8, show.legend = F, hjust = 0.5, vjust = 0.75) +
  scale_fill_viridis(name = "rank", option = "magma", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_size_continuous(name = "rank", limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1))+
  theme_proj() +
  coord_flip()+
  guides(size = guide_legend(title = "% Rank", 
                             order = 2, 
                             byrow = T, 
                             theme = theme(legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"), 
                                           legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in"), family = "Helvetica"),
                                           legend.margin = margin(l = -0.1, unit = "in"))), 
         fill = guide_colorbar(title = "% Rank", 
                               order = 1, 
                               barheight = unit(1, "in"),  
                               barwidth = unit(0.2, "in"),  
                               theme = theme(legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
                                             legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in"), family = "Helvetica"),
                                             legend.margin = margin(l = -0.1, unit = "in"))))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 12, angle = 45, hjust = 1, family = "Helvetica", face = "italic"),
        axis.text.y = element_markdown(size = 12, family = "Helvetica"),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75))

cairo_pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 6, width = 9.5)
print(p)
dev.off()
