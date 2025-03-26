#!/usr/bin/Rscript

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib")

library(Seurat)
library(SeuratObject)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/path/to"
res_path <- file.path(home_path, "merscope/merged")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_helpers.R"))

# DATA
#-------------------------------------------------------------------------------
run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))

data <- read_rds(file.path(res_path, "merged_clustered.rds"))
DefaultAssay(data) <- "RNA"
for(i in run_metadata$batch_code){
  data[[i]] <- NULL
}
data[["pca"]] <- NULL
data[["harmony"]] <- NULL
data[["SCT"]] <- NULL
gc()

# SELECT RESOLUTION
#-------------------------------------------------------------------------------
data <- rename_clusters(data, rename_from = "SCT_snn_res.0.8", rename_to = "final_clusters")
Idents(data) <- "final_clusters"

# EXCLUDE SMALL CLUSTERS
#-------------------------------------------------------------------------------
rm_clust <- names(which(table(data$final_clusters) <= 50))

keep_cells <- data@meta.data %>% 
  filter(!is.na(final_clusters) & !final_clusters %in% rm_clust) %>% 
  pull(cell_names)
data <- subset(data, subset = cell_names %in% keep_cells)

# ANNOTATE CELLS
#-------------------------------------------------------------------------------
cell_types <- list(`Astrocytes` = c(2,10,29,23,36,39),
                  `Endothelial` = c(7,26),
                  `Inhibitory`= c(9,17,22,24),
                  `Excitatory` = c(1,8,12,16,19,25,28,34,37),
                  `Microglia` = c(6),
                  `Monocytes/Macrophages`= c(35),
                  `Oligodendrocytes` = c(3:5,18,23,20,21,27,32,30),
                  `OPCs` = c(14,15),
                  `Pericytes` = c(11),
                  `Adaptive Immune` = c(31,33),
                  `Ambiguous`= c(13,23,38))

data@meta.data <- data@meta.data %>%
  mutate(cell_type = ifelse(final_clusters %in% cell_types[["Oligodendrocytes"]], "Oligodendrocytes", final_clusters),
         cell_type = ifelse(final_clusters %in% cell_types[["OPCs"]], "OPCs", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Astrocytes"]], "Astrocytes", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Microglia"]], "Microglia", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Monocytes/Macrophages"]], "Monocytes/Macrophages", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Endothelial"]], "Endothelial", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Pericytes"]], "Pericytes", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Excitatory"]], "Excitatory", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Inhibitory"]], "Inhibitory", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Adaptive Immune"]], "Adaptive Immune", cell_type),
         cell_type = ifelse(final_clusters %in% cell_types[["Ambiguous"]], "Ambiguous", cell_type),
         cell_type = factor(cell_type, levels = names(cell_types)))

subtype_labs <- data@meta.data %>%
  group_by(cell_type, final_clusters) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(abb = recode(cell_type, !!!cell_type_abb),
         subtype = paste0(abb, ".", row_number())) %>%
  arrange(abb, final_clusters) %>%
  mutate(subtype = factor(subtype, levels = subtype))%>% 
  dplyr::select(cell_type, final_clusters, subtype)

# SAVE
#-------------------------------------------------------------------------------
data@meta.data <- data@meta.data %>% left_join(subtype_labs , by = c("cell_type","final_clusters"))
rownames(data@meta.data) <- data@meta.data$cell_names

saveRDS(data@meta.data %>% dplyr::select(orig.ident, batch_name, cell_names, final_clusters, cell_type, subtype), 
        file.path(res_path, "cell_type_annotations.rds"))
