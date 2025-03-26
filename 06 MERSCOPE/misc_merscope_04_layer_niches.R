#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE: Define 3 layer niches.
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib_spat")

library(logger)
library(patchwork)
library(Seurat); options(Seurat.object.assay.version = "v5")
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))
cell_type_annotations <- read_rds(file.path(home_path, "merscope/merged/cell_type_annotations.rds")) %>%
  dplyr::select(cell_names, orig.ident, cell_type, subtype)
keep_vars <- ls()

#===============================================================================
# 
# DATA
#
#===============================================================================
ref <- read_rds(file.path(home_path, "merscope/merged/merged_clustered.rds"))
DefaultAssay(ref) <- "RNA"
for(i in run_metadata$batch_code){ ref[[i]] <- NULL }
ref[["pca"]] <- NULL
ref[["harmony"]] <- NULL
ref[["SCT"]] <- NULL
gc()

ref@meta.data <- ref@meta.data %>% left_join(cell_type_annotations, by = c("orig.ident","cell_names"))
rownames(ref@meta.data) <- ref@meta.data$cell_names

keep_cells <- rownames(ref@meta.data %>% filter(cell_type != "Ambiguous"))
ref <- subset(ref, subset = cell_names %in% keep_cells)

ref@meta.data$niche_wm <- rep(NA, nrow(ref@meta.data))

#===============================================================================
# 
# NICHES: GRAY MATTER VS. WHITE MATTER
#
#===============================================================================
coords <- data.frame()
for(i in 1:nrow(run_metadata)){
  
  res_path <- file.path(home_path, "merscope", run_metadata$run_name[i])
  
  # DATA
  #-----------------------------------------------------------------------------
  data <- read_rds(file.path(home_path, "merscope", run_metadata$run_name[i], "data_proc.rds"))
  data[["RNA"]] <- data[["Vizgen"]]
  DefaultAssay(data) <- "RNA"
  data[["Vizgen"]] <- NULL
  gc()
  
  # CELL TYPE ANNOTATIONS
  #-----------------------------------------------------------------------------
  data@meta.data <- data@meta.data %>% left_join(cell_type_annotations, by = c("fov.y"="orig.ident","cell_names"))
  rownames(data@meta.data) <- data@meta.data$cell_names
  
  # ARRANGE PHYSICAL SPACE
  #-------------------------------------------------------------------------------
  tmp_coords <- as.data.frame(data@images[[run_metadata$batch_code[i]]]@boundaries$centroids@coords)
  tmp_coords$batch <- run_metadata$batch_code[i]
  tmp_coords$donor <- run_metadata$donor[i]
  tmp_coords$cell_names <- data@images[[run_metadata$batch_code[i]]]@boundaries$centroids@cells
  
  if(i > 1){  tmp_coords[,1] = tmp_coords[,1] + max(coords[,1]) + 2000 }
  coords <- rbind(coords, tmp_coords)
  
  # SUBSET
  #-------------------------------------------------------------------------------
  keep_cells <- rownames(data@meta.data %>% filter(!is.na(cell_type) & cell_type != "Ambiguous"))
  data <- subset(data, subset = cell_names %in% keep_cells)
  data$cell_type <- factor(data$cell_type, levels = names(cell_type_cols))
  
  data <- UpdateSeuratObject(data)
  data <- NormalizeData(data) %>% ScaleData()
  
  # NICHES: WM vs GM
  #-----------------------------------------------------------------------------
  set.seed(1234)
  data <- BuildNicheAssay(data,
                         fov = run_metadata$batch_code[i],
                         group.by = "cell_type",
                         neighbors.k = 50,
                         niches.k = 2,
                         assay = "niches_wm",
                         cluster.name = "wm")
  
  Idents(data) <- "wm"
  
  de <- FindAllMarkers(data,
                       assay = "RNA",
                       slot = "data",
                       logfc.threshold = 0.5,
                       test.use = "MAST",
                       min.pct = 0.1,
                       only.pos = T,
                       max.cells.per.ident = 2500,
                       min.cells.group = 500,
                       random.seed = 1234,
                       latent.vars = "volume") %>%
    group_by(cluster) %>%
    top_n(n = -10, wt = p_val) %>%
    summarise(markers = paste0(gene, collapse = ","))
  
  gm_cluster <- de %>% filter(grepl("MOG", markers)==F)
  
  wm_labels <- ifelse(data@meta.data$wm == gm_cluster$cluster[1], "GM", "WM")
  names(wm_labels) <- rownames(data@meta.data)
  
  ref <- AddMetaData(ref, wm_labels, "niche_wm")
  
  rm(list = c("data", "keep_cells", "de", "wm_labels", "gm_cluster", "res_path", "tmp_coords")); gc()
}

#===============================================================================
# 
# NICHES: GRAY MATTER LAYERS
#
#===============================================================================
coords <- coords %>% column_to_rownames("cell_names")
ref@images[["merged"]] <- CreateFOV(coords[,1:2], type = "centroids", name = "merged")

keep_cells <- rownames(ref@meta.data %>% filter(niche_wm == "GM"))
ref_gm <- subset(ref, subset = cell_names %in% keep_cells)

set.seed(1234)
ref_gm <- BuildNicheAssay(ref_gm, 
                         fov = "merged", 
                         group.by = "cell_type",
                         neighbors.k = 75,
                         niches.k = 2,
                         assay = "niches_layers",
                         cluster.name = "layers")

Idents(ref_gm) <- "layers"
DefaultAssay(ref_gm) <- "downsampled"

set.seed(1234)
keep_cells <- ref_gm@meta.data %>% group_by(orig.ident, layers) %>% slice_sample(n = 2500) %>% pull(cell_names)
tmp <- subset(ref_gm, subset = cell_names %in% keep_cells)

de <- FindAllMarkers(tmp, 
                     logfc.threshold = 0.25,
                     test.use = "MAST",
                     min.pct = 0.1,
                     only.pos = T,
                     max.cells.per.ident = 2500,
                     min.cells.group = 500,
                     random.seed = 1234,
                     latent.vars = c("volume", "orig.ident"))

write_csv(de, file.path(home_path, "merscope/merged/niche_DE.csv"))

layer_labels <- ref_gm@meta.data %>% 
  dplyr::select(cell_names, layers) %>%
  left_join(de %>% dplyr::select(cluster, top_marker), by = c("layers"="cluster"))

tmp <- ref@meta.data %>%
  left_join(layer_labels) %>%
  mutate(niches = ifelse(niche_wm == "WM", "WM", NA),
         niches = ifelse(niche_wm == "GM" & !is.na(layers), paste0("GM.", top_marker), niches),
         niches = recode(niches, "GM.NA"="GM.Other"))

layer_niches <- tmp$niches
names(layer_niches) <- tmp$cell_names

ref <- AddMetaData(ref, layer_niches, "layer_niches")

saveRDS(ref@meta.data %>% dplyr::select(orig.ident, cell_names, layer_niches), 
        file.path(home_path, "merscope/merged/niche_assignments.rds"))

saveRDS(ref[["merged"]],
        file.path(home_path, "merscope/merged/merged_spatial_coords.rds"))
