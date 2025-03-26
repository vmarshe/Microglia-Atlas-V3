#!/usr/bin/env Rscript

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib2")

library(harmony)
library(logger)
library(Matrix)
library(matrixStats)
library(patchwork)
library(progressr)
library(sctransform)
library(Seurat)
library(SeuratObject)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
data_path <- file.path(home_path, "merscope")
res_path <- file.path(data_path, "merscope/merged")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/funcs_integration.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_helpers.R"))

# GET DATA
#-------------------------------------------------------------------------------
run_metadata <- read_csv(file.path(data_path, "batch_names.csv"))
cols <- run_metadata$cols
names(cols) <- run_metadata$batch_name
run_cols <- cols

for(i in 1:nrow(run_metadata)){
  if(file.exists(file.path(data_path, run_metadata$run_name[i], "data_proc.rds"))){
    log_info("Reading: ", run_metadata$batch_code[i])
    tmp <- read_rds(file.path(data_path, run_metadata$run_name[i], "data_proc.rds"))
    tmp$orig.ident <- run_metadata$batch_code[i]
    tmp[["RNA"]] <- tmp[["Vizgen"]]
    DefaultAssay(tmp) <- "RNA"
    tmp[["Vizgen"]] <- NULL
    
    assign(run_metadata$batch_code[i], tmp)
    rm(tmp); gc()
  }
}

# DOWNSAMPLE
#-------------------------------------------------------------------------------
downsample_prop <- do.call(rbind.data.frame, lapply(run_metadata$batch_code, function(x) eval(parse(text = x))@meta.data)) %>%
  group_by(orig.ident) %>%
  summarise(median = median(nCount_Vizgen)) %>%
  mutate(down_prop = min(median)/median)

for(i in 1:nrow(run_metadata)){

  if(any(grepl(run_metadata$batch_code[i], ls()))){
    tmp <- eval(parse(text = run_metadata$batch_code[i]))
    counts <- as.matrix(GetAssayData(subset(tmp, subset = orig.ident %in% run_metadata$batch_code[i]), assay = "RNA", slot = "counts"))
    prop <- downsample_prop %>% filter(orig.ident==run_metadata$batch_code[i]) %>% pull(down_prop)
    if(prop != 1){
      set.seed(1234)
      downsampled_counts <- DropletUtils::downsampleMatrix(counts, prop = prop)
      tmp[["downsampled"]] <- CreateAssayObject(downsampled_counts)
      assign(run_metadata$batch_code[i], tmp)

    } else {
      tmp[["downsampled"]] <- CreateAssayObject(counts)
      assign(run_metadata$batch_code[i], tmp)
    }
  }
}
write_csv(downsample_prop, file.path(res_path, "downsample_prop.csv"))

rm(list = c("downsample_prop", "tmp", "counts", "prop"))

# MERGE
#-------------------------------------------------------------------------------
merged <- merge(eval(parse(text = run_metadata$batch_code[1])), eval(parse(text = run_metadata$batch_code[2])))
for(i in run_metadata$batch_code[3:length(run_metadata$batch_code)]){
  log_info("Merging: ", i)
  merged <- merge(merged, eval(parse(text = i)))
  rm(list = c(i))
}
gc()

merged@meta.data$cell_names <- rownames(merged@meta.data)
merged@meta.data <- merged@meta.data %>% left_join(run_metadata, by = c("orig.ident" = "batch_code"))
rownames(merged@meta.data) <- merged@meta.data$cell_names

#===============================================================================
#
# NORMALIZATION
#
#===============================================================================
log_info("Normalizing...")
DefaultAssay(merged) <- "RNA"

merged <- NormalizeData(merged, assay = "RNA")
merged <- ScaleData(merged, assay = "RNA", features = rownames(merged))

merged <- SCTransform(merged,
                      assay = "RNA",
                      method = "glmGamPoi",
                      seed.use = 1234,
                      variable.features.n = nrow(merged),
                      vars.to.regress = "volume")

DefaultAssay(merged) <- "SCT"

#===============================================================================
#
# INTEGRATION
#
#===============================================================================

# PCA
#-------------------------------------------------------------------------------
merged <- RunPCA(merged, assay = "SCT", npcs = 40, seed.use = 1234, features = rownames(data))

# INTEGRATION w/ HARMONY
#-------------------------------------------------------------------------------
set.seed(1234)
pdf(file.path(res_path, "harmony_convergence.pdf"), height = 3, width = 5)
merged <- RunHarmony(merged, 
                     group.by.vars = "orig.ident",
                     dims.use = 1:40,
                     reduction.use = "pca",
                     verbose = T,
                     plot_convergence = T)
dev.off()

# UMAP w/o Integration
#-------------------------------------------------------------------------------
set.seed(1234)
tmp <- FindNeighbors(object = merged,
                     assay = "SCT",
                     reduction = "pca",
                     dims = 1:40,
                     k.param = 30)
set.seed(1234)
tmp <- RunUMAP(object = tmp,
               assay = "SCT",
               reduction = "pca",
               dims = 1:40,
               umap.method = "uwot",
               reduction.key = "UMAP_",
               seed.use = 1234,
               n.neighbors = 30,
               return.model = TRUE,
               reduction.name = "umap")

tmp <- calculate_lisi(data = tmp, split_by = "batch_name", reduction = "umap")

saveRDS(tmp, file.path(res_path, "unintegrated.rds"))

# UMAP w/ Harmony
#-------------------------------------------------------------------------------
set.seed(1234)
merged <- FindNeighbors(object = merged,
                        assay = "SCT",
                        reduction = "harmony",
                        dims = 1:40,
                        k.param = 30)
set.seed(1234)
merged <- RunUMAP(object = merged,
                  assay = "SCT",
                  reduction = "harmony",
                  dims = 1:40,
                  umap.method = "uwot",
                  reduction.key = "UMAP_",
                  seed.use = 1234,
                  n.neighbors = 30,
                  return.model = TRUE,
                  reduction.name = "umap")

merged <- calculate_lisi(data = merged, split_by = "batch_name", reduction = "umap")

#===============================================================================
#
# CLUSTERING
#
#===============================================================================
res <- seq(0.3, 1, 0.1)

merged <- FindClusters(object = merged,
                       graph.name = "SCT_snn",
                       resolution = res,
                       algorithm = 1,
                       random.seed = 1234,
                       method = "igraph",
                       group.singletons = F)

res_plots <- list()
for(i in 1:length(res)){
  
  set.seed(1234)
  cols <- pal()(length(unique(merged@meta.data %>% pull(paste0("SCT_snn_res.", res[i])))))
  cols <- sample(cols)
  
  order <- merged@meta.data %>%
    group_by_at(paste0("SCT_snn_res.", res[i])) %>%
    summarise(n = n()) %>%
    arrange(n) %>%
    mutate(cols = cols)
  
  res_plots[[i]] <- DimPlot(merged,  
                            reduction = "umap",  
                            group.by = paste0("SCT_snn_res.", res[i]), 
                            order = T, 
                            label = T,
                            cols = order$cols) +
    labs(title = paste0("resolution: ", res[i]))+
    theme_umap()+
    theme(legend.position = "none",
          legend.key = element_rect(color = NA),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(size = 18, face = "plain"),
          plot.margin = margin(t = 0.3, r = 0.3, b = 0.3, l = 0.1, unit = "cm"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
          axis.line.x = element_line(linewidth = 0.2, color = "black"),
          axis.line.y = element_line(linewidth = 0.2, color = "black"))
}

p <- wrap_plots(res_plots) + plot_layout(nrow = 2)
pdf(file.path(res_path, "resolution_comparison.pdf"), height = 9, width = 16)
print(p)
dev.off()

merged@meta.data$seurat_clusters <- NULL
saveRDS(merged, file.path(res_path, "merged_clustered.rds"))

# SELECT RESOLUTION
#-------------------------------------------------------------------------------
merged <- rename_clusters(merged, rename_from = "SCT_snn_res.0.8", rename_to = "final_clusters")
Idents(merged) <- "final_clusters"

#===============================================================================
#
# DEGS
#
#===============================================================================
Idents(merged) <- "final_clusters"

de <- FindAllMarkers(merged, 
                     assay = "downsampled",
                     slot = "data",
                     logfc.threshold = 0.25,
                     test.use = "MAST",
                     min.pct = 0.1,
                     only.pos = T,
                     max.cells.per.ident = 2000,
                     min.cells.group = 50,
                     random.seed = 1234,
                     latent.vars = "orig.ident")

de_markers <- de %>% 
  filter(p_val_adj < 0.001 & avg_log2FC > 0.5) %>%
  group_by(cluster) %>% 
  top_n(n = -20, wt = p_val_adj) %>% 
  summarise(genes = paste0(gene, collapse = ","),
            top_marker = gene[1])

saveRDS(de_markers, file.path(res_path, "major_cluster_markers.rds"))
saveRDS(de, file.path(res_path, "major_cluster_de.rds"))
