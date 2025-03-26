#!/usr/bin/env Rscript

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib")

library(ggtext)
library(logger)
library(optparse)
library(patchwork)
library(progressr)
library(Seurat)
library(SeuratObject)
library(tidyverse)

opts <- list(make_option(c("--dir"), action = "store", type = "character"))
opts <- parse_args(OptionParser(option_list=opts))

# PATHS
#-------------------------------------------------------------------------------
data_path <- "~/data/MERSCOPE"
out_path <- "~/projects/microglia_atlas/merscope"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_merscope.R"))

run <- unlist(str_split(opts$dir, "/"))[[9]]
run_metadata <- read_csv(file.path(out_path, "batch_names.csv"))
run_metadata <- run_metadata %>% filter(run_name %in% run)

segmentation <- "cp_2D_allLayers_nuclei"
#===============================================================================
#
# QC & DATA PREP
#
#===============================================================================
print(paste0("Processing: ", run_metadata$batch_code))
dir_path <- file.path(data_path, run_metadata$run_name, segmentation)
res_path <- file.path(out_path, run_metadata$run_name)
fov <- run_metadata$batch_code

if(dir.exists(dir_path)){
  if(file.exists(file.path(dir_path, "detected_transcripts.csv"))){
    
    
    if(file.exists(file.path(res_path, "get_data.log"))) {file.remove(file.path(res_path, "get_data.log"))}
    log_threshold(TRACE)
    log_appender(appender_tee(file = file.path(res_path, "get_data.log")))
    
    data <- load_merscope(folder_path = dir_path, fov_name = fov)
    saveRDS(data, file = file.path(res_path, "data_preqc.rds"))
    
    # MAIN QC
    #-------------------------------------------------------------------------
    data <- qc_spatial(data = data,
                       min_umi = max(quantile(data$nCount_Vizgen, 0.025), 25),
                       max_umi = min(quantile(data$nCount_Vizgen, 0.975), 2500),
                       min_genes = max(quantile(data$nFeature_Vizgen, 0.025), 5),
                       min_volume = max(quantile(data$volume, 0.025), 50),
                       max_volume = min(quantile(data$volume, 0.975), 2500),
                       fov_name = fov,
                       assay = "Vizgen",
                       size_feature = "volume",
                       res_path = res_path)
    
    # FORM FACTOR
    #-------------------------------------------------------------------------
    geometry <- get_cell_geometry(data, fov_names = fov)
    
    data@meta.data <- data@meta.data %>% left_join(geometry, by = "cell_names")
    rownames(data@meta.data) <- data@meta.data$cell_names
    
    data@meta.data$solidity <- data@meta.data$solidity
    data@meta.data$form_factor <- data@meta.data$form_factor
    
    p <- ImageFeaturePlot(data, features = "form_factor") +
      scale_fill_gradient(name = "Form\nfactor", low = "lightgrey", high = proj_cols("blue"))+
      theme_umap()+
      theme(plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
            plot.title = element_blank(),
            legend.key.width = unit(0.3, "in"),
            legend.key.height = unit(0.3, "in"),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 12))
    
    pdf(file.path(res_path, "form_factor.pdf"),  height = 4, width = 4)
    print(p)
    dev.off()
    
    lab <- paste0("<span style='color:", proj_cols("red"), "'>Median:</span> ", round(median(data@meta.data$form_factor), 2), "<br>",
                  "<span style='color:", proj_cols("yellow"), "'>Mean:</span> ", round(mean(data@meta.data$form_factor), 2),  "<br>",
                  "SD: ", round(sd(data@meta.data$form_factor), 2))
    
    p <- ggplot(data@meta.data, aes(x = form_factor)) +
      geom_histogram(color = proj_cols("blue"), fill =  proj_cols("blue"), alpha = 0.8) +
      geom_vline(xintercept = median(data@meta.data$form_factor), lty = 2, color = proj_cols("red"))+
      geom_vline(xintercept = mean(data@meta.data$form_factor), lty = 2, color = proj_cols("yellow"))+
      annotate(geom='richtext', x = Inf, y =Inf, hjust= 1, vjust = 1, label=lab) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "Form factor (higher is worse)", y = "Number of cells")+
      theme_proj() +
      theme(axis.title = element_text(size = 12))
    
    pdf(file.path(res_path, "form_factor_hist.pdf"),  height = 3, width = 5)
    print(p)
    dev.off()
    
    # QC on FORM FACTOR
    #-------------------------------------------------------------------------
    form_factor_thresh <- quantile(data@meta.data$form_factor, 0.02)
    keep_cells <- data@meta.data$cell_names[which(data@meta.data$form_factor > form_factor_thresh)]
    data <- subset(data, subset = cell_names %in% keep_cells)
    
    # AMBIENCE
    #-------------------------------------------------------------------------
    detected_transcripts <- read_csv(file.path(dir_path, "detected_transcripts.csv"))
    
    tot_mols <- nrow(detected_transcripts)
    prop_assigned_mols <- nrow(detected_transcripts %>% filter(cell_id > 0))/tot_mols
    
    tot_cells <- ncol(data)
    tot_volume_occupied_cells <- sum(data$volume)
    density_mols <- nrow(detected_transcripts %>% filter(cell_id > 0))/tot_volume_occupied_cells*100
    
    tissue_metrics <- c(tot_mols, prop_assigned_mols, tot_cells, tot_volume_occupied_cells, density_mols)
    names(tissue_metrics) <- c("tot_mols", "prop_assigned_mols","tot_cells", "tot_vol_occupied",  "mol_density")
    
    write_csv(as.data.frame(tissue_metrics), file.path(res_path, "tissue_metrics.csv"), col_names=F)
    
    DefaultFOV(data) <- fov
    saveRDS(data, file = file.path(res_path, "data_proc.rds"))
    invisible(lapply(capture.output(sessionInfo()), log_info))
    
  }
}