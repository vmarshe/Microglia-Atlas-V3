#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(ComplexHeatmap)
library(circlize)
library(logger)
library(miloR)
library(patchwork)
library(scCustomize)
library(schex)
library(scran)
library(scater)
library(SeuratObject)
library(Seurat)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "proj_snuc/analysis")
fig <- "fig_s18"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_miloR.R"))

#===============================================================================
#
# DATA
#
#===============================================================================

# ROSMAP METADATA
#-------------------------------------------------------------------------------
rosmap_meta <- read_delim(file.path(home_path, "data/dataset_707_basic_11-01-2023.txt"), delim="\t")

vars <- grep("tangles_|nft_|amyloid_|plaq_n_|plaq_d_", names(rosmap_meta), v = T)
vars <- vars[vars != "amyloid_mes_t"]

rosmap_meta <- rosmap_meta %>%
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  arrange(projid) %>%
  dplyr::select(projid, msex, age_death, educ, pmi, study, amyloid, tangles, cogng_demog_slope, niareagansc, cogdx, all_of(vars)) %>%
  distinct() %>%
  mutate(pathoAD = factor(ifelse(niareagansc %in% c(3, 4), 0, 1), levels = c(0, 1)), 
         amyloid_sqrt = sqrt(amyloid), 
         tangles_sqrt = sqrt(tangles)) %>%
  mutate_at(all_of(vars), .funs = function(x) sqrt(x))

# SNUC
#-------------------------------------------------------------------------------
snuc <- readRDS(file.path(home_path, "proj_snuc/projection/snuc.rds"))
snuc@meta.data$cell_names <- rownames(snuc@meta.data)
snuc@meta.data <- snuc@meta.data %>% 
  left_join(snuc@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(snuc@meta.data) <- snuc@meta.data$cell_names

snuc@meta.data <- snuc@meta.data %>% mutate(projid = as.character(str_pad(projid, 8, pad = "0")))

for(i in c("projid", "study", "msex", "pmi",  "age_death", "educ", "pathoAD", "amyloid", "tangles",  "amyloid_sqrt", "tangles_sqrt")){
  
  tmp <- snuc@meta.data %>% dplyr::select(projid, cell_names) %>%
    left_join(rosmap_meta %>% dplyr::select(projid, all_of(i)))
  var <- tmp[,i]
  names(var) <- tmp$cell_names
  
  if(any(grepl(i, names(snuc@meta.data)))){ snuc@meta.data[,i] = NULL }
  
  snuc <- AddMetaData(snuc, var, i)
}

cell_meta_data <- readRDS(file.path(home_path, "proj_snuc/input/full_cell_metadata.rds")) %>%
  filter(cell.type == "Microglia") %>%
  dplyr::select(projid, annotation)
cell_meta_data$cell_names <- rownames(cell_meta_data)
snuc@meta.data <- snuc@meta.data %>% left_join(cell_meta_data)
rownames(snuc@meta.data) <- snuc@meta.data$cell_names

snuc_agg <- readRDS(file.path(home_path, "proj_snuc/projection/snuc_aggregated.rds"))%>%
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  left_join(rosmap_meta)

#===============================================================================
# 
# PANELS A-C, E-G
# 
#===============================================================================
for(trait in c("amyloid", "tangles")){
  
  k <- 30
  deltaFC <- 2
  
  # PATHS
  #-------------------------------------------------------------------------------
  res_path <- file.path(home_path, "proj_snuc", "analysis", paste0(trait, "_k", k, "_deltaFC", deltaFC))
  if(!dir.exists(res_path)) dir.create(res_path)
  log_appender(appender_tee(file.path(res_path, "log.txt")))
  
  # DATA
  #-----------------------------------------------------------------------------
  data <- DietSeurat(snuc, assays = "RNA", dimreducs = c("scHPF", "schpf.umap"))
  Idents(data) <- "cell_names"
  
  # CREATE OBJECT
  #-----------------------------------------------------------------------------
  log_info("Subsetting cells...")
  set.seed(1234)
  keep_cells <- data@meta.data %>%
    dplyr::select(projid, pathoAD, amyloid_sqrt, tangles_sqrt, pmi, msex, age_death, cell_names) %>%
    drop_na() %>%
    group_by(projid) %>%
    pull(cell_names)
  
  data@meta.data <- data@meta.data %>%
    mutate(amyloid = ifelse(amyloid_sqrt > median(amyloid_sqrt, na.rm=T), 1, 0)) %>%
    mutate(tangles = ifelse(tangles_sqrt > median(tangles_sqrt, na.rm=T), 1, 0))
  
  data <- subset(data, subset = cell_names %in% keep_cells)
  rownames(data@meta.data) <- data@meta.data$cell_names
  data@meta.data <- data@meta.data %>% rename_at(all_of(vars(trait)), ~"trait")
  data@meta.data$new_id <- paste0(data@meta.data$batch, "_", data@meta.data$projid)
  invisible(lapply(capture.output(data), log_info))
  
  # BUILD GRAPH
  #-----------------------------------------------------------------------------
  log_info("Creating Milo object...")
  milo <- Milo(as.SingleCellExperiment(data))
  invisible(lapply(capture.output(milo), log_info))
  
  set.seed(1234)
  milo <- buildGraph(milo, k = k, d = 23, reduced.dim = "SCHPF")
  
  set.seed(1234)
  milo <- makeNhoods(milo, prop = 0.1, k = k, d = 23,
                     refined = TRUE, reduced_dim = "SCHPF", refinement_scheme = "graph")
  
  plot_milo_nhood_size(milo, bins = 50, res_path = res_path)
  
  # GET NEIGHBORHOODS
  #-----------------------------------------------------------------------------
  log_info("Calculating neighborhoods...")
  set.seed(1234)
  milo <- countCells(milo,
                     meta.data = data@meta.data[,c("new_id", "trait", "study","msex", "age_death", "pmi", "educ")],
                     samples = "new_id")
  set.seed(1234)
  milo <- calcNhoodDistance(milo, d = 23, reduced.dim = "SCHPF")
  
  # EXPERIMENTAL DESIGN
  #-----------------------------------------------------------------------------
  log_info("Setting up experimental design...")
  design <- data@meta.data[,c("new_id", "trait", "study", "msex", "age_death", "pmi", "educ")]
  design$trait <- as.factor(design$trait)
  rownames(design) <- NULL
  design <- distinct(design)
  rownames(design) <- design$new_id
  
  # RUN DA
  #-----------------------------------------------------------------------------
  log_info("Performing DA testing...")
  formula <- as.formula("~ study + age_death + msex + pmi + trait")
  
  log_info("Formula: ", paste0(as.character(formula), collapse = ""))
  res <- testNhoods(milo, 
                    design = formula, 
                    design.df = design,
                    reduced.dim = "SCHPF", 
                    fdr.weighting = "graph-overlap")
  saveRDS(res, file.path(res_path, "testNhoods.rds"))
  
  log_info("Calculating neighborhood adjacency...")
  set.seed(1234)
  milo <-  buildNhoodGraph(milo, overlap = 3)
  saveRDS(milo@nhoods, file.path(res_path, "milo_Nhood_graph.rds"))
  
  log_info("Grouping neighborhoods...")
  da_results <- groupNhoods(x = milo, 
                            da.res = res, 
                            da.fdr = 0.05,
                            max.lfc.delta = deltaFC, 
                            overlap = 3, 
                            compute.new = F)
  saveRDS(da_results, file.path(res_path, "milo_grouped_Nhoods.rds"))
  
  # GET NEIGHBORHOOD COMPOSITION
  #-----------------------------------------------------------------------------
  log_info("Getting neighborhood composition...")
  get_nhood_comp(milo, res_path)
  
  if(!is.null(da_results)){
    
    log_info("Aggregating scHPF scores across neighborhoods...")
    
    all_agg <- get_nhood_schpf_scores(milo,
                                      da_results = da_results,
                                      cell_scores = colData(milo) %>% as.data.frame(),
                                      sum_cols = unname(factors),
                                      res_path = res_path)
    
  }
  # PLOTS - QC
  #-----------------------------------------------------------------------------
  log_info("Plotting QC metrics...")
  plot_milo_qc(res = res, res_path = res_path)
  
  # PLOTS - NETWORKS
  #-----------------------------------------------------------------------------
  if(!is.null(da_results)){
    log_info("Plotting graphs...")
    plot_milo_res(milo, da_results, label_size = 3, res_path = res_path, spatial.fdr = 0.05)
  }
  
  # PLOTS - FACTORS
  #-------------------------------------------------------------------------------
  milo_agg_factors <- read_rds(file.path(res_path, "milo_agg_factors.rds"))
  p <- plot_nhood_dotplot(milo_agg_factors, res_path = res_path, factors = factors, height = 4.5)
  write_xlsx(milo_agg_factors, file.path(res_path, "milo_agg_factors.xlsx"))
  
  saveRDS(milo, file.path(res_path, "milo.rds"))
  
  rm(list = c("res", "milo", "da_results", "res_path", "all_agg", "data", 
              "design", "formula", "keep_cells", "p")); gc()
}
#===============================================================================
#
# PANELS D, H
#
#===============================================================================
vars <- c('amyloid', 'tangles')
names(vars) <- c('High Amyloid', 'High Tangles')

for(i in 1:length(vars)){
  compared_nhoods(run1_path = file.path(home_path, "proj_snuc/analysis/pathoAD_k30_deltaFC2"),
                  name1 = "Pathological AD",
                  run2_path = file.path(home_path, "proj_snuc/analysis", paste0(vars[i], "_k30_deltaFC2")),
                  name2 = names(vars)[i],
                  file_name = paste0(fig, "_pathoAD_vs_", vars[i], ".pdf"),
                  res_path = file.path(home_path, "proj_snuc/analysis/"))
}
