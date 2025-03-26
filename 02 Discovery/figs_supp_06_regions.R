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
library(ComplexHeatmap)
library(grid)
library(logger)
library(miloR)
library(patchwork)
library(scater)
library(scCustomize)
library(schex)
library(scran)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "miloR_REGIONS")
fig <- "fig_supp06_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_miloR.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))

for(i in c("AWS", "BA20_21", "BA9_46", "H")){
  tmp <- data@meta.data %>% mutate(tmp = ifelse(REGIONS == i, 1, 0)) %>% pull(tmp)
  names(tmp) <- data@meta.data$cell_names
  data <- AddMetaData(data, tmp, i)
}

#===============================================================================
#
# PANEL A
#
#===============================================================================
regions <- c("AWS", "BA20_21", "BA9_46", "H")
names(regions) <- c("AWS", "BA20/21", "BA9/46", "H")

plots <- list()
for(i in regions){
  
  keep_cells <- data@meta.data %>%filter(REGIONS %in% i & DX_CAT == "AD")  %>% pull(cell_names)
  tmp <- make_hexbin(subset(data, subset = cell_names %in% keep_cells), nbins = 100, dimension_reduction = "schpf.umap")
  
  p <- ggplot() +
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(tmp@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = expression(italic("N")*" cells"))+
    labs(title = names(regions)[regions %in% i])+
    theme_umap()+
    theme(plot.title = element_text(size = 14, family = "Helvetica", margin = margin(b = 0.05, unit = "in")),
          plot.margin = margin(t = 0.05, r = 0, b = 0, l = 0.05, unit = "in"),
          legend.title = element_text(size = 12, family = "Helvetica", margin = margin(b = 0.1, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.ticks = element_line(linewidth = 0.5, color = "white"),
          legend.ticks.length = unit(0.05, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 1)

png(file.path(res_path, paste0(fig, "panel_A.png")), height = 2, width = 8, res = 600, units = "in")
print(p)
dev.off()
#===============================================================================
#
# PANEL B - DEMOGRAPHICS
#
#===============================================================================
regions <- c("AWS", "BA20_21", "BA9_46", "H")
names(regions) <- c("AWS", "BA20/21", "BA9/46", "H")

sex_cols <- unname(proj_cols("yellow", "blue"))
names(sex_cols) <- c("M", "F")

p1 <- data@meta.data %>% 
  filter(DX_CAT == "AD" & REGIONS %in% c("BA9_46","AWS", "BA20_21","H")) %>%
  group_by(REGIONS, SEX) %>%
  summarise(n=n()) %>%
  group_by(REGIONS) %>%
  summarise(SEX, n, pct=n/sum(n)) %>%
  ggplot(aes(y = pct, x = REGIONS, fill = SEX)) +
  geom_col()+
  geom_text(aes(label = formatC(n, big.mark = ",")), position = position_stack(vjust = .5), color = "white") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     label = function(x) formatC(x, format = "g"))+
  scale_fill_manual(name = "Sex", values = sex_cols) +
  scale_x_discrete(breaks = regions, labels = names(regions))+
  labs(y = "Proportion of cells (%)")+
  theme_proj()+
  theme(panel.grid.major = element_blank(),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'top',
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "in"),
        legend.margin = margin(b = -0.1, unit = "in"))

age_cols <- unname(proj_cols("purple", "blue", "teal", "yellow"))
names(age_cols) <- c("<=75", "76-85", "86-95", ">95")

data@meta.data <- data@meta.data %>%
  mutate(AGE_CAT = ifelse(AGE <=75, "<=75", NA),
         AGE_CAT = ifelse(AGE >75 & AGE <=85, "76-85", AGE_CAT),
         AGE_CAT = ifelse(AGE >85 & AGE <=95, "86-95", AGE_CAT),
         AGE_CAT = ifelse(AGE >95, ">95", AGE_CAT),
         AGE_CAT = factor(AGE_CAT, levels = c("<=75", "76-85", "86-95", ">95")))

p2 <- data@meta.data %>% 
  filter(DX_CAT == "AD" & REGIONS %in% c("BA9_46","AWS", "BA20_21","H")) %>%
  group_by(REGIONS, AGE_CAT) %>%
  summarise(n=n()) %>%
  group_by(REGIONS) %>%
  summarise(AGE_CAT, n, pct=n/sum(n)) %>%
  ggplot(aes(y = pct, x = REGIONS, fill = AGE_CAT)) +
  geom_col()+
  geom_text(aes(label = formatC(n, big.mark = ",")), position = position_stack(vjust = .5), color = "white") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     label = function(x) formatC(x, format = "g"))+
  scale_fill_manual(name = "Age\n(years)", values = age_cols) +
  scale_x_discrete(breaks = regions, labels = names(regions))+
  labs(y = "Proportion of cells (%)")+
  theme_proj()+
  theme(panel.grid.major = element_blank(),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = 'top',
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "in"),
        legend.margin = margin(b = -0.1, unit = "in"))

p <- wrap_plots(p1, p2)

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 3, width = 7.5)
print(p)
dev.off()
#===============================================================================
#
# miloR
#
#===============================================================================
data@meta.data$SEX <- ifelse(data@meta.data$SEX == "M", 1, 0)

for(trait in c("AWS", "BA20_21", "H")){ 
  
  log_info("Trait: ", trait)
  k <- 30
  deltaFC <- 1
  
  res_path <- file.path(home_path, "miloR_REGIONS", paste0(trait, "_k", k, "_deltaFC", deltaFC))
  if(!dir.exists(res_path)){dir.create(res_path)}
  if(file.exists(file.path(res_path, "log.txt"))) {file.remove(file.path(res_path, "log.txt"))}
  log_appender(appender_tee(file.path(res_path, "log.txt")))
  
  # DATA
  #-----------------------------------------------------------------------------
  milo <- data
  
  # CREATE OBJECT
  #-----------------------------------------------------------------------------
  log_info("Subsetting cells...")
  
  set.seed(1234)
  keep_cells <- milo@meta.data %>% 
    filter(DX_CAT == "AD" & REGIONS %in% c("BA9_46", "AWS", "BA20_21", "H")) %>%
    dplyr::select(cell_names, Donor.ID, REGIONS, AGE, SEX) %>% 
    drop_na() %>%
    group_by(Donor.ID, REGIONS) %>%
    slice_sample(n = 500) %>%
    pull(cell_names)
  
  keep_cells <- milo@meta.data %>%
    filter(cell_names %in% keep_cells & REGIONS %in% c("BA9_46", trait)) %>%
    pull(cell_names)
  
  milo <- subset(milo, subset = cell_names %in% keep_cells)
  milo@meta.data$new_id <- paste0(gsub(" ", "", milo@meta.data$Donor.ID), "_", milo@meta.data$REGIONS)
  
  milo@meta.data <- milo@meta.data %>% mutate(trait = ifelse(REGIONS == trait, 1, 0))
  meta <- milo@meta.data[,c("new_id", "SEX", "AGE", "ORIGIN", "trait")]
  
  log_info(trait, " counts: ", paste0(names(table(meta$trait)), " (", table(meta$trait), ")", collapse = ", "))
  
  invisible(lapply(capture.output(milo), log_info))
  
  # CREATE MILO OBJECT
  #-------------------------------------------------------------------------------
  log_info("Creating Milo object...")
  milo <- Milo(as.SingleCellExperiment(milo))
  invisible(lapply(capture.output(milo), log_info))
  
  set.seed(1234)
  milo <- buildGraph(milo, k = k, d = 23, reduced.dim = "SCHPF")
  
  set.seed(1234)
  milo <- makeNhoods(milo, prop = 0.1, k = k, d = 23, refined = TRUE, reduced_dim = "SCHPF", refinement_scheme = "graph")
  
  plot_milo_nhood_size(milo, bins = 50, res_path = res_path)
  
  # GET NEIGHBORHOODS
  #-------------------------------------------------------------------------------
  log_info("Calculating neighborhoods...")
  
  set.seed(1234)
  milo <- countCells(milo, meta.data = meta, samples = "new_id")
  
  set.seed(1234)
  milo <- calcNhoodDistance(milo, d = 23, reduced.dim = "SCHPF")
  
  # EXPERIMENTAL DESIGN
  #-------------------------------------------------------------------------------
  log_info("Setting up experimental design...")
  design <- meta
  design$trait <- as.factor(design$trait) 
  rownames(design) <- NULL
  design <- distinct(design)
  rownames(design) <- design$new_id
  
  # RUN DA
  #-------------------------------------------------------------------------------
  log_info("Performing DA testing...")
  formula <- as.formula(paste0("~", paste0(names(meta)[-1], collapse = " + ")))
  
  log_info("Formula: ", paste0(gsub("trait", trait, as.character(formula)), collapse = ""))
  
  res <- testNhoods(milo, 
                    design = formula, 
                    design.df = design,
                    reduced.dim = "SCHPF", 
                    fdr.weighting = "graph-overlap")
  
  saveRDS(res, file.path(res_path, "testNhoods.rds"))
  
  log_info("Calculating neighborhood adjacency...")
  set.seed(1234)
  milo <- buildNhoodGraph(milo, overlap = 3)
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
  #-------------------------------------------------------------------------------
  log_info("Getting neighborhood composition...")
  
  get_nhood_comp(milo, res_path)
  
  if(!is.null(da_results)){
    
    log_info("Aggregating scHPF scores across neighborhoods...")
    
    all_agg <- get_nhood_schpf_scores(milo,
                                      da_results = da_results,
                                      cell_scores = data@meta.data,
                                      sum_cols = unname(factors),
                                      res_path = res_path)
    
  }
  
  # PLOTS - QC
  #-------------------------------------------------------------------------------
  log_info("Plotting QC metrics...")
  plot_milo_qc(res = res, res_path = res_path)
  
  # PLOTS - NETWORKS
  #-------------------------------------------------------------------------------
  if(!is.null(da_results)){
    log_info("Plotting graphs...")
    plot_milo_res(milo, da_results, label_size = 3, res_path = res_path, spatial.fdr = 0.05)
  }
  
  # PLOTS - FACTORS
  #-------------------------------------------------------------------------------
  milo_agg_factors <- read_rds(file.path(res_path, "milo_agg_factors.rds"))
  plot_nhood_dotplot(milo_agg_factors, res_path = res_path, factors = factors)
  write_xlsx(milo_agg_factors, file.path(res_path, "milo_agg_factors.xlsx"))
  
  saveRDS(milo, file.path(res_path, "milo.rds"))
  invisible(lapply(capture.output(sessionInfo()), log_info))
  rm(list = ls()[!ls() %in% keep_vars])
  
}
