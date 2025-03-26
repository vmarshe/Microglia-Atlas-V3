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
library(scran)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(viridis)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "miloR_PMI_NYBB")
fig <- "fig_supp07_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_miloR.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))
data@meta.data$SEX <- ifelse(data@meta.data$SEX == "M", 1, 0)

#===============================================================================
#
# PANEL A - TIME TO DISSECTION
#
#===============================================================================
# save umap boundaries
minx <- min(data@reductions$schpf.umap@cell.embeddings[,1])
maxx <- max(data@reductions$schpf.umap@cell.embeddings[,1])
miny <- min(data@reductions$schpf.umap@cell.embeddings[,2])
maxy <- max(data@reductions$schpf.umap@cell.embeddings[,2])

set.seed(1234)
keep_cells <- data@meta.data %>% 
  filter(ORIGIN == "CUMC/NYBB" & !is.na(PMI)) %>% 
  pull(cell_names)

tmp <- subset(data, subset = cell_names %in% keep_cells)

p <- FeaturePlot(tmp, feature = "PMI", reduction = 'schpf.umap', raster = T) +
  scale_colour_gradientn(name = "Time to tissue\ndissection\n(hours)", colors = viridis_pal(option = "F")(10))+
  theme_umap() +
  theme(legend.text = element_text(size = 12, family = "Helvetica", margin = margin(0, 0, 0, 0.05, 'in')),
        legend.title = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, unit = "in"),
        plot.title = element_blank(),
        panel.border = element_rect(linewidth = 1, color = "black"),
        axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.25, color = "black"),
        legend.key.width  = unit(0.25, "in"), 
        legend.key.height = unit(0.3, "in"))

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 3, width = 3.25)
print(p)
dev.off()

#===============================================================================
#
# PANEL B - HISTOGRAM
#
#===============================================================================
set.seed(1234)
med <- median(data@meta.data %>% group_by(`Donor.ID`) %>% sample_n(1) %>% filter(ORIGIN == "CUMC/NYBB") %>% pull(PMI), na.rm=T)

set.seed(1234)
p <- data@meta.data%>% 
  filter(ORIGIN == "CUMC/NYBB" & !is.na(PMI)) %>% 
  group_by(`Donor.ID`) %>% 
  sample_n(1) %>%
  ggplot(aes(x = PMI)) + 
  geom_histogram(bins = 15, color = "black", fill = proj_cols("orange2")) + 
  geom_vline(xintercept = med, lty = 2, color = proj_cols('red')) +
  labs(x = "NYBB time-to dissection (hours)", y = "Frequency (N)")+
  theme_proj() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 2.5, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL C, D, E, F - MILO
#
#===============================================================================
trait <- "PMI"
k <- 30
deltaFC <- 2

if(file.exists(file.path(res_path, "log.txt"))) {file.remove(file.path(res_path, "log.txt"))}
log_appender(appender_tee(file.path(res_path, "log.txt")))
log_info("Trait: ", trait)

# CREATE OBJECT
#-----------------------------------------------------------------------------
log_info("Subsetting cells...")
milo <- data

set.seed(1234)
med <- milo@meta.data %>% 
  filter(ORIGIN == "CUMC/NYBB" & DX_CAT == "AD" & !is.na(PMI)) %>%
  group_by(`Donor.ID`) %>% 
  sample_n(1) %>% 
  pull(PMI) %>%
  median()

set.seed(1234)
keep_cells <- milo@meta.data %>%
  filter(ORIGIN == "CUMC/NYBB" & DX_CAT == "AD") %>%
  dplyr::select(Donor.ID, REGIONS, trait, cell_names, SEX, AGE) %>%
  drop_na() %>%
  group_by(Donor.ID, REGIONS) %>%
  slice_sample(n = 1000) %>%
  pull(cell_names)

milo <- subset(milo, subset = cell_names %in% keep_cells)
milo@meta.data$new_id <- paste0(gsub(" ", "", milo@meta.data$Donor.ID), "_", milo@meta.data$REGIONS)
milo@meta.data$trait <- ifelse(milo@meta.data$PMI > med, 1, 0)

meta <- milo@meta.data[,c("new_id", "SEX", "AGE", "REGIONS", "trait")]

invisible(lapply(capture.output(milo), log_info))
log_info(trait, " counts: ", paste0(names(table(meta$trait)), " (", table(meta$trait), ")", collapse = ", "))

# CREATE MILO OBJECT
#-------------------------------------------------------------------------------
log_info("Creating Milo object...")
milo <- Milo(as.SingleCellExperiment(milo))
invisible(lapply(capture.output(milo), log_info))

set.seed(1234)
milo <- buildGraph(milo, k = k, d = 23, reduced.dim = "SCHPF")

set.seed(1234)
milo <- makeNhoods(milo, prop = 0.1, k = k, d = 23, refined = TRUE, 
                   reduced_dim = "SCHPF", refinement_scheme = "graph")

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
formula <- as.formula(paste0("~ ", paste0(names(meta)[-1], collapse = " + ")))

log_info("Formula: ", paste0(gsub("trait", trait, as.character(formula)), collapse = ""))

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
#===============================================================================
#
# PANEL G - DENSITY PLOTS
#
#===============================================================================
set.seed(1234)
med <- data@meta.data %>% 
  filter(ORIGIN == "CUMC/NYBB" & DX_CAT == "AD" & !is.na(PMI)) %>%
  group_by(`Donor.ID`) %>% 
  sample_n(1) %>% 
  pull(PMI) %>%
  median()

keep_cells <- data@meta.data %>% filter(ORIGIN == "CUMC/NYBB" & DX_CAT == "AD" & PMI <= med) %>% pull(cell_names)
pmi_low <- subset(data, subset = cell_names %in% keep_cells)

keep_cells <- data@meta.data %>% filter(ORIGIN == "CUMC/NYBB" & DX_CAT == "AD" & PMI > med)  %>% pull(cell_names)
pmi_high <- subset(data, subset = cell_names %in% keep_cells)

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

keep_factors <- paste0("scHPF_", c(1, 23, 17, 21))

plots <- list()
for(i in keep_factors){
  
  p1 <- Plot_Density_Custom(pmi_low, features = i, reduction = "schpf.umap", pt.size = 0.1)
  p1$data$facet <- names(factors)[factors %in% i]
  p1 <- p1 + facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(limits = c(ymin, ymax))+
    theme_umap()+
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.title = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          strip.placement = "inside", 
          strip.background = element_rect(fill = NA, color = NA, linewidth = NA),
          strip.text = element_text(size = 8, family = "Helvetica"),
          plot.margin  = margin(0, 0, 0, 0, "in"))
  
  p2 <- Plot_Density_Custom(pmi_high, features = i, reduction = "schpf.umap", pt.size = 0.1)+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(limits = c(ymin, ymax))+
    theme_umap()+
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.title = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"))
  
  p <- wrap_plots(p1, p2) +  theme(plot.margin  = margin(0, 0, 0, 0, "in")) + plot_layout(ncol = 1)
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 1)

cairo_pdf(file.path(res_path,  paste0(fig, "panel_G.pdf")), height = 3.5, width = 8)
print(p)
dev.off()

#===============================================================================
#
# PANEL I - ROSMAP
#
#===============================================================================
res_path <- file.path(home_path, "miloR_PMI_ROSMAP")

minx <- min(data@reductions$schpf.umap@cell.embeddings[,1])
maxx <- max(data@reductions$schpf.umap@cell.embeddings[,1])
miny <- min(data@reductions$schpf.umap@cell.embeddings[,2])
maxy <- max(data@reductions$schpf.umap@cell.embeddings[,2])

set.seed(1234)
keep_cells <- data@meta.data %>% 
  filter(ORIGIN == "RUSH" & !is.na(PMI)) %>% 
  pull(cell_names)

tmp <- subset(data, subset = cell_names %in% keep_cells)

p <- FeaturePlot(tmp, feature = "PMI", reduction = 'schpf.umap', raster = T) +
  scale_colour_gradientn(name = "PMI\n(hours)", colors = viridis_pal(option = "F")(10))+
  theme_umap() +
  theme(legend.text = element_text(size = 12, family = "Helvetica", margin = margin(0, 0, 0, 0.05, 'in')),
        legend.title = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, unit = "in"),
        plot.title = element_blank(),
        panel.border = element_rect(linewidth = 1, color = "black"),
        axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.25, color = "black"),
        legend.key.width  = unit(0.25, "in"), 
        legend.key.height = unit(0.3, "in"))

pdf(file.path(res_path, paste0(fig, "panel_H.pdf")), height = 3, width = 3.25)
print(p)
dev.off()

#===============================================================================
#
# PANEL I
#
#===============================================================================
set.seed(1234)
med <- median(data@meta.data %>% group_by(`Donor.ID`) %>% sample_n(1) %>% filter(ORIGIN == "RUSH") %>% pull(PMI), na.rm=T)

set.seed(1234)
p <- data@meta.data%>% 
  filter(ORIGIN == "RUSH" & !is.na(PMI)) %>% 
  group_by(`Donor.ID`) %>% 
  sample_n(1) %>%
  ggplot(aes(x = PMI)) + 
  geom_histogram(bins = 15, color = "black", fill = proj_cols("orange2")) + 
  geom_vline(xintercept = med, lty = 2, color = proj_cols('red')) +
  labs(x = "ROSMAP PMI (hours)", y = "Frequency (N)")+
  theme_proj() + 
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

pdf(file.path(res_path, paste0(fig, "panel_I.pdf")), height = 2.5, width = 3.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL K, L, M - MILO
#
#===============================================================================
trait <- "PMI"
k <- 30
deltaFC <- 2

if(file.exists(file.path(res_path, "log.txt"))) {file.remove(file.path(res_path, "log.txt"))}
log_appender(appender_tee(file.path(res_path, "log.txt")))
log_info("Trait: ", trait)

# CREATE OBJECT
#-----------------------------------------------------------------------------
log_info("Subsetting cells...")
milo <- data

set.seed(1234)
med <- milo@meta.data %>% 
  filter(ORIGIN == "RUSH" & DX_CAT == "AD" & !is.na(PMI)) %>%
  group_by(`Donor.ID`) %>% 
  sample_n(1) %>% 
  pull(PMI) %>%
  median()

set.seed(1234)
keep_cells <- milo@meta.data %>%
  filter(ORIGIN == "RUSH" & DX_CAT == "AD") %>%
  dplyr::select(Donor.ID, REGIONS, trait, cell_names, SEX, AGE) %>%
  drop_na() %>%
  group_by(Donor.ID, REGIONS) %>%
  slice_sample(n = 1000) %>%
  pull(cell_names)

milo <- subset(milo, subset = cell_names %in% keep_cells)
milo@meta.data$new_id <- paste0(gsub(" ", "", milo@meta.data$Donor.ID), "_", milo@meta.data$REGIONS)
milo@meta.data$trait <- ifelse(milo@meta.data$PMI > med, 1, 0)

meta <- milo@meta.data[,c("new_id", "SEX", "AGE", "REGIONS", "trait")]

invisible(lapply(capture.output(milo), log_info))
log_info(trait, " counts: ", paste0(names(table(meta$trait)), " (", table(meta$trait), ")", collapse = ", "))

# CREATE MILO OBJECT
#-------------------------------------------------------------------------------
log_info("Creating Milo object...")
milo <- Milo(as.SingleCellExperiment(milo))
invisible(lapply(capture.output(milo), log_info))

set.seed(1234)
milo <- buildGraph(milo, k = k, d = 23, reduced.dim = "SCHPF")

set.seed(1234)
milo <- makeNhoods(milo, prop = 0.1, k = k, d = 23, refined = TRUE, 
                   reduced_dim = "SCHPF", refinement_scheme = "graph")

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
formula <- as.formula(paste0("~ ", paste0(names(meta)[-1], collapse = " + ")))

log_info("Formula: ", paste0(gsub("trait", trait, as.character(formula)), collapse = ""))

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
  
  p <- plot_nhood_feature(nhood_meta_data = all_agg, features = factors)
  p <- wrap_plots(p) + 
    plot_layout(nrow = 4) & 
    theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, 'in'))
  
  pdf(file.path(res_path,  "milo_NhoodGroup_factors.pdf"), height = 7, width = 20)
  print(p)
  dev.off()
  
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

#===============================================================================
#
# PANEL N
#
#===============================================================================
set.seed(1234)
med <- data@meta.data %>% 
  filter(ORIGIN == "RUSH" & DX_CAT == "AD" & !is.na(PMI)) %>%
  group_by(`Donor.ID`) %>% 
  sample_n(1) %>% 
  pull(PMI) %>%
  median()

keep_cells <- data@meta.data %>% filter(ORIGIN == "RUSH" & DX_CAT == "AD" & PMI <= med) %>% pull(cell_names)
pmi_low <- subset(data, subset = cell_names %in% keep_cells)

keep_cells <- data@meta.data %>% filter(ORIGIN == "RUSH" & DX_CAT == "AD" & PMI > med)  %>% pull(cell_names)
pmi_high <- subset(data, subset = cell_names %in% keep_cells)

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

keep_factors <- paste0("scHPF_", c(9, 1, 23))

plots <- list()
for(i in keep_factors){
  
  p1 <- Plot_Density_Custom(pmi_low, features = i, reduction = "schpf.umap", pt.size = 0.1)
  p1$data$facet <- names(factors)[factors %in% i]
  p1 <- p1 + facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(limits = c(ymin, ymax))+
    theme_umap()+
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.title = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          strip.placement = "inside", 
          strip.background = element_rect(fill = NA, color = NA, linewidth = NA),
          plot.margin  = margin(0, 0, 0, 0, "in"))
  
  p2 <- Plot_Density_Custom(pmi_high, features = i, reduction = "schpf.umap", pt.size = 0.1)+
    scale_x_continuous(limits = c(xmin, xmax))+
    scale_y_continuous(limits = c(ymin, ymax))+
    theme_umap()+
    theme(legend.position = "right",
          legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.2, unit = "in"),
          legend.title = element_blank(),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"))
  
  p <- wrap_plots(p1, p2) +  theme(plot.margin  = margin(0, 0, 0, 0, "in")) + plot_layout(ncol = 1)
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 1)

pdf(file.path(res_path,  paste0(fig, "panel_N.pdf")), height = 3.5, width = 6.5)
print(p)
dev.off()
