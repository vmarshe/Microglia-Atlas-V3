#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(colorspace)
library(ggraph)
library(ggrastr)
library(ggtext)
library(igraph)
library(logger)
library(miloR)
library(patchwork)
library(scater)
library(scCustomize)
library(schex)
library(scran)
library(Seurat)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "/mnt/vast/hpc/homes/vsm2116/projects/microglia_atlas"
res_path <- file.path(home_path, "miloR_SEX")
fig <- "fig_supp04_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_miloR.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))

#===============================================================================
#
# PANEL A
#
#===============================================================================
plots <- list()
for(i in c("M", "F")){
  
  keep_cells <- data@meta.data %>% 
    filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H") & SEX == i) %>% 
    pull(cell_names)
  
  tmp <- make_hexbin(subset(data, subset = cell_names %in% keep_cells), nbins = 100, dimension_reduction = "schpf.umap")
  
  p <- ggplot() +
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(tmp@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "Number\nof cells per\nhexbin")+
    theme_umap()+
    theme(legend.title = element_text(size = 12, family = "Helvetica", margin = margin(b = 0.05, unit = "in")),
          plot.margin = margin(0, 0, 0, 0),
          legend.margin = margin(l = -0.1, unit = "in"),
          legend.key.height = unit(0.3, "in"),
          legend.key.width = unit(0.25, "in"),
          legend.ticks = element_line(linewidth = 1, color = "white"),
          legend.ticks.length = unit(0.05, "in"),
          legend.text = element_text(size = 14, margin = margin(l = 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(ncol = 2)

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 3, width = 7)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
sex_cols <- unname(proj_cols("yellow", "blue"))
names(sex_cols) <- c("M", "F")

p <- data@meta.data %>% 
  filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H")) %>%
  dplyr::select(cell_names, Donor.ID, REGIONS, SEX, DX_CAT, AGE) %>% 
  drop_na() %>%
  group_by(SEX, REGIONS) %>% 
  summarise(n=n()) %>% 
  group_by(REGIONS) %>% 
  summarise(SEX, n, pct=n/sum(n)) %>% 
  mutate(REGIONS = recode(REGIONS, "BA9_46"="BA9/46", "BA20_21"="BA20/21"),
         REGIONS = factor(REGIONS, levels = rev(c("AWS", "BA9/46", "BA20/21", "H")))) %>%
  ggplot(aes(x = pct, y = REGIONS, fill = SEX)) +
  geom_col(color = "black", alpha = 0.9) +
  geom_text(aes(label = formatC(n, big.mark = ",")), position = position_stack(vjust = .5)) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) formatC(x, format = "g"))+
  scale_fill_manual(values = sex_cols) + 
  labs(y = "Tissue", x = "Proportion")+
  theme_proj() + 
  theme(panel.grid.major =   element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.2, 'in'),
        legend.direction = "horizontal", 
        legend.position = "top", 
        legend.margin = margin(b = -0.15, unit = "in"))

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 3.5, width = 4)
print(p)
dev.off()

#===============================================================================
#
# PANELS C-F
#
#===============================================================================
trait <- "SEX"
k <- 30
deltaFC <- 1

# PATHS
#-------------------------------------------------------------------------------
res_path <- file.path(home_path, "miloR_SEX", paste0(trait, "_k", k, "_deltaFC", deltaFC))
if(!dir.exists(res_path)) dir.create(res_path)
if(file.exists(file.path(res_path, "log.txt"))) {file.remove(file.path(res_path, "log.txt"))}
log_appender(appender_tee(file.path(res_path, "log.txt")))

# DATA
#-------------------------------------------------------------------------------
log_info("Trait: ", trait)
log_info("Subsetting cells...")

milo <- data
milo@meta.data$SEX <- ifelse(milo@meta.data$SEX == "M", 1, 0)
milo@meta.data$trait <- milo@meta.data$SEX

set.seed(1234)
keep_cells <- milo@meta.data %>% 
  filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H")) %>%
  dplyr::select(cell_names, Donor.ID, REGIONS, trait, DX_CAT, AGE) %>% 
  drop_na() %>%
  group_by(Donor.ID, REGIONS) %>%
  slice_sample(n = 500) %>%
  pull(cell_names)

milo <- subset(milo, subset = cell_names %in% keep_cells)
milo@meta.data$new_id <- paste0(gsub(" ", "", milo@meta.data$Donor.ID), "_", milo@meta.data$REGIONS)
milo@meta.data <- milo@meta.data %>% mutate(trait = SEX)
meta <- milo@meta.data[,c("new_id", "AGE", "REGIONS", "ORIGIN", "trait")]

log_info(trait, " counts: ", paste0(names(table(meta$trait)), " (", table(meta$trait), ")", collapse = ", "))

invisible(lapply(capture.output(data), log_info))

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
res_path <- file.path(home_path, "miloR_SEX")

milo <- read_rds(file.path(res_path, "SEX_k30_deltaFC2/milo.rds"))
milo_res <- read_rds(file.path(res_path, "SEX_k30_deltaFC2/milo_grouped_Nhoods.rds"))
comparison <- read_rds(file.path(res_path, "SEX_k30_deltaFC2/comparison_test.rds")) %>% 
  mutate(dir = ifelse(t > 0, "up", "down"))

colData(milo)[unlist(nhoodIndex(milo)[milo_res$Nhood]),"NhoodGroup"] <- milo_res$NhoodGroup
colData(milo)[unlist(nhoodIndex(milo)[milo_res$Nhood]),"SpatialFDR"] <- milo_res$SpatialFDR
colData(milo)[unlist(nhoodIndex(milo)[milo_res$Nhood]),"logFC"] <- milo_res$logFC

nh_graph <- nhoodGraph(milo)

layout <- reducedDim(milo, "SCHPF.UMAP")[as.numeric(vertex_attr(nh_graph)$name),] %>% as.matrix()

order <- sort(as.numeric(unique(colData(milo)[, "NhoodGroup"])))
col_vals <- as.character(colData(milo)[as.numeric(vertex_attr(nh_graph)$name), "NhoodGroup"])
col_vals <- factor(col_vals, levels = order)
V(nh_graph)$NhoodGroup <- col_vals

V(nh_graph)$SpatialFDR <- colData(milo)[as.numeric(vertex_attr(nh_graph)$name), "SpatialFDR"]
V(nh_graph)$logFC <- colData(milo)[as.numeric(vertex_attr(nh_graph)$name), "logFC"]

g.df <- create_layout(nh_graph, layout = layout)
g.df[, "NhoodGroup"] <- vertex_attr(nh_graph)[["NhoodGroup"]]

pl <- ggraph(simplify(nh_graph), layout = layout) 
pl$data$dir <- ifelse(pl$data$logFC > 0, "up", "down" )

set.seed(1234)
group_cols <- lighten(sample(pal()(length(order))), 0.2, space = "HLS")

keep_nhood <- comparison  %>% 
  filter(p.adj < 0.05) %>% 
  arrange(dir, desc(abs(t))) %>% 
  group_by(dir) %>% 
  top_n(2, wt = abs(diff)) %>%
  pull(nhood)

nhood_plots <- list()
for(i in as.numeric(keep_nhood)){
  tmp <- comparison %>% filter(nhood == i) 
  tmp_data <- pl$data %>% filter(SpatialFDR < 0.05 & NhoodGroup == i) %>% filter(dir == tmp$dir)
  
  n_tot <- nrow(pl$data %>% filter(NhoodGroup == i))
  n_sig <- nrow(tmp_data)
  lab <- paste0(n_sig, "/", n_tot, " (", round(n_sig/n_tot*100, 1), "%)")
  
  p <- ggplot() + 
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_node_point(data = pl$data %>% filter(SpatialFDR > 0.05 & NhoodGroup == i), 
                    aes(fill = NhoodGroup, size = size), 
                    alpha = 0.8, 
                    shape = 21, 
                    stroke = 0.3, 
                    color = proj_cols_grey("med grey"),
                    fill = proj_cols_grey("light grey")) +
    geom_node_point(data = tmp_data,
                    aes(fill = NhoodGroup, size = size), 
                    alpha = 1, 
                    shape = 21, 
                    stroke= 0.3,
                    color = "black", 
                    fill = group_cols[i]) +
    annotate("text", label = lab, x = Inf, y = -Inf, vjust = -0.2, hjust = 1, family = "Helvetica", size = 3.5)+
    annotate("text", label = paste0("M-Nhood ", i), x = -Inf, y = Inf, vjust = 1, hjust = 0, family = "Helvetica", size = 3.5, 
             fontface = "bold")+
    scale_size(range =c(0.2, 2), name="Nhood size") +
    theme_umap()+ 
    theme(legend.position = "none")  
  
  nhood_plots <- c(nhood_plots, list(p))
}

male <- subset(data, subset = SEX == "M" & DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H"))
female <-  subset(data, subset = SEX == "F" & DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H"))

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

plots <- list()
for(i in c(9, 1, 10, 3)){
  
  ind <- which(factors %in% paste0("scHPF_", i))
  
  p1 <- Plot_Density_Custom(male, features = paste0("scHPF_", i), reduction = "schpf.umap")
  p1$data$facet <- "M"
  p1 <- ggplot()+
    geom_point_rast(data = p1$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5)+
    facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(name = names(factors)[ind], limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    scale_color_gradientn(colors = viridis_magma_light_high)+
    theme_umap()+
    theme(strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in"), family = "Helvetica"),
          strip.background = element_rect(fill = lighten(sex_cols["M"], 0.6) , color = sex_cols["M"], linewidth = 1.5),
          axis.title.y = element_textbox_simple(family = "Helvetica", face = "bold", size = 12, orientation = "left-rotated",
                                                hjust = 1, halign = 0.5, box.color = "#F0F0F0", fill = "#F0F0F0", padding = margin(2.5, 2.5, 2.5, 2.5),
                                                margin = margin(0, 0, -0.05, 0, unit = "in")),
          plot.title = element_blank(),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_blank(),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0, 0, 0, 0, unit = "in"),
          legend.key.height = unit(0.25, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(t=0, r=0, b=0, l=0.05, 'in'), family = "Helvetica"),
          legend.title = element_blank(),
          legend.margin = margin(t=0, r=0, b=0, l=-0.2, 'in'),
          legend.ticks = element_line(linewidth = 0.5, color = 'white'),
          legend.ticks.length = unit(0.05, "in"))
  
  p2 <- Plot_Density_Custom(female, features = paste0("scHPF_", i), reduction = "schpf.umap")
  p2$data$facet <- "F"
  p2 <- ggplot()+
    geom_point_rast(data = p2$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5)+
    facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(name = names(factors)[ind], limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    scale_color_gradientn(colors = viridis_magma_light_high)+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    theme_umap()+
    theme(strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in"), family = "Helvetica"),
          strip.background = element_rect(fill = lighten(sex_cols["F"], 0.6) , color = sex_cols["F"], linewidth = 1.5),
          plot.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "in"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.title.y = element_blank(),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          legend.key.height = unit(0.25, "in"),
          legend.key.width = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(t=0, r=0, b=0, l=0.05, 'in'), family = "Helvetica"),
          legend.title = element_blank(),
          legend.margin = margin(t=0, r=0, b=0, l=-0.2, 'in'),
          legend.ticks = element_line(linewidth = 0.5, color = 'white'),
          legend.ticks.length = unit(0.05, "in"))
  
  plots <- c(plots, list(wrap_plots(p1, p2) + plot_layout(nrow = 1)))
}
cairo_pdf(file.path(res_path, paste0(fig, "panel_H.pdf")), height = 5.5, width = 15)
wrap_plots(nhood_plots[[1]], plots[[1]], 
           nhood_plots[[2]], plots[[2]],
           nhood_plots[[3]], plots[[3]],
           nhood_plots[[4]], plots[[4]]) + 
  plot_layout(widths = c(0.15, 0.35, 0.15, 0.35), nrow = 2, ncol = 4)
dev.off()
