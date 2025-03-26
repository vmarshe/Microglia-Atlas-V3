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
library(colorspace)
library(ComplexHeatmap)
library(ggraph)
library(ggrastr)
library(igraph)
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

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "proj_snuc/analysis")
fig <- "fig_s17"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))

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

for(i in c("projid", "study", "msex", "pmi",  "age_death", "educ", "pathoAD", 
           "amyloid", "tangles",  "amyloid_sqrt", "tangles_sqrt", "niareagansc")){
  
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
# PANEL A-D
#
#===============================================================================
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

trait <- "pathoAD"
k <- 30
deltaFC <- 2

# PATHS
#-------------------------------------------------------------------------------
res_path <- file.path(home_path, "proj_snuc", "analysis", paste0(trait, "_k", k, "_deltaFC", deltaFC))

# DATA
#-------------------------------------------------------------------------------
log_info("Reading data...")
data <- DietSeurat(snuc, assays = "RNA", dimreducs = c("scHPF", "schpf.umap"))
Idents(data) <- "cell_names"

# CREATE OBJECT
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------
log_info("Calculating neighborhoods...")
set.seed(1234)
milo <- countCells(milo,
                   meta.data = data@meta.data[,c("new_id", "trait", "study","msex", "age_death", "pmi", "educ")],
                   samples = "new_id")
set.seed(1234)
milo <- calcNhoodDistance(milo, d = 23, reduced.dim = "SCHPF")

# EXPERIMENTAL DESIGN
#-------------------------------------------------------------------------------
log_info("Setting up experimental design...")
design <- data@meta.data[,c("new_id", "trait", "study", "msex", "age_death", "pmi", "educ")]
design$trait <- as.factor(design$trait)
rownames(design) <- NULL
design <- distinct(design)
rownames(design) <- design$new_id

# RUN DA
#-------------------------------------------------------------------------------
log_info("Performing DA testing...")
formula <- as.formula("~ study + age_death + msex + pmi + trait")

log_info("Formula: ", paste0(as.character(formula), collapse = ""))
res <- testNhoods(milo, design = formula, design.df = design,
                  reduced.dim = "SCHPF", fdr.weighting = "graph-overlap")
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
                                    cell_scores = colData(milo) %>% as.data.frame(),
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
p <- plot_nhood_dotplot(milo_agg_factors, res_path = res_path, factors = factors, height = 4.5)
write_xlsx(milo_agg_factors, file.path(res_path, "milo_agg_factors.xlsx"))

saveRDS(milo, file.path(res_path, "milo.rds"))
#===============================================================================
#
# PANEL E
#
#===============================================================================
res_path <- "~/projects/microglia_atlas/proj_snuc/analysis/pathoAD_k30_deltaFC2"
hex <- make_hexbin(snuc, nbins = 100, dimension_reduction = "schpf.umap")
milo <- read_rds(file.path(res_path, "milo.rds"))
milo_res <- read_rds(file.path(res_path, "milo_grouped_Nhoods.rds"))
comparison <- read_rds(file.path(res_path, "comparison_test.rds")) %>% mutate(dir = ifelse(t > 0, "up", "down"))

colData(milo)[unlist(nhoodIndex(milo)[milo_res$Nhood]),"NhoodGroup"] <- milo_res$NhoodGroup
colData(milo)[unlist(nhoodIndex(milo)[milo_res$Nhood]),"SpatialFDR"] <- milo_res$SpatialFDR
colData(milo)[unlist(nhoodIndex(milo)[milo_res$Nhood]),"logFC"] <- milo_res$logFC

nh_graph <- nhoodGraph(milo)

layout <- reducedDim(milo, "SCHPF.UMAP")[as.numeric(vertex_attr(nh_graph)$name),] %>% as.matrix()

order <- sort(as.numeric(unique(colData(milo)[, "NhoodGroup"])))
col_vals <- as.character(colData(milo)[as.numeric(vertex_attr(nh_graph)$name), "NhoodGroup"])
col_vals <- factor(col_vals, levels=order)
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
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
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
  
  nhood_plots = c(nhood_plots, list(p))
}

ad <- subset(snuc, subset = pathoAD == 1)
non_ad <- subset(snuc, subset = pathoAD == 0)

xmin <- min(ref@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(ref@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(ref@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(ref@reductions$schpf.umap@cell.embeddings[,2])

cols <- proj_cols("yellow", "purple")
names(cols) <- c("non-AD", "AD")

plots <- list()
for(i in c(1, 10, 26, 25)){
  
  ind <- which(factors %in% paste0("scHPF_", i))
  
  p1 <- Plot_Density_Custom(non_ad, features = paste0("scHPF_", i), pt.size = 0.5)
  p1$data$facet = "non-AD"
  p1 <- ggplot()+
    geom_point_rast(data = p1$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5)+
    facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(name = names(factors)[ind], limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    scale_color_gradientn(colors = viridis_magma_light_high)+
    theme_umap()+
    theme(strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in"), family = "Helvetica"),
          strip.background = element_rect(fill = lighten(cols["non-AD"], 0.6) , color = cols["non-AD"], linewidth = 1.5),
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
  
  p2 <- Plot_Density_Custom(ad, features = paste0("scHPF_", i), pt.size = 0.5)
  p2$data$facet = "AD"
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
          strip.background = element_rect(fill = lighten(cols["AD"], 0.6) , color = cols["AD"], linewidth = 1.5),
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

cairo_pdf(file.path(res_path, "panel_E.pdf"), height = 5, width = 14)
wrap_plots(nhood_plots[[1]], plots[[1]], 
           nhood_plots[[2]], plots[[2]],
           nhood_plots[[3]], plots[[3]],
           nhood_plots[[4]], plots[[4]]) + 
  plot_layout(widths = c(0.15, 0.35, 0.15, 0.35), nrow = 2, ncol = 4)
dev.off()


#===============================================================================
#
# PANEL H
#
#===============================================================================

# GENES
#-------------------------------------------------------------------------------
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

snuc <- NormalizeData(snuc, normalization.method = "LogNormalize", assay = "RNA") %>%
  ScaleData(assay = "RNA", features = rownames(schpf@feature.loadings))

genes <- data.frame()
for(i in factors){
  genes <- rbind(genes, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}

keep_factors <- paste0("scHPF_", c(1, 19, 25, 26))

keep_highlights <- data.frame()
for(i in keep_factors){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}

# PSEUDOBULK
#-------------------------------------------------------------------------------
snuc@meta.data <- snuc@meta.data %>%
  mutate(var = recode(as.character(niareagansc), "3"="3/4", "4"="3/4"),
         var = factor(var, levels = c("1", "2", "3/4")))

projids <- snuc@meta.data %>% dplyr::select(projid, var) %>% group_by(var) %>% arrange(projid) %>% distinct()
averages <- AverageExpression(snuc,
                              features = genes$name,
                              group.by = "projid",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data
averages <- t(averages)[projids$projid, genes$name]

keep_highlights <- keep_highlights %>% filter(name %in% colnames(averages))

var_cols <- unname(proj_cols("blue", "yellow", "pink"))
names(var_cols) <- levels(projids$var)

col_split <- genes %>% filter(name %in% colnames(averages)) %>% pull(factor_name)
col_split <- factor(col_split, levels = unique(genes$factor_name))

row_split <- projids$var
row_split <- factor(row_split, levels = c("1", "2", "3/4"), labels = c("High\n(1)", "Intermediate\n(2)", "Low/No AD\n(3/4)"))

bottom_ha <- HeatmapAnnotation(labs = anno_mark(at = which(colnames(averages) %in% keep_highlights$name),
                                               labels = keep_highlights$name,
                                               side ="bottom",
                                               labels_gp = gpar(fontface = "italic", family = "Helvetica", fontsize = 8)),
                              which = "column")

fill_cols <- unname(proj_cols("blue", "yellow", "pink"))

cairo_pdf(file.path(home_path, "proj_snuc", "analysis", paste0(fig, "panelF.pdf")), height = 5, width = 9)
set.seed(1234)
p <- Heatmap(averages, name = "mat",
             col = colorRamp2(breaks = c(min(averages), 0, max(averages)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             left_annotation = rowAnnotation(row = anno_empty(border = FALSE, width = unit(0.1, "in"))),
             cluster_rows = T,
             show_row_names = FALSE,
             show_row_dend = FALSE,
             row_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 8),
             row_split = row_split,
             row_gap = unit(0.05, "in"),
             top_annotation = HeatmapAnnotation(foo = anno_empty(border = FALSE, height = unit(0.1, "in")), which = "column"),
             column_split = col_split,
             column_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 8),
             column_title_rot = 90,
             column_gap = unit(0.02, "in"),
             bottom_annotation = bottom_ha,
             cluster_columns = T,
             show_column_dend = FALSE,
             show_column_names = FALSE,
             heatmap_legend_param = list(
               title = expression('Avg. scaled exp.'),
               legend_width = unit(1.25, "in"),
               direction = "horizontal",
               title_position = "lefttop",
               border = "black"
             ))

draw(p,
     background = "transparent",
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     padding = unit(c(b=0.2, l=0.2, t=0.2, r=1), "in"),
     legend_title_gp = gpar(fontsize = 10, hjust = 0.5, fontfamily = "Helvetica"))

for(i in 1:length(factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, height = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}
for(i in 1:3){
  decorate_annotation("row", slice = i, {
    grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = fill_cols[i], col = NA), just = "left")
  })
}
dev.off()

#===============================================================================
#
# PANEL G
#
#===============================================================================
genes <- data.frame()
for(i in factors){
  genes <- rbind(genes, cbind(name =  names(sort(schpf@feature.loadings[,i], decreasing = T))[1:10],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}

# High AD vs. Low/No AD
#-------------------------------------------------------------------------------
keep_ids <- unique(snuc@meta.data %>% filter(var %in% c("1", "3/4")) %>% pull(projid))

pseudobulk <- AggregateExpression(subset(snuc, subset = projid %in% keep_ids),
                                  features = genes$name,
                                  group.by = "projid",
                                  slot = "count",
                                  assays = "RNA",
                                  return.seurat = T)

pseudobulk@meta.data$projid <- rownames(pseudobulk@meta.data)

meta <- pseudobulk@meta.data %>% 
  left_join(snuc@meta.data %>% distinct(projid, msex, study, var, age_death)) %>%
  mutate(var = factor(ifelse(var == "3/4", "0", "1"), levels = c("0", "1")),
         study = ifelse(study %in% "ROS ", 1, 0)) %>%
  mutate_at(vars(c("msex", "study", "age_death")), ~as.numeric(as.character(.))) %>%
  mutate(age_death = scale(age_death))
rownames(meta) = pseudobulk@meta.data$projid

counts <- as.matrix(pseudobulk@assays$RNA@counts)
colnames(counts) <- meta$projid

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta,  
                              design = ~var)
dds <- DESeq(dds)
res_high <- results(dds, contrast = c("var", "1", "0"), alpha = 0.05) %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") 

# Intermediate AD vs. Low/No AD
#-------------------------------------------------------------------------------
keep_ids <- unique(snuc@meta.data %>% filter(var %in% c("2", "3/4")) %>% pull(projid))

pseudobulk <- AggregateExpression(subset(snuc, subset = projid %in% keep_ids),
                                  features = genes$name,
                                  group.by = "projid",
                                  slot = "count",
                                  assays = "RNA",
                                  return.seurat = T)

pseudobulk@meta.data$projid <- rownames(pseudobulk@meta.data)

meta <- pseudobulk@meta.data %>% 
  left_join(snuc@meta.data %>% distinct(projid, msex, study, var, age_death)) %>%
  mutate(var = factor(ifelse(var == "3/4", "0", "1"), levels = c("0", "1")),
         study = ifelse(study %in% "ROS ", 1, 0)) %>%
  mutate_at(vars(c("msex", "study", "age_death")), ~as.numeric(as.character(.))) %>%
  mutate(age_death = scale(age_death))
rownames(meta) <- pseudobulk@meta.data$projid

counts <- as.matrix(pseudobulk@assays$RNA@counts)
colnames(counts) <- meta$projid

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta,  
                              design = ~var)
dds <- DESeq(dds)
res_med <- results(dds, contrast = c("var", "1", "0"), alpha = 0.05) %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") 

# PLOT
#-------------------------------------------------------------------------------
cols <- c(unname(proj_cols("blue", "pink")), "black", unname(proj_cols_grey("med grey")))
names(cols) <- c("1", "2", "3", "4")

plot_data <- rbind(res_high %>% mutate(var = 1), res_med %>% mutate(var = 2)) %>%
  filter(!is.na(padj)) %>%
  mutate(lab = ifelse(padj < 0.05 & abs(log2FoldChange) > 0.5, genes, NA), 
         col = ifelse(padj < 0.05 & log2FoldChange < -0.5, "1", NA),
         col = ifelse(padj < 0.05 & log2FoldChange > 0.5, "2", col),
         col = ifelse(padj < 0.05 & abs(log2FoldChange) < 0.5, "3", col),
         col = ifelse(padj > 0.05, "4", col),
         var = factor(var, levels = c(1, 2), labels = c("High (1)", "Intermediate (2)")))

p <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = col), alpha = 0.7)+
  geom_hline(yintercept = -log10(0.05), color = proj_cols_grey("dark grey"), lty = 2, linewidth = 0.5)+
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.5)+
  geom_vline(xintercept = -0.5, color = proj_cols_grey("dark grey"), lty = 2, linewidth = 0.5)+
  geom_vline(xintercept = 0.5, color = proj_cols_grey("dark grey"), lty = 2, linewidth = 0.5)+
  geom_text_repel(aes(label = lab, color = col), show.legend = F, max.overlaps = 20) +
  scale_color_manual(values = cols,
                     labels = c("Downregulated", "Upregulated", "NS (FC < 0.5)", "NS (FDR *p* > 0.05)"))+
  guides(color = guide_legend(override.aes = aes(size = 4)))+
  facet_wrap(~var) +
  labs(x = "log<sub>2</sub>FC (compared to low-likelihood/no AD)", y = "-log<sub>10</sub>(FDR *p*-value)")+
  theme_proj()+
  theme(axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
        panel.grid.major = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(b = -0.1, unit= "in"),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.key.spacing.x = unit(0.05, "in"),
        axis.title.y = element_markdown(size = 14),
        axis.title.x = element_markdown(size = 14),
        strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in")))

pdf(file.path(home_path, "proj_snuc", "analysis", paste0(fig, "_panelH_volcano2.pdf")), height = 4, width = 6)
print(p)
dev.off()
