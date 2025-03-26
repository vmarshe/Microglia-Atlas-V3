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
library(ggtext)
library(igraph)
library(logger)
library(metafor)
library(miloR)
library(patchwork)
library(scater)
library(scCustomize)
library(schex)
library(scran)
library(Seurat)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "/mnt/vast/hpc/homes/vsm2116/projects/microglia_atlas"
res_path <- file.path(home_path, "miloR_AGE")
fig <- "fig_supp05_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_miloR.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))

data@meta.data <- data@meta.data %>%
  mutate(AGE_CAT = ifelse(AGE <=45, "<=45", NA),
         AGE_CAT = ifelse(AGE >45 & AGE <=65, "46-65", AGE_CAT),
         AGE_CAT = ifelse(AGE >65 & AGE <=75, "66-75", AGE_CAT),
         AGE_CAT = ifelse(AGE >75 & AGE <=85, "76-85", AGE_CAT),
         AGE_CAT = ifelse(AGE >85 & AGE <=95, "86-95", AGE_CAT),
         AGE_CAT = ifelse(AGE >95, ">95", AGE_CAT),
         AGE_CAT = factor(AGE_CAT, levels = c("<=45", "46-65", "66-75", "76-85", "86-95", ">95")))

data@meta.data$AGE_75 <- ifelse(data@meta.data$AGE <=75, 0, 1)
data@meta.data$AGE_85 <- ifelse(data@meta.data$AGE <=85, 0, 1)
data@meta.data$AGE_95 <- ifelse(data@meta.data$AGE <=95, 0, 1)

data@meta.data$SEX <- ifelse(data@meta.data$SEX == "M", 1, 0)

#===============================================================================
#
# PANEL A
#
#===============================================================================
dx_cols <- c(proj_cols("purple", "teal"),
             proj_cols_grey("dark grey"),
             proj_cols("yellow", "light blue", "red", "dark blue", "orange1"),
             proj_cols_grey("light grey"),
             proj_cols("blue", "pink", "orange2"))
names(dx_cols) <- c("AD", "ALS/FTD", "CNTRL","DNET", "ET", "HD", "MCI",
                    "MS", "Other","PD", "Stroke", "TLE")

p1 <- data@meta.data %>%
  group_by(AGE_CAT) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = AGE_CAT, y = n)) +
  geom_col(fill = proj_cols_grey("light grey"))+
  geom_text(aes(label = formatC(n, big.mark = ",")), vjust = 0) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 170000),
                     breaks = c(0, 50000, 100000, 150000),
                     label = function(x) format(x, format = "g", big.mark = ","))+
  labs(y = "N microglia")+
  theme_proj()+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, vjust = 1),
        plot.margin = margin(t = 0.1, b = 0.05, unit = "in"))

p2 <- data@meta.data %>%
  group_by(AGE_CAT, DX_CAT) %>%
  summarise(n = n()) %>%
  group_by(AGE_CAT) %>%
  summarise(DX_CAT, n, pct = n/sum(n)) %>%
  ggplot(aes(x = AGE_CAT, y = pct, fill = DX_CAT)) +
  geom_col()+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     label = function(x) formatC(x, format = "g"))+
  scale_fill_manual(name = "Diagnosis", values = dx_cols) +
  labs(y = "Proportion of cells\n by diagnostic category", x = "Age (years)")+
  theme_proj()+
  theme(panel.grid.major = element_blank(),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.key.size = unit(0.2, "in"),
        legend.margin = margin(l = -0.1, unit = "in"))

p <- wrap_plots(p1, p2) + plot_layout(nrow = 2, heights = c(0.3, 0.7))

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 5, width = 6)
print(p)
dev.off()

#===============================================================================
#
# PANELS B-E
#
#===============================================================================


trait <- "AGE_85"
k <- 30
deltaFC <- 1

# PATHS
#-------------------------------------------------------------------------------
res_path <- file.path(home_path, "miloR_AGE", paste0(trait, "_k", k, "_deltaFC", deltaFC))
if(!dir.exists(res_path)) dir.create(res_path)
if(file.exists(file.path(res_path, "log.txt"))) {file.remove(file.path(res_path, "log.txt"))}
log_appender(appender_tee(file.path(res_path, "log.txt")))

# DATA
#-------------------------------------------------------------------------------
log_info("Subsetting cells...")
milo <- data

milo@meta.data <- milo@meta.data %>% rename_at(all_of(vars(trait)), ~"trait")

set.seed(1234)
keep_cells <- data@meta.data %>% 
  filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H")) %>%
  dplyr::select(cell_names, Donor.ID, REGIONS, trait, DX_CAT, AGE) %>% 
  drop_na() %>%
  group_by(Donor.ID, REGIONS) %>%
  slice_sample(n = 500) %>%
  pull(cell_names)

milo <- subset(milo, subset = cell_names %in% keep_cells)
milo@meta.data$new_id <- paste0(gsub(" ", "", milo@meta.data$Donor.ID), "_", milo@meta.data$REGIONS)
meta <- milo@meta.data[,c("new_id", "SEX", "REGIONS", "ORIGIN", "trait")]

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
milo <- makeNhoods(milo, prop = 0.1, k = k, d = 23,
                   refined = TRUE, reduced_dim = "SCHPF", refinement_scheme = "graph")

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

#===============================================================================
#
# PANEL M - Density Plots
#
#===============================================================================
res_path <- file.path(home_path, "miloR_AGE")

xmin <- min(data@reductions$schpf.umap@cell.embeddings[,1])
xmax <- max(data@reductions$schpf.umap@cell.embeddings[,1])
ymin <- min(data@reductions$schpf.umap@cell.embeddings[,2])
ymax <- max(data@reductions$schpf.umap@cell.embeddings[,2])

keep_cells <- data@meta.data %>%
  filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H")) %>%
  pull(cell_names)

data <- subset(data, subset = cell_names %in% keep_cells)

data@meta.data <- data@meta.data %>%
  mutate(AGE =  as.numeric(gsub("\\+", "", as.character(AGE))),
         AGE_CAT = ifelse(AGE <=75, "<=75", NA),
         AGE_CAT = ifelse(AGE >75 & AGE <=85, "76-85", AGE_CAT),
         AGE_CAT = ifelse(AGE >85 & AGE <=95, "86-95", AGE_CAT),
         AGE_CAT = ifelse(AGE >95, ">95", AGE_CAT),
         AGE_CAT = factor(AGE_CAT, levels = c("<=75", "76-85", "86-95", ">95")))

age75 <- subset(data, subset = AGE_CAT %in% "<=75")
age85 <- subset(data, subset = AGE_CAT %in% "76-85")
age95 <- subset(data, subset = AGE_CAT %in% "86-95")
age95_plus <- subset(data, subset = AGE_CAT %in% ">95")

milo <- read_rds(file.path(res_path, "AGE_85_k30_deltaFC1/milo.rds"))
milo_res <- read_rds(file.path(res_path, "AGE_85_k30_deltaFC1/milo_grouped_Nhoods.rds"))
comparison <- read_rds(file.path(res_path, "AGE_85_k30_deltaFC1/comparison_test.rds")) %>% 
  mutate(dir = ifelse(t > 0, "up", "down"))

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
    scale_size(range = c(0.2, 2), name="Nhood size") +
    theme_umap()+ 
    theme(legend.position = "none")  
  
  nhood_plots <- c(nhood_plots, list(p))
}


age_cols <- proj_cols("teal", "yellow", "orange2", "purple")
names(age_cols) <- c("age75", "age85", "age95", "age95_plus")

age_annot <- c("<=75", "76-85", "86-95", ">95")
names(age_annot) <- c("age75", "age85", "age95", "age95_plus")

plots <- list()
for(i in c(1, 9, 17, 24)){
  
  ind <- which(factors %in% paste0("scHPF_", i))
  
  tmp <- list()
  for(j in c("age75", "age85", "age95", "age95_plus")){
    
    p <- Plot_Density_Custom(eval(parse(text = j)), features = paste0("scHPF_", i), reduction = "schpf.umap")
    p$data$facet = age_annot[j]
    p <- ggplot()+
      geom_point_rast(data = p$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5)+
      scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
      scale_y_continuous(name = names(factors)[ind], limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
      scale_color_gradientn(colors = viridis_magma_light_high)+
      theme_umap()+
      theme(strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in"), family = "Helvetica"),
            strip.background = element_rect(fill = lighten(age_cols[j], 0.6, space = "HLS") , color = age_cols[j], linewidth = 1.5),
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
    
    if(length(plots) < 1){
      p <- p +facet_wrap(~facet, strip.position = "top")
    }
    
    if(j == "age75"){
      p <- p + theme(axis.title.y = element_textbox_simple(family = "Helvetica", face = "bold", size = 10, orientation = "left-rotated",
                                                          hjust = 1, halign = 0.5, box.color = "#F0F0F0", fill = "#F0F0F0", padding = margin(2.5, 2.5, 2.5, 2.5),
                                                          margin = margin(0, 0, -0.05, 0, unit = "in")))
    } else {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    tmp <- c(tmp, list(p))
    
  }
  
  plots <- c(plots, list(wrap_plots(tmp) + plot_layout(nrow = 1)))
}
cairo_pdf(file.path(res_path, paste0(fig, "panel_F.pdf")), height = 9, width = 13)
wrap_plots(nhood_plots[[1]], plots[[1]], 
           nhood_plots[[2]], plots[[2]],
           nhood_plots[[3]], plots[[3]],
           nhood_plots[[4]], plots[[4]]) + 
  plot_layout(widths = c(0.175, 0.825), nrow = 4, ncol = 2)
dev.off()

#===============================================================================
#
# PANEL G
#
#===============================================================================
aggregated_scores <- data@meta.data %>%
  filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H"))  %>%
  dplyr::select(Donor.ID, REGIONS, SEX, ORIGIN, scHPF_1:scHPF_26) %>%
  group_by(Donor.ID, REGIONS) %>%
  summarise_at(.vars = vars(all_of(factors)), .funs = function(x) mean(x))

set.seed(1234)
meta <- data@meta.data %>% 
  filter(DX_CAT == "AD" & REGIONS %in% c("AWS", "BA20_21", "BA9_46", "H"))  %>%
  dplyr::select(Donor.ID, AGE_85, REGIONS, SEX, ORIGIN) %>%
  group_by(Donor.ID, REGIONS) %>% 
  sample_n(1) %>%
  ungroup() %>%
  mutate(lab = paste0(REGIONS, "_", AGE_85))

aggregated_scores <- aggregated_scores %>%
  left_join(meta, by = c("Donor.ID", "REGIONS"))

res <- data.frame()
for(i in factors){
  for(j in c("AWS", "BA9_46", "BA20_21", "H")){
    
    tmp <- aggregated_scores %>%
      rename_at(all_of(i), ~"factor") %>%
      filter(REGIONS == j)
    
    out <- lm(factor ~ AGE_85 + SEX, data = tmp)
    coef <- coef(summary(out))
    confint <- confint(out)
    res <- rbind(res, cbind(factor = i,
                            region = j,
                            n = nrow(tmp), 
                            coef = coef["AGE_85", "Estimate"], 
                            se = coef["AGE_85", "Std. Error"], 
                            lowerCI = confint["AGE_85", "2.5 %"],
                            upperCI = confint["AGE_85", "97.5 %"],
                            p = coef["AGE_85", "Pr(>|t|)"]))
  }
}
res <- res %>%
  mutate_at(all_of(c("n","coef", "se", "p")), ~as.numeric(as.character(.))) %>%
  group_by(region) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

meta.res <- data.frame()
for(i in factors){
  
  tmp <- res %>% filter(factor == i)
  meta <- rma(yi = tmp$coef, sei = tmp$se, weights = tmp$n, data = tmp, method = "REML")
  
  meta.res <- rbind(meta.res, cbind(factor = i, 
                                    meta.beta = c(meta$beta), 
                                    meta.se = meta$se, 
                                    meta.ci.lb = c(meta$ci.lb),
                                    meta.ci.ub = c(meta$ci.ub),
                                    meta.z = meta$zval,
                                    meta.p = c(meta$pval),
                                    meta.tau = meta$tau2,
                                    meta.se.tau = meta$se.tau2,
                                    meta.QE = meta$QE,
                                    meta.QEp = meta$QEp))
  
}
meta.res$p.adj <- p.adjust(meta.res$meta.p, method = "BH")

mat <- res %>%
  mutate(log10padj = ifelse(coef > 0, -log10(p.adj), log10(p.adj))) %>%
  dplyr::select(region, factor, log10padj) %>%
  spread(factor, log10padj) %>%
  column_to_rownames("region") %>%
  as.matrix()

meta.mat <- meta.res %>%
  mutate(log10padj = ifelse(meta.beta > 0, -log10(p.adj), log10(p.adj))) %>%
  dplyr::select(factor, log10padj) %>%
  spread(factor, log10padj) %>%
  as.matrix()
rownames(meta.mat) <- "meta"

mat <- rbind(mat, meta.mat) 
mat <- mat[, factors] 
colnames(mat) <- names(factors)
rownames(mat) <- c("AWS", "BA20/21", "BA9/46", "H", "Meta-analysis")

row_split <- c(0, 0, 0, 0, 1)
row_split <- factor(row_split, levels = c("0", "1"))

cairo_pdf(file.path(res_path, paste0(fig, "panel_G.pdf")), height = 3.5, width = 9.5)
set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(min(mat, na.rm=T), 0, max(mat, na.rm=T)),
                              colors = c(unname(proj_cols("blue")),"white", unname(proj_cols("red")))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & abs(mat[i, j]) > -log10(0.05)) {
                 grid.points(pch = 8, x, y, size = unit(0.065, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             row_names_side = "left",
             row_names_gp = gpar(fontfamily = "Helvetica", fontsize = 12),
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 0),
             row_gap = unit(0.1, "in"),
             right_annotation = rowAnnotation(foo = anno_empty(border = FALSE, width = unit(0.3, "in"))),
             row_split = row_split,
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 0),
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(FDR '*italic(p)*')'),
               legend_height = unit(1, "in"),
               border = "black",
               title_gp = gpar(fontfamily = "Helvetica", fontsize = 12),
               labels_gp = gpar(fontfamily = "Helvetica", fontsize = 10)
             ))

draw(p,
     align_heatmap_legend = "heatmap_center",
     background = "transparent", 
     padding = unit(c(0.2, 0.2, 0.2, 0.2), "in"))
for(i in 1:2){
  decorate_heatmap_body("mat", slice = i, { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2)) })
  decorate_annotation("foo", slice = i, { grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left") })
}

dev.off()

write_xlsx(res, file.path(res_path, "age_assoc_lm.xlsx"))
write_xlsx(meta.res, file.path(res_path, "age_assoc_meta.xlsx"))
