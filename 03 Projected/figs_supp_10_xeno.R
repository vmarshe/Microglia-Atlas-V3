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
library(ggrepel)
library(ggridges)
library(ggsignif)
library(ggtext)
library(grid)
library(igraph)
library(lmerTest)
library(logger)
library(miloR)
library(patchwork)
library(scater)
library(scCustomize)
library(schex)
library(scran)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
dataset <- "xeno"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, paste0("proj_", dataset), "analysis")
fig <- "fig_supp10_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/funcs_helpers.R"))
source(file.path(home_path, "src/00 Helpers/funcs_miloR.R"))

# REFERENCE DATA
#-------------------------------------------------------------------------------
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

# DATA
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  as.data.frame() %>% 
  mutate(treatment = ifelse(grepl("WT", cell_names), 0, 1),
         treatment_name = ifelse(grepl("WT", cell_names), "WT", "5X"),
         sex = ifelse(grepl("Female", cell_names), 1, 0)) %>%
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern ="^MT-|^MTRN")

id_labs <- factor(c("WTFemale", "WTMale", "5XFemale", "5XMale"))
names(id_labs) <- c("WT-F", "WT-M", "5X-F", "5X-M")

# HEX
#-------------------------------------------------------------------------------
hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")

p1 <- ggplot() + 
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N Cells")+
  labs(title = "Murine-xenografted<br>human microglia (scRNA-seq)")+
  theme_umap()+
  theme(legend.key.height = unit(0.3, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 12, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
        legend.margin = margin(0, 0, 0, -0.1, 'in'),
        plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 16, family = "Helvetica"),
        plot.margin = margin(0.05, 0.2, 0.05, 0.05))

pdf(file.path(res_path, paste0(fig, "panel_A1.pdf")), height = 3.5, width = 3.5)
print(p1)
dev.off()

plot_list <- list()
for(i in c("WT", "5X")){
  hex <- make_hexbin(subset(data, subset = treatment_name %in% i), nbins = 50, dimension_reduction = "schpf.umap")
  p <- ggplot() + 
    geom_point(data = ref_hex, aes(x = x, y = y), color = "#F0F0F0") +
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    ggtitle(paste0(i, " (", format(nrow(hex@meta.data), big.mark = ","), ")"))+
    theme_umap()+
    theme(legend.key.height = unit(0.2, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
          legend.title = element_blank(),
          legend.margin = margin(t = 0, r = 0, b = 0, l = -0.2, 'in'),
          plot.title = element_text(margin = margin(0, 0, 0.01, 0, 'in'), size = 14, family = "Helvetica"),
          panel.border = element_rect(linewidth = 1.5, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          plot.margin = margin(0.05, 0.05, 0.05, 0.05))
  plot_list = c(plot_list, list(p))
}

p <- wrap_plots(c(list(p1), plot_list)) + 
  plot_layout(
    design = 
      "111122
       111133")

pdf(file.path(res_path, paste0(fig, "panel_A2.pdf")), height = 5, width = 7)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
id_cols <- c(lighten(proj_cols("blue"), 0.5, space = "HLS"),
             lighten(proj_cols("yellow"), 0.5, space = "HLS"),
             proj_cols("blue"), proj_cols("yellow"))
names(id_cols) <- c("WT-F", "WT-M", "5X-F", "5X-M")

plots <- list()
for(i in factors){
  
  tmp <- data@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(facet = names(factors)[factors %in% i], 
           treatment_name = ifelse(grepl("WT", cell_names), "WT", "5X"),
           name = paste0(treatment_name, "-", ifelse(grepl("Female", cell_names), "F", "M")),
           name = factor(name, levels = rev(names(id_cols))))
  
  p <- ggplot(tmp, aes(x = factor, y = name, fill = name)) + 
    geom_density_ridges()+
    scale_x_continuous(trans = "log2", labels = function(x) formatC(x, digits = 1,  format = "g")) + 
    geom_vline(xintercept = log2(0.01), lty = 2) + 
    scale_fill_manual(values = id_cols)+
    facet_wrap(~facet)+
    labs(x = "log(Cell score)")+
    theme_proj() +
    theme(strip.text = element_text(size = 10, margin = margin(b = 0.05, t = 0.05, unit = "in"), family = "Helvetica"),
          axis.text = element_text(size = 10),
          axis.title.x =  element_blank(),
          axis.title.y =  element_blank(),
          legend.position = "none")
  
  plots <- c(plots, list(p))
}

plots[[21]] <- plots[[21]] + theme(axis.title.x = element_text(size = 14, family = "Helvetica"))
p <- wrap_plots(plots) + plot_layout(guides = "collect", nrow = 4)

cairo_pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 7, width = 14)
print(p)
dev.off()
#===============================================================================
#
# PANEL C
#
#===============================================================================
prop <- data.frame()
for(i in factors){
  
  tmp <- data@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

keep_factors <- prop %>% filter(as.numeric(as.character(prop)) >= 0.25) %>% pull(i)

p <- prop %>%
  mutate(prop = as.numeric(as.character(prop)), 
         factor = factor(i, levels = factors, labels = names(factors)),
         lab = paste0(round(as.numeric(prop)*100, 1), "%")) %>%
  arrange(factor) %>%
  ggplot(aes(x = factor, y = prop)) +
  geom_col(aes(fill = ifelse(prop > 0.25, "1", "0"))) +
  geom_text(aes(label = lab, hjust = ifelse(prop > 0.28, 1, 0)), vjust = 0.5, family = "Helvetica", size = 3, angle = 90)+
  geom_hline(yintercept = 0, linewidth = 1)+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1.1), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))+
  scale_fill_manual(values = c(`0`="#969696", `1`="#ffb200"))+
  labs(y = "Proportion of\nexpressing cells")+
  theme_proj()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        plot.margin = margin(t=0.1, r=0.1, b=0.1, l=0.7, unit = "in"),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 11, margin = margin(-0.3, 0, 0, 0, unit = "in"), family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 11, family = "Helvetica", angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5))

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 3, width = 8)
print(p)
dev.off()
#===============================================================================
#
# PANEL D
#
#===============================================================================
res <- data.frame()
for(i in keep_factors){
  
  tmp <- data@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(condition = as.numeric(as.character(treatment)),
           scale_UMI = c(scale(nCount_RNA)), 
           scale_MT = c(scale(percent.mt))) %>%
    dplyr::select(orig.ident, scale_UMI, scale_MT, condition, factor) %>%
    group_by(orig.ident)
  
  out <- lmer(log2(factor) ~ scale_UMI + scale_MT + condition + (1 | orig.ident),
              data = tmp, 
              na.action = na.omit)
  
  coef <- coef(summary(out))
  
  log2FC <- logfc(tmp$factor[tmp$condition==1], tmp$factor[tmp$condition==0])
  
  rand_eff <- VarCorr(out)$orig.ident[1]
  sing <- any(grepl("singular", unlist(out@optinfo$conv)))
  
  if(rand_eff == 0 | sing == T){
    
    d <- datadist(tmp); options(datadist = "d")
    out <- ols(log2(factor) ~ scale_UMI + scale_MT + condition, data = tmp)
    coef <- summary.lm(out)$coefficients
    
  }
  
  confint <- confint(out, parm = "condition")
  
  res <- rbind(res, cbind(factor = i,
                          n = nrow(tmp), 
                          log2FC,
                          coef = coef["condition", "Estimate"], 
                          se = coef["condition", "Std. Error"], 
                          lowerCI = confint["condition", "2.5 %"],
                          upperCI = confint["condition", "97.5 %"],
                          rand_eff,
                          sing,
                          p = coef["condition", "Pr(>|t|)"]))
  
  rm(list = ls()[ls() %in% c("out", "coef", "tmp", "p")]) 
  
}

res <- res %>%
  mutate_at(all_of(c("log2FC", "coef", "se", "lowerCI", "upperCI","p", "rand_eff")), ~as.numeric(as.character(.))) %>%
  mutate(factor = factor(factor, levels = factors, labels = names(factors)),
         p.adj = p.adjust(p, method = "BH"),
         col = ifelse(coef > 0 & p < 0.05, "Upregulated (Unadj. *p* < 0.05)", "NS (Unadj. *p* > 0.05)"),
         col = ifelse(coef < 0 & p < 0.05, "Downregulated (Unadj. *p* < 0.05)", col),
         col = factor(col, levels = c("Downregulated (Unadj. *p* < 0.05)", "Upregulated (Unadj. *p* < 0.05)", "NS (Unadj. *p* > 0.05)")))

write_xlsx(res, file.path(res_path, paste0(fig, dataset, "_res.xlsx")))

cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("Downregulated (Unadj. *p* < 0.05)", "Upregulated (Unadj. *p* < 0.05)", "NS (Unadj. *p* > 0.05)")

p <- ggplot(res, aes(y = log2FC, x = factor, fill = col)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)+
  geom_col()+
  scale_fill_manual(values = cols)+
  geom_hline(yintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
  geom_hline(yintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
  labs(y = "Log<sub>2</sub>FC<br>(5X vs. WT)")+
  scale_y_continuous(labels = function(x) formatC(x, format = "g"), limits = c(-2, 2), breaks = seq(-2, 2, 1))+
  theme_proj() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_markdown(size = 11, family = "Helvetica"),
        axis.text.x = element_text(size = 11, family = "Helvetica", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11, family = "Helvetica"),
        plot.margin = margin(t=0.1, r=0.1, b=0.1, l=0.3, unit = "in"),
        panel.grid.major.y = element_line(linewidth = 0.25, color = "#F0F0F0"),
        panel.grid.major.x = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5), 
        axis.ticks.x = element_line(linewidth = 0.5), 
        panel.border = element_rect(linewidth = 1, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        legend.margin = margin(l = -0.15, unit = "in"),
        legend.title = element_blank(),
        legend.text = element_markdown(size = 10, family = "Helvetica", margin = margin(l = 0.05, unit = "in")),
        legend.key.size = unit(0.15, "in"),
        legend.spacing.x = unit(0.02, "in"))

cairo_pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 3, width = 10)
print(p)
dev.off()
#===============================================================================
#
# PANEL E-H
#
#===============================================================================

# OPTIONS
#-------------------------------------------------------------------------------
trait <- "treatment"
k <- 30
deltaFC <- 3

# SUBSET
#-------------------------------------------------------------------------------
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  as.data.frame() %>% 
  mutate(treatment = ifelse(grepl("WT", cell_names), 0, 1),
         treatment_name = ifelse(grepl("WT", cell_names), "WT", "5X"),
         sex = ifelse(grepl("Female", cell_names), 1, 0)) %>%
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names

data@meta.data <- data@meta.data %>% rename_at(all_of(vars(trait)), ~"trait")

log_info("Subsetting cells...")
set.seed(1234)
keep_cells <- data@meta.data %>% 
  group_by(orig.ident) %>%
  slice_sample(n = 5000) %>%
  pull(cell_names)

log_info("Keeping ", length(keep_cells), " cells.")

milo <- subset(data, subset = cell_names %in% keep_cells)
milo$new_id <- paste0(milo$orig.ident, "_", milo$trait)
meta <- milo@meta.data[,c("new_id", "trait")]

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
milo <-  buildNhoodGraph(milo, overlap = 1)
saveRDS(milo@nhoods, file.path(res_path, "milo_Nhood_graph.rds"))

log_info("Grouping neighborhoods...")
da_results <- groupNhoods(x = milo, 
                          da.res = res, 
                          da.fdr = 0.05,
                          max.lfc.delta = deltaFC, 
                          overlap = 1, 
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
plot_nhood_dotplot(milo_agg_factors, res_path = res_path, factors = factors, height = 3.5)
write_xlsx(milo_agg_factors, file.path(res_path, "milo_agg_factors.xlsx"))

saveRDS(milo, file.path(res_path, "milo.rds"))
invisible(lapply(capture.output(sessionInfo()), log_info))

#===============================================================================
#
# PANEL I - DENSITY PLOTS
#
#===============================================================================
data <- readRDS(file.path(home_path, paste0("proj_", dataset), "projection", paste0(dataset, ".rds")))

data@meta.data <- data@meta.data %>% 
  as.data.frame() %>% 
  mutate(treatment = ifelse(grepl("WT", cell_names), 0, 1),
         treatment_name = ifelse(grepl("WT", cell_names), "WT", "5X"),
         sex = ifelse(grepl("Female", cell_names), 1, 0)) %>%
  left_join(data@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(data@meta.data) <- data@meta.data$cell_names

hex <- make_hexbin(data, nbins = 100, dimension_reduction = "schpf.umap")
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
  top_n(1, wt = abs(diff)) %>%
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
  
  nhood_plots <- c(nhood_plots, list(p))
}

wt <- subset(data, subset = treatment == 0)
non_wt <- subset(data, subset = treatment == 1)

xmin <- min(ref_hex$x)
xmax <- max(ref_hex$x)
ymin <- min(ref_hex$y)
ymax <- max(ref_hex$y)

treat_cols <- proj_cols("teal", "purple")
names(treat_cols) <- c("WT", "5X")

plots <- list()
for(i in c(20, 26)){
  
  ind <- which(factors %in% paste0("scHPF_", i))
  
  p1 <- Plot_Density_Custom(wt, features = paste0("scHPF_", i), pt.size = 0.5)
  p1$data$facet <- "WT"
  p1 <- ggplot()+
    geom_point_rast(data = p1$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature), size = 0.1, scale = 0.5)+
    facet_wrap(~facet, strip.position = "top")+
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(floor(xmin), ceiling(xmax), 2))+
    scale_y_continuous(name = names(factors)[ind], limits = c(ymin, ymax), breaks = seq(floor(ymin), ceiling(ymax), 2))+
    scale_color_gradientn(colors = viridis_magma_light_high)+
    theme_umap()+
    theme(strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in"), family = "Helvetica"),
          strip.background = element_rect(fill = lighten(treat_cols["WT"], 0.6) , color = treat_cols["WT"], linewidth = 1.5),
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
  
  p2 <- Plot_Density_Custom(non_wt, features = paste0("scHPF_", i), pt.size = 0.5)
  p2$data$facet <- "5X"
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
          strip.background = element_rect(fill = lighten(treat_cols["5X"], 0.6) , color = treat_cols["5X"], linewidth = 1.5),
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
cairo_pdf(file.path(res_path, paste0(fig, "panel_I.pdf")), height = 5, width = 7)
wrap_plots(nhood_plots[[1]], plots[[1]], 
           nhood_plots[[2]], plots[[2]]) + 
  plot_layout(widths = c(0.15, 0.35, 0.15, 0.35), nrow = 2, ncol = 2)
dev.off()

#===============================================================================
#
# PANEL J
#
#===============================================================================
prop <- data.frame()
for(i in factors){
  
  tmp <- data@meta.data %>% rename_at(all_of(i), ~"factor")
  prop_express <- nrow(tmp %>% mutate(class = ifelse(tmp$factor > 0.01, 1, 0)) %>% filter(class == 1))/nrow(tmp)
  prop <- rbind(prop, cbind(i, prop = prop_express))
}

diff_factors <- prop %>% 
  filter(as.numeric(as.character(prop)) >= 0.25) %>% 
  pull(i)

genes <- data.frame()
for(i in diff_factors){
  genes <- rbind(genes, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                              factor = factors[factors %in% i],
                              factor_name = names(factors[factors %in% i])))
}

keep_factors <- paste0("scHPF_", c(20, 11, 26))

keep_highlights <- data.frame()
for(i in keep_factors){
  keep_highlights <- rbind(keep_highlights, cbind(name = names(sort(schpf@feature.loadings[,i], decreasing = T))[1:5],
                                                  factor = factors[factors %in% i],
                                                  factor_name = names(factors[factors %in% i])))
}


id_cols <- c(lighten(proj_cols("blue"), 0.5, space = "HLS"),
             lighten(proj_cols("yellow"), 0.5, space = "HLS"),
             proj_cols("blue"), proj_cols("yellow"))
names(id_cols) <- c("WT-F", "WT-M", "5X-F", "5X-M")

averages <- AverageExpression(data,
                              features = genes$name,
                              group.by = "orig.ident",
                              return.seurat = T,
                              assays = "RNA",
                              slot = "data")$RNA@scale.data
colnames(averages) <- gsub("Female", "-F", colnames(averages))
colnames(averages) <- gsub("Male", "-M", colnames(averages))

averages <- averages[genes$name, names(id_cols)]
averages <- t(averages)

ha <- rowAnnotation(mod = anno_simple(names(id_cols), 
                                      col = id_cols, 
                                      height = unit(0.15, "in"), 
                                      gp = gpar(col = "white", lwd = 1.5)),
                    annotation_name_gp = gpar(fontsize = 0),
                    show_legend = F,
                    border = F)

col_split <- genes %>% filter(name %in% colnames(averages)) %>% pull(factor_name)
col_split <- factor(col_split, levels = unique(genes$factor_name))

Idents(data) <- "treatment_name"
sig <- FindMarkers(data, 
                   ident.1 = "5X", 
                   ident.2 ="WT",
                   features = keep_highlights$name,
                   assay = "RNA",
                   test.use = "MAST", 
                   latent.vars = c("nCount_RNA", "sex"),
                   logfc.threshold = 0,
                   min.pct = 0) %>%
  rownames_to_column("gene") %>%
  mutate(gene = factor(gene, levels = keep_highlights$name),
         cols = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, proj_cols("red"), 
                       ifelse(avg_log2FC < 0 & p_val_adj < 0.05, proj_cols("blue"), "black")))%>%
  rowwise() %>%
  mutate(at = which(colnames(averages) %in% gene)) %>%
  ungroup() %>%
  arrange(at)

bottom_ha <- HeatmapAnnotation(labs = anno_mark(at = sig$at, 
                                                labels = sig$gene,
                                                side ="right",
                                                labels_gp = gpar(fontface = "italic", 
                                                                 fontfamily = "Helvetica", 
                                                                 fontsize = 10,
                                                                 col = sig$cols)),
                               which = "column")

cairo_pdf(file.path(res_path, paste0(fig, "panel_J.pdf")), height = 5, width = 9)
set.seed(1234)
p <- Heatmap(averages, name = "mat",
             col = colorRamp2(breaks = c(min(averages), 0, max(averages)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             left_annotation = ha,
             top_annotation = columnAnnotation(foo = anno_empty(border = FALSE,  height = unit(0.1, "in"))),
             cluster_rows = T,
             row_names_side = "left",
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_dend_side = "right",
             row_dend_width = unit(0.2, "in"),
             column_split = col_split,
             row_title_side = "right",
             row_title_rot = 270,
             row_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 0),
             column_title_rot = 90,
             column_title_gp = gpar(fontface = "bold", fontfamily = "Helvetica", fontsize = 8, hjust = 0),
             column_names_side = "bottom",
             column_names_gp = gpar(fontface = "italic", fontfamily = "Helvetica", fontsize = 0),
             column_gap = unit(0.02, "in"),
             bottom_annotation = bottom_ha,
             cluster_columns = T,
             show_column_dend = FALSE,
             heatmap_legend_param = list(
               title = expression('Average scaled expression level'), 
               legend_width = unit(1, "in"),
               direction = "horizontal",
               title_position = "topcenter",
               border = "black",
               title_gp = gpar(size = 12, fontfamily = "Helvetica"),
               labels_gp = gpar(size = 12, fontfamily = "Helvetica")
             ))


draw(p, 
     align_heatmap_legend = "heatmap_center", 
     heatmap_legend_side = "top", 
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))

for(i in 1:length(diff_factors)){
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, height = unit(0.15, "in"), gp = gpar(fill = "#282828", col = NA), just = "left")
  })
}

dev.off()
