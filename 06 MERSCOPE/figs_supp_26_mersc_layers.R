#!/usr/bin/Rscript

#===============================================================================
# 
# MICROGLIA ATLAS V3
# Figure S26.	Characterization of tissue layer niches and microglia factors. 
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib")

library(circlize)
library(ComplexHeatmap)
library(data.table)
library(ggsignif)
library(ggtext)
library(grid)
library(Hmisc)
library(lmerTest)
library(patchwork)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(viridisLite)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope")
fig <- "fig_s26"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_seurat.R"))
source(file.path(home_path, "src/00 Helpers/funcs_viz.R"))

# DATA
#-------------------------------------------------------------------------------
run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

keep_cells <- rownames(data@meta.data %>% filter(!is.na(cell_type)))
data <- subset(data, subset = cell_names %in% keep_cells)
data <- NormalizeData(data, assay = "RNA") %>% ScaleData(assay = "RNA")

keep_cells <- rownames(data@meta.data %>% filter(cell_type == "Microglia"))
mglia <- subset(data, subset = cell_names %in% keep_cells)
mglia <- NormalizeData(mglia, assay = "RNA") %>% ScaleData(assay = "RNA")

mglia@reductions[["umap"]] = readRDS(file.path(home_path, "merscope/mglia/mglia.umap.rds"))
mglia@reductions[["schpf.umap"]] = readRDS(file.path(home_path, "merscope/mglia/mglia.schpf.umap.rds"))

data@images[["merged"]] <- readRDS(file.path(home_path, "merscope", "merged", "integrated_coords.rds"))

aracne <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
enrichment <- read_rds(file.path(home_path, "aracne_sc/TF_scHPF_enrichment.rds"))

#===============================================================================
#
# PANEL A - DISTRIBUTION OF MICROGLIA ACROSS LAYERS
#
#===============================================================================

# PROPORTION OF MICROGLIA
#-------------------------------------------------------------------------------
mglia@meta.data$wm <- ifelse(mglia@meta.data$layer_niches == "WM", "WM", "GM")

p1 <- mglia@meta.data %>%
  group_by(batch_name) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = batch_name, y = n)) +
  geom_col(fill = proj_cols_grey("light grey"))+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 9000),
                     breaks = c(0, 4000, 8000),
                     label = function(x) format(x, format = "g", big.mark = ","))+
  labs(y = "N microglia")+
  theme_proj()+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, vjust = 1),
        plot.margin = margin(t = 0.1, b = 0.05, unit = "in"))

layer_cols <- unname(proj_cols("yellow","teal", 'pink'))
names(layer_cols) <- sort(unique(mglia@meta.data$layer_niches))
layer_labels <- c("GM (LI)", "GM (LII-VI)", "WM")
names(layer_labels) <- sort(unique(mglia@meta.data$layer_niches))

p2 <- mglia@meta.data %>%
  group_by(batch_name, layer_niches) %>%
  summarise(n = n()) %>%
  group_by(batch_name) %>%
  summarise(layer_niches, n, pct = n/sum(n)) %>%
  ggplot(aes(x = batch_name, y = pct, fill = layer_niches)) +
  geom_col()+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     label = function(x) formatC(x, format = "g"))+
  scale_fill_manual(name = "Layer", values = layer_cols, labels = layer_labels)+
  guides(fill = guide_legend(nrow = 1))+
  labs(y = "Proportion of cells", x = "Tissue")+
  theme_proj()+
  theme(panel.grid.major = element_blank(),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.key.size = unit(0.2, "in"),
        legend.title.position = "left",
        legend.direction = "horizontal",
        legend.margin = margin(b = -0.1, unit = "in"))

p <- wrap_plots(p1, p2) + 
  plot_layout(nrow = 2, heights = c(0.3, 0.7), guides = "collect") & 
  theme(legend.position = "bottom")

pdf(file.path(res_path, paste0(fig, "_panelA_mglia_prop.pdf")), height = 6, width = 5)
print(p)
dev.off()

#===============================================================================
#
# PANEL B - ARE THERE DIFFERENCES BETWEEN WHITE AND GRAY MATTER?
#
#===============================================================================
genes <- c("IRF7", names(sort(schpf@feature.loadings[,"scHPF_20"], decreasing = T))[1:20], 
           "ARID5B", "CEBPA", "MITF", "PPARG",
           names(sort(schpf@feature.loadings[,"scHPF_26"], decreasing = T))[1:20])
genes <- unique(genes[genes %in% rownames(data)])

mglia@meta.data <- mglia@meta.data %>%
  mutate(nCount_RNA_scaled = scale(nCount_RNA),
         dx = ifelse(donor_name %in% c("AD.26", "AD.33"), 1, 0),
         layer = recode(layer_niches, "GM.SLC17A7.GAS7.HSP90AB1"="GM", "GM.AHNAK.FOS.ACTB"="GM"),
         layer_lab = paste0(donor_name, ".", layer))

Idents(mglia) <- "layer"

res <- FindMarkers(mglia, 
                   ident.1 = "WM", 
                   ident.2 = "GM", 
                   features = genes,
                   test.use = "MAST", 
                   assay = "RNA",
                   slot = "data",
                   latent.vars = c("nCount_RNA_scaled", "orig.ident", "dx"), 
                   re.var ="orig.ident",
                   ebayes = F,
                   max.cells.per.ident = 2500,
                   min.cells.group = 50,
                   random.seed = 1234, 
                   min.pct = 0.1,
                   logfc.threshold = 0.25, 
                   verbose = F) %>% 
  filter(p_val_adj < 0.05) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

pdf(file.path(res_path, paste0(fig, "_panelB.pdf")), height = 5.5, width = 4.5, onefile = T)
for(i in paste0("scHPF_", c(20, 26))){
  
  genes <- if(i == "scHPF_20"){c("IRF7")} else {c("ARID5B", "CEBPA", "MITF", "PPARG")}
  genes <- c(genes, names(sort(schpf@feature.loadings[,i], decreasing = T))[1:20])
  genes <- unique(genes[genes %in% rownames(data)])
  
  set.seed(1234)
  Clustered_DotPlot_Single_Group_v2(mglia, 
                                    identity_lab = "", 
                                    assay = "RNA", 
                                    colors_use_exp = viridis_plasma_light_high, 
                                    group.by = "layer_lab", 
                                    plot_padding = unit(c(0.2, 0.3, 0.2, 0.2), "in"),
                                    features = genes,
                                    column_label_size = 9, 
                                    row_label_size = 9,
                                    max_col_size = 2.5,
                                    x_lab_rotate = T)
  
}
dev.off()

#===============================================================================
#
# PANEL C - FACTOR EXPRESSION ACROSS WM and GM
#
#===============================================================================
mglia@meta.data <- mglia@meta.data %>% mutate(layer = ifelse(layer_niches == "WM", "WM", "GM"))

layer_cols <- proj_cols("pink", "teal")
names(layer_cols) <- c("WM", "GM")

plots <- list()
for(i in paste0("scHPF_", c(20, 26))){
  
  tmp <- mglia@meta.data %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(DX = ifelse(donor %in% c(26, 33), 1, 0))
  
  out <- lmer(factor ~ scale(nCount_RNA) + layer + DX + (1|batch_name), data = tmp, na.action = na.omit)
  out <- coefficients(summary(out))
  anno.p <- ifelse(out[3,5] < 0.05, "*p* < 0.001", paste0("*p* = ", formatC(out[3,5], format = "e", digits = 2)))
  anno.beta <- paste0("*β* (SE) = ", formatC(out[3,1], digits = 2), " (", formatC(out[3,2], digits = 2), ")")
  anno <- data.frame(label = paste0(anno.beta, ", ", anno.p))
  
  p <-  tmp %>% 
    mutate(lab = gsub('-high', paste0("<sup>high</sup>"), names(factors)[factors %in% i])) %>%
    ggplot(aes(x = donor_name, y = factor, fill = layer)) + 
    geom_boxplot(outlier.alpha = 0.7)+
    geom_richtext(data = anno, aes(label = label, x = Inf, y = Inf), hjust = 1, 
                  vjust = 1, inherit.aes = F, fill = alpha("#F0F0F0", 0.5),
                  label.color = NA, size = 2.5,
                  label.margin = unit(c(t=0.02, 0, 0, 0), "in"),
                  label.padding = unit(c(0.01, 0.01, 0.01, 0.01), "in"))+
    labs(y = "Factor expression level")+
    facet_wrap(~lab)+
    scale_fill_manual(values = layer_cols)+
    theme_proj()+
    theme(panel.grid.major = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          strip.text = element_markdown(size = 10, margin = margin(b = 0.05, t= 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + 
  plot_layout(nrow = 1, guides = "collect", axis_titles = "collect", axes = "collect_x") &
  theme(legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.margin = margin(b = -0.1, unit = "in"), 
        legend.position = "top",
        legend.direction = "horizontal",
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in"))

cairo_pdf(file.path(res_path, paste0(fig, "_panelC.pdf")), height = 3, width = 5)
print(p)
dev.off()

#===============================================================================
#
# PANEL D - IN SITU DISTRIBUTION OF LOW/HIGH EXPRESSING MICROGLIA
#
#===============================================================================
cols <- proj_cols("pink", 'yellow', "teal")
names(cols) <- c("WM", "GM\n(L1)", "GM\n(L2-6)")

plots <- list()
for(j in c(20, 26)){
  tmp <-  mglia@meta.data %>%
    left_join(niches %>% dplyr::select(cell_names, orig.ident, layer_niches)) %>%
    rename_at(all_of(paste0("scHPF_", j, "_class")), ~"class") %>%
    mutate(class = str_extract(class, "high|low"),
           class = ifelse(class == "high", 1, 0),
           layer_niches = factor(layer_niches, levels = c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1"),
                                 labels =  c("WM", "GM\n(L1)", "GM\n(L2-6)"))) %>%
    group_by(orig.ident, layer_niches) %>%
    summarise(pct = mean(class), nCount_RNA = sum(nCount_RNA)) %>%
    mutate(DX = ifelse(orig.ident %in% ( run_metadata %>% filter(donor %in% c(1, 9)) %>% pull(batch_code)), 0, 1))
  
  out <- lmer(pct ~ layer_niches + scale(nCount_RNA) + DX + (1|orig.ident), data = tmp)
  out <- anova(out)["layer_niches","Pr(>F)"]
  out <- ifelse(out < 0.001, "***", ifelse(out < 0.01, "**", ifelse(out < 0.05, "*", "")))
  
  anno <- TukeyHSD(aov(pct~layer_niches, tmp), conf.level=0.95)$layer_niches
  anno <- unname(ifelse(anno[,4] < 0.05, formatC(unname(anno[,4]), digits = 1), "ns"))
  
  p <- tmp %>%
    mutate(facet = paste0(names(factors)[which(factors %in% paste0("scHPF_", j))], out), 
           lab = ifelse(pct*100 > 1, paste0(round(pct*100, 0), "%"), paste0(round(pct*100, 1), "%"))) %>%
    ggplot(aes(x = layer_niches, y = pct, fill = layer_niches)) +
    geom_boxplot() +
    geom_point()+
    geom_signif(comparisons = list(c("WM", "GM\n(L1)"), c("WM", "GM\n(L2-6)"), c("GM\n(L1)", "GM\n(L2-6)")),
                annotation = anno,
                step_increase = 0.1, 
                tip_length = 0)+
    scale_fill_manual(values = cols) +
    facet_wrap(~facet)+
    theme_proj()+
    theme(axis.text.x = element_text(size = 12), 
          panel.grid.major = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(size = 14),
          axis.text.y = element_text(size = 12),
          legend.position = "none",
          strip.text = element_text(size = 10, margin = margin(t = 0.05, b = 0.05, unit = "in")))
  
  p <- p + 
    scale_y_continuous(name = "Proportion of<br>high-expressing microglia",
                       limits = c(0, layer_scales(p)$y$range$range[2] + 0.01))
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + 
  plot_layout(ncol = 2, axis_titles = "collect", axes = "collect") &
  theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"))

pdf(file.path(res_path, paste0(fig, "_panelD.pdf")), height = 3, width = 5)
print(p)
dev.off() 
