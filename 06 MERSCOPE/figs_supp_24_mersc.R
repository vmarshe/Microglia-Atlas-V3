#!/usr/bin/Rscript

#===============================================================================
# 
# MERSCOPE
# Figure S24. IFN-I Response (scHPF_20) and GPNMB-High (scHPF_26) factors can be identified in situ. 
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(ggrastr)
library(ggrepel)
library(ggtext)
library(grid)
library(logger)
library(parallel)
library(patchwork)
library(scCustomize)
library(schex)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(viridisLite)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope")
fig <- "fig_s24"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

#===============================================================================
#
# DATA
#
#===============================================================================
run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))

# REFERENCE DATA
#-------------------------------------------------------------------------------
ref_hex <- readRDS(file.path(home_path, "~/data/ref_hex.rds"))

# MERSCOPE
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

keep_cells <- rownames(data@meta.data %>% filter(cell_type == "Microglia"))
mglia <- subset(data, subset = cell_names %in% keep_cells)
mglia <- NormalizeData(mglia, assay = "RNA") %>% ScaleData(assay = "RNA")

mglia@reductions[["umap"]] <- readRDS(file.path(home_path, "merscope/mglia/mglia.umap.rds"))
mglia@reductions[["schpf.umap"]] <- readRDS(file.path(home_path, "merscope/mglia/mglia.schpf.umap.rds"))

aracne <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
enrichment <- read_rds(file.path(home_path, "aracne_sc/TF_scHPF_enrichment.rds"))

#===============================================================================
#
# PANEL A - SCHEX
#
#===============================================================================

# HEX
#-------------------------------------------------------------------------------
plots = list()
for(i in 1:nrow(run_metadata)){
  tmp <- readRDS(file.path(home_path, "merscope", run_metadata$run_name[i], "projection/data/data.rds"))
  hex <- make_hexbin(tmp, nbins = 50, dimension_reduction = "schpf.umap")
  
  p <- ggplot() + 
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
    scale_fill_viridis_c(name = "N Cells")+
    labs(title = run_metadata$batch_name[i])+
    theme_umap()+
    theme(legend.key.height = unit(0.15, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
          legend.title = element_text(size = 8, margin = margin(b = 0.1, unit = "in"), family = "Helvetica"),
          legend.margin = margin(l=-0.2, unit='in'),
          plot.margin = margin(0, 0, 0, 0),
          plot.title = element_markdown(margin = margin(0, 0, 0.01, 0, 'in'), size = 12, face = "bold", family = "Helvetica"))
  plots = c(plots, list(p))
  rm(list = c("tmp", "p")); gc()
}

p <- wrap_plots(plots) + plot_layout(nrow = 4)

png(file.path(res_path, paste0(fig, "_panelA_schex_by_tissue.png")), height = 7, width = 6, units = "in", res = 600)
print(p)
dev.off()

#===============================================================================
#
# PANEL B - TF and TARGET GENE EXPRESSION DOTPLOT
#
#===============================================================================

# FACTOR 20 GENES
genes <- names(sort(schpf@feature.loadings[,"scHPF_20"], decreasing = T))[1:20]
targets <- aracne %>% filter(Regulator == "IRF7" & Target %in% genes) %>% pull(Target)
genes <- sort(unique(c(genes, targets, "IRF7")))
genes20 <- genes[genes %in% rownames(data)]

# FACTOR 26 GENES
genes <- names(sort(schpf@feature.loadings[,"scHPF_26"], decreasing = T))[1:20]
targets <- aracne %>% filter(Regulator %in% c("ARID5B", "MITF", "PPARG", "CEBPA") & Target %in% genes) %>% pull(Target)
genes <- sort(unique(c("ARID5B", "MITF", "PPARG", "CEBPA", genes, targets)))
genes26 <- genes[genes %in% rownames(data)]

Idents(mglia) <- "batch_name"

p <- DotPlot(object = mglia, features = c(genes20, genes26), assay = "RNA") +
  geom_point(aes(size = pct.exp),
             shape = 21,
             colour = "black",
             stroke = 0.5) +
  scale_colour_gradientn(colors = viridis_plasma_light_high) + 
  guides(color = guide_colorbar(title = "Avg. Exp", theme = theme(legend.key.height = unit(0.8, "in"), 
                                                                  legend.key.width = unit(0.2, "in"))),
         size = guide_legend(title = "% Exp")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")))

pdf(file.path(res_path, paste0(fig, "_panelB_TF_target_dotplot.pdf")), height = 3, width = 9)
print(p)
dev.off()

#===============================================================================
#
# PANEL C - UMAP OF FACTOR EXPRESSION
#
#===============================================================================
plots <- list()
for(i in c(20, 26)){
  
  p <- FeaturePlot(mglia, features =  paste0("scHPF_",i))
  p$data$feature <- p$data[,4]
  p <- ggplot() +
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_point_rast(data = p$data %>% arrange(feature), aes(x = umap_1, y = umap_2, color = feature), size = 0.1)+
    scale_color_gradient(low = unname(proj_cols_grey("light grey")), high = "blue")+
    labs(title = names(factors)[which(factors %in% paste0("scHPF_", i))])+
    theme_umap()+
    theme(legend.margin = margin(l=-0.2, unit = "in"),
          legend.key.height = unit(0.15, "in"),
          legend.key.width = unit(0.15, "in"),
          legend.text = element_text(size = 10, margin = margin(l=0.05, unit="in")),
          plot.title = element_text(size = 10, margin = margin(b = 0.05, unit="in"), face = "bold"),
          legend.title = element_blank())
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + 
  plot_layout(nrow = 1) & 
  theme(plot.margin = margin(t = 0.05, b = 0.05, unit = "in"))

pdf(file.path(res_path, paste0(fig, "_panelC.pdf")), height = 2.5, width = 5)
print(p)
dev.off()

#===============================================================================
#
# PANEL C - UMAP FACTOR EXPRESSION
#
#===============================================================================
cols <- c(proj_cols_grey("light grey"), proj_cols("pink"))
names(cols) <- c("Low", "High")

plots <- list()
for(i in c(20, 26)){
  
  p <- DimPlot(mglia, group.by = paste0("scHPF_",i, "_class"))
  
  p$data$class <- factor(str_extract(p$data[,3], "high|low"), levels = c("low", "high"), labels = c("Low", "High"))
  
  p <- ggplot() +
    geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
    geom_point_rast(data = p$data %>% arrange(class), aes(x = umap_1, y = umap_2, color = class), size = 0.1)+
    scale_color_manual(values = cols)+
    guides(color = guide_legend(override.aes = aes(size = 4)))+
    labs(title = names(factors)[which(factors %in% paste0("scHPF_", i))])+
    theme_umap()+
    theme(legend.key = element_rect(color = NA),
          legend.margin = margin(l = -0.2, unit="in"),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit="in")),
          plot.title = element_text(size = 10, margin = margin(b = 0.05, unit="in"), face = "bold"),
          legend.title = element_blank())
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + 
  plot_layout(nrow = 1) & 
  theme(plot.margin = margin(t = 0.05, b = 0.05, unit = "in"))

pdf(file.path(res_path, paste0(fig, "_panelD.pdf")), height = 2.5, width = 5)
print(p)
dev.off()

#===============================================================================
#
# PANEL D - LOW/HIGH EXPRESSING MICROGLIA DEGS
#
#===============================================================================
mglia@meta.data$nCount_RNA_scaled = base::scale(mglia@meta.data$nCount_RNA)

all_res <- data.frame()
for(j in c(20, 26)){
  
  Idents(mglia) <- paste0("scHPF_", j, "_class")
  
  for(niche in  c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    if(!(j==26 & niche=="GM.AHNAK.FOS.ACTB")){
      tmp <- subset(mglia, subset = layer_niches %in% niche)
      tmp <- NormalizeData(tmp) %>% ScaleData()
      
      res <- FindMarkers(tmp,
                         ident.1 = paste0("Mic.", j, ".high"),
                         ident.2 = paste0("Mic.", j, ".low"), 
                         assay = "RNA",
                         slot = "data",
                         test.use = "MAST",
                         latent.vars = c("nCount_RNA_scaled", "orig.ident"),
                         re.var = 'orig.ident',
                         ebayes = F,
                         min.pct = 0.1,
                         logfc.threshold = 0.25,
                         max.cells.per.ident = 2500,
                         random.seed = 1234, 
                         verbose = F) %>% 
        rownames_to_column("gene") %>% 
        mutate(factor =  paste0("scHPF_", j), niche)
      all_res = rbind(all_res, res)
    }
  }
}
write_xlsx(all_res, file.path(res_path, paste0(fig, "_tab_microglia_DEGs.xlsx")))


cols <- c(unname(proj_cols("blue", "pink")),  "black", unname(proj_cols_grey("med grey")))
names(cols) <- c("1", "2", "3", "4")
col_labs <- c("Downregulated", "Upregulated", "FC < 0.5", "NS (Adj. p > 0.05)")

layer_labs <- c("WM", "GM (L1)", "GM (L2-6)")
names(layer_labs) <- c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")

plots <- list()
for(j in c(20, 26)){
  
  genes <-  names(sort(schpf@feature.loadings[,paste0("scHPF_", j)], decreasing=T)[1:10])
  
  plot_data <- all_res %>%
    filter(factor == paste0("scHPF_", j)) %>%
    mutate(facet = names(factors)[which(factors %in% paste0("scHPF_", j))]) %>%
    mutate(log10.p.adj = -log10(p_val_adj),
           col = ifelse(avg_log2FC > 0.5 & p_val_adj < 0.05, "2", "4"),
           col = ifelse(avg_log2FC < -0.5 & p_val_adj < 0.05, "1", col),
           col = ifelse(abs(avg_log2FC) < 0.5 & p_val_adj < 0.05, "3", col),
           dir = ifelse(avg_log2FC > 0, "up", "down")) %>%
    group_by(dir) %>%
    arrange(p_val_adj) %>%
    mutate(lab = ifelse(row_number() <=5 & p_val_adj < 0.05 & abs(avg_log2FC)>=0.5, gene, NA),
           lab = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.5 & gene %in% genes, gene, lab),
           niche = recode(niche, !!!layer_labs))
  
  p <- ggplot(plot_data, aes(x = avg_log2FC, y = log10.p.adj, color = col)) +
    geom_hline(yintercept = -log10(0.05), lty = 2, col = proj_cols_grey("light grey"))+
    geom_vline(xintercept = 0.5, lty = 2, col = proj_cols_grey("light grey"))+
    geom_vline(xintercept = -0.5, lty = 2, col = proj_cols_grey("light grey"))+
    geom_point(size = 1, alpha = 0.6)+
    geom_vline(xintercept = 0,col = "black")+
    scale_color_manual(values = cols,
                       labels = col_labs)+
    guides(color = guide_legend(override.aes = aes(size = 4)))+
    scale_x_continuous(labels = function(x) formatC(x, format = "g"))+
    scale_y_continuous(labels = function(x) formatC(x, format = "g"))+
    geom_text_repel(data = plot_data %>% drop_na(), aes(label = lab), family = "Helvetica", seed = 1234, size = 3, show.legend = F) +
    facet_wrap(~niche, ncol = 1)+
    labs(x = "logFC (high- vs. low-expressing)", y = "log<sub>10</sub>(Adj. <i>p</i>-value)")+
    theme_proj()+
    theme(legend.title = element_blank(), 
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in"), 
          panel.border = element_rect(linewidth = 1, color = "black"),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_markdown(size = 12),
          axis.title.y = element_markdown(size = 12),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          strip.text = element_text(size = 12, margin = margin(t = 0.05, b = 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + 
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") &
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")), 
        legend.margin = margin(l = -0.1, unit = "in"))

pdf(file.path(res_path, paste0(fig, "_panelE_microglia_DEGS.pdf")), height = 6, width = 7)
print(p)
dev.off() 


