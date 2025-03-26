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
library(fgsea)
library(ggtext)
library(patchwork)
library(rms)
library(SeuratObject)
library(Seurat)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "proj_snuc/analysis")
fig <- "fig_s16"

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
# PANEL A
#
#===============================================================================
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

# LISTS
#-------------------------------------------------------------------------------
fgsea_lists <- list()
for(i in 1:length(factors)){
  fgsea_lists[[i]] <- names(sort(schpf@feature.loadings[,factors[i]], decreasing = T))[1:100]
  names(fgsea_lists)[i] <- factors[i]
}

# FGSEA
#-------------------------------------------------------------------------------
keep_cells <- snuc@meta.data %>% filter(celltype == "Micr" & !is.na(annotation)) %>% pull(cell_names)
snuc_mini <- subset(snuc, subset = cell_names %in% keep_cells)

cell_types <- na.omit(unique(snuc_mini@meta.data$annotation))

snuc_mini <- NormalizeData(snuc_mini, assay = "RNA") %>%
  ScaleData(assay = "RNA", features = rownames(schpf@feature.loadings))

res <- data.frame()
n_sig <- data.frame()
for(i in cell_types){
  
  snuc_mini@meta.data$col <- ifelse(snuc_mini@meta.data$annotation == i, i, "Other")
  Idents(snuc_mini) <- "col"
  
  markers <- FindMarkers(snuc_mini,
                         ident.1 = i,
                         ident.2 = "Other",
                         only.pos = F,
                         assay = "RNA",
                         slot = "data",
                         max.cells.per.ident = 2000,
                         min.cells.group = 50,
                         test.use = "MAST",
                         features = rownames(schpf@feature.loadings),
                         random.seed = 1234,
                         latent.vars = c("batch", "nCount_RNA"))
  
  tmp <- markers %>%
    filter(p_val_adj < 0.05) %>%
    arrange(avg_log2FC)
  
  genes <- tmp$avg_log2FC
  names(genes) <- rownames(tmp)
  
  n_sig <- rbind(n_sig, cbind(cell_type = i, n = nrow(tmp), na_fail = sum(is.infinite(genes))))
  
  genes = genes[!is.infinite(genes)]
  
  set.seed(123)
  fgsea_res <- fgsea(pathways = fgsea_lists,
                     stats = genes,
                     minSize = 3,
                     maxSize = 700) %>%
    mutate(cell_type = i)
  
  res <- rbind(res, fgsea_res)
}

res <- res %>%
  mutate(pathway = factor(pathway, levels = factors, labels = names(factors)))

write_csv(res, file.path(home_path, "proj_snuc/analysis/tab_rosmap_fgsea.csv"))

# HEATMAP
#-------------------------------------------------------------------------------
mat <- res %>%
  mutate(log10padj = ifelse(NES > 0, -log10(padj), log10(padj))) %>%
  dplyr::select(pathway, cell_type, log10padj) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c("Homeostatic", "Homeostatic-Redox", "Stress-Response", "Reactive" ,
                                       "Reactive-Redox", "Proliferating", "Interferon-Response",
                                       "Disease-Elevated", "Mic.15", "Mic.16")),
         pathway = factor(pathway, levels = names(factors))) %>%
  spread(cell_type, log10padj) %>%
  column_to_rownames("pathway") %>%
  mutate_all(~as.numeric(.)) %>%
  as.matrix() %>%
  t()

cairo_pdf(file.path(home_path, "proj_snuc/analysis", paste0(fig, "_panelA.pdf")), height = 4.5, width = 8)
set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(min(mat, na.rm = T), 0, max(mat, na.rm = T)),
                              colors = c(proj_cols("blue"), "white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & abs(mat[i, j]) >-log10(0.05)) {
                 grid.points(pch = 8, x, y, size = unit(0.1, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             cluster_columns = F,
             cluster_rows = T,
             row_title = "CUIMC1-derived microglia subtypes",
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             row_names_side = "left",
             row_dend_side = "right",
             column_title = "Expression programs",
             column_title_side = "bottom",
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             heatmap_legend_param = list(
               title = expression('NES-signed -log'[10]*'(FDR '*italic(p)*'-value)'),
               legend_width = unit(1.5, "in"),
               direction = "horizontal",
               title_position = "lefttop",
               border = "black",
               labels_gp = gpar(font = 14, fontfamily = "Helvetica"),
               title_gp = gpar(font = 14,  fontfamily = "Helvetica")
             ))
draw(p, 
     align_heatmap_legend = "heatmap_center", 
     heatmap_legend_side = "top",
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "white", lwd = 2)) })
dev.off()

#===============================================================================
#
# PANEL B - RIDGE PLOT
#
#===============================================================================
plots <- list()
for(i in factors){
  
  tmp <- snuc_agg %>% 
    rename_at(all_of(i), ~"factor") %>%
    mutate(facet = names(factors)[factors %in% i],
           facet = gsub('-high', paste0("<sup>high</sup>"), facet))
  
  p <- ggplot(tmp, aes(x = factor)) + 
    geom_density(fill = proj_cols("light blue"))+
    scale_x_continuous(labels = function(x) formatC(x, digits = 1,  format = "g")) + 
    geom_vline(xintercept = 0.01, lty = 2) + 
    facet_wrap(~facet)+
    labs(x = "Donor-aggregated scores", y = "Density")+
    theme_proj() +
    theme(strip.text = element_markdown(size = 8, margin = margin(b = 0.05, t = 0.05, unit = "in"), family = "Helvetica"),
          axis.text = element_text(size = 10, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          panel.grid.major = element_blank())
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(guides = "collect", ncol = 4)

cairo_pdf(file.path(res_path, "panel_B.pdf"), height = 9, width = 7.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
res <- data.frame()
for(trait in vars){
  for(i in paste0("scHPF_", c(2, 4, 5, 10, 11, 19, 25, 26))){

    tmp <- snuc_agg %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      drop_na()

    d <- datadist(tmp); options(datadist = "d")
    mod <- ols(trait ~ study + age_death + msex + pmi + factor, data = tmp)

    res <- rbind(res,
                 cbind(factor = i,
                       trait,
                       n = nrow(tmp),
                       mod_p = unname(1-pchisq(mod$stats[2], mod$stats[3])),
                       coef = coefficients(summary.lm(mod))["factor","Estimate"],
                       se = coefficients(summary.lm(mod))["factor","Std. Error"],
                       p = summary.lm(mod)$coefficients["factor", "Pr(>|t|)"]))

    rm(list = c("tmp", "d", "mod"))
  }
}

res <- res %>%
  mutate_at(all_of(c("mod_p","coef", "se", "p")), ~as.numeric(as.character(.))) %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni"),
         log10.padj = ifelse(coef > 0, -log10(p.adj), log10(p.adj)))

write_csv(res, file.path(home_path, "proj_snuc/analysis/tab_rosmap_detailed_neuropathology.csv"))

mat <- res %>%
  mutate(log10.padj = -log10(p.adj),
         log10.padj = ifelse(coef > 0, log10.padj, -log10.padj)) %>%
  mutate(trait = factor(trait, levels = vars)) %>%
  dplyr::select(factor, var = trait, log10.padj) %>%
  mutate(factor = factor(factor, levels = factors, labels = names(factors))) %>%
  spread(var, log10.padj) %>%
  column_to_rownames("factor") %>%
  mutate_all(~as.numeric(as.character(.))) %>%
  as.matrix()

col_split <- data.frame(vars = colnames(mat)) %>%
  mutate(trait = str_extract(vars, "amyloid|tangles|plaq_d|plaq_n|nft")) %>%
  mutate(trait = factor(trait, levels = c("amyloid", "plaq_d", "plaq_n", "tangles", "nft"),
                        labels = c("Amyloid", "Diffuse Plaque", "Neuritic Plaque", "Tangles", "NFT")))

pdf(file.path(home_path, "proj_snuc/analysis", paste0(fig, "_panel_C.pdf")), height = 3.5, width = 10.5)
set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(min(mat, na.rm=T), 0, max(mat, na.rm=T)),
                              colors = c(proj_cols("blue"),"white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & abs(mat[i, j]) >-log10(0.05)) {
                 grid.points(pch = 8, x, y, size = unit(0.1, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             cluster_rows = T,
             row_names_side = "left",
             row_dend_side = "right",
             top_annotation = HeatmapAnnotation(foo = anno_empty(border = FALSE, height = unit(0.5, "mm"))),
             column_split = col_split$trait,
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_names_rot = 45,
             column_names_gp = gpar(fontfamily = "Helvetica", fontsize = 12),
             column_title_gp = gpar(fontfamily = "Helvetica", fontsize = 12),
             cluster_columns = F,
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('Coefficient-signed -log'[10]*'(Adj. '*italic(p)*'-value)'),
               legend_width = unit(1.5, "in"),
               direction = "horizontal",
               title_position = "lefttop",
               border = "black",
               labels_gp = gpar(font = 14, fontfamily = "Helvetica"),
               title_gp = gpar(font = 14,  fontfamily = "Helvetica")
             ))
draw(p, 
     align_heatmap_legend = "heatmap_center", 
     heatmap_legend_side = "top",
     background = "transparent", 
     padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))

for(i in 1:5){

  decorate_annotation("foo", slice = i, {
    grid.rect(x = 1, height = unit(2, "mm"), gp = gpar(fill = "#282828", col = NA), just = "right")
  })
}

dev.off()
