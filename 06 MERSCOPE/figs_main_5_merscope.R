#!/usr/bin/Rscript

#===============================================================================
# 
# MICROGLIA ATLAS V3
# Figure 5
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(circlize)
library(ComplexHeatmap)
library(data.table)
library(data.table)
library(ggrastr)
library(ggrepel)
library(ggsignif)
library(ggtext)
library(grid)
library(harmony)
library(Hmisc)
library(lmerTest)
library(logger)
library(metafor)
library(parallel)
library(patchwork)
library(readxl)
library(scCustomize)
library(schex)
library(Seurat)
library(SeuratObject)
library(sf)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
data_path <- "~/data/MERSCOPE"
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope")
fig <- "fig_5_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_viz.R"))

#===============================================================================
#
# DATA
#
#===============================================================================
run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))

# REFERENCE DATA
#-------------------------------------------------------------------------------
ref_hex <- readRDS(file.path(home_path, "data", "ref_hex.rds"))
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

# MERSCOPE
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

keep_cells <- rownames(data@meta.data %>% filter(cell_type == "Microglia"))
mglia <- subset(data, subset = cell_names %in% keep_cells)
mglia <- NormalizeData(mglia, assay = "RNA") %>% ScaleData(assay = "RNA")

mglia@reductions[["umap"]] <- readRDS(file.path(home_path, "merscope/mglia/mglia.umap.rds"))
mglia@reductions[["schpf.umap"]] <- readRDS(file.path(home_path, "merscope/mglia/mglia.schpf.umap.rds"))

adjacent_counts <- readRDS(file.path(home_path, "merscope/spatial_stats/adjacent_counts.rds"))

#===============================================================================
#
# PANEL B - CELL TYPE UMAP
#
#===============================================================================
p <- DimPlot(data,
             reduction = "umap",
             group.by = "cell_type",
             order = T,
             cols = cell_type_cols) +
  guides(color = guide_legend(override.aes = aes(size = 5)))+
  theme_umap()+
  theme(plot.title = element_blank(),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.margin = margin(l = -0.2, unit = "in"),
        legend.key = element_blank())

p <- LabelClusters(p, id = "cell_type")

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 5, width = 6.5)
print(p)
dev.off()

#===============================================================================
#
# PANEL C/D - TF-GENE ASSOCIATION
#
#===============================================================================
aracne <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
enrichment <- read_rds(file.path(home_path, "aracne_sc/TF_scHPF_enrichment.rds"))

combos <- data.frame(factor = c("scHPF_20", rep("scHPF_26", 4)),
                     reg = c("IRF7", "ARID5B", "MITF", "PPARG", "CEBPA"))

res <- data.frame()
all_eff <- data.frame()
for(i in 1:nrow(combos)){
  
  genes <- names(sort(schpf@feature.loadings[, combos$factor[i]], decreasing = T))[1:20]
  genes <- genes[genes %in% rownames(data)]
  genes <- genes[!genes %in% combos$reg[i]]
  exp_targets <- aracne %>% filter(Regulator == combos$reg[i]) %>% pull(Target)
  genes <- genes[genes %in%  exp_targets]
  weight <- schpf@feature.loadings[genes, combos$factor[i]]
  
  for(niche in c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "all")){
    
    if(niche != "all"){
      keep_cells <- mglia@meta.data$cell_names[mglia@meta.data$layer_niches==niche]
    } else {
      keep_cells <- mglia@meta.data$cell_names
    }
    
    tmp <- subset(mglia, subset = cell_names %in% keep_cells)
    
    counts <- as.data.frame(t(tmp@assays$RNA@data)) %>%
      rownames_to_column("cell_names") %>%
      left_join(mglia@meta.data) %>%
      dplyr::select(cell_names, orig.ident, nCount_RNA, all_of(genes), all_of(combos$reg[i]))
    
    eff <-  mclapply(genes, mc.set.seed = 1234, FUN = function(x){
      
      tmp <- counts %>% rename_at(all_of(x), ~"gene") %>% rename_at(all_of(combos$reg[i]), ~"reg")
      
      out <- lmer(gene ~ reg + scale(nCount_RNA) + (1|orig.ident), data = tmp, na.action = na.omit)
      sum <- coefficients(summary(out))
      conf <- confint(out)["reg",]
      
      cbind.data.frame(niche, 
                       factor = combos$factor[i],
                       gene = x, 
                       pct_exp = sum(counts[,x] > 0)/nrow(counts),
                       reg = combos$reg[i], 
                       coef = sum["reg", "Estimate"], 
                       se = sum["reg", "Std. Error"], 
                       p = sum["reg", "Pr(>|t|)"],
                       lower = unname(conf[1]),
                       upper = unname(conf[2]),
                       singular = any(grepl("singular", unlist(out@optinfo$conv))))
      
    })
    
    eff <- do.call(rbind.data.frame, eff)
    eff <- left_join(eff, data.frame(gene = names(weight), weight = weight))
    m_reg <- rma(yi = eff$coef, sei = eff$se, weights = eff$weight, data = eff, method = "REML")
    
    res <- rbind(res, cbind(niche, 
                            reg = combos$reg[i],
                            factor = combos$factor[i],
                            n_cells = nrow(counts),
                            n = nrow(eff),
                            n_failed = length(genes) - nrow(eff),
                            meta.se = c(m_reg$se),
                            meta.beta = c(m_reg$beta),
                            meta.p = c(m_reg$pval),
                            meta.ci.lb = c(m_reg$ci.lb),
                            meta.ci.ub = c(m_reg$ci.ub),
                            qep = m_reg$QEp))
    
    all_eff <- rbind(all_eff, eff)
    
  }
}

res <- res %>%
  mutate_at(all_of(vars(contains("meta"))), ~as.numeric(as.character(.))) %>%
  mutate(category = ifelse(reg %in% c("IRF7","ARID5B", "MITF", "PPARG", "CEBPA"), "target", "neg")) %>%
  group_by(category, niche) %>%
  mutate(meta.p.adj = p.adjust(meta.p, method = "BH")) %>%
  ungroup()

all_eff <- all_eff %>%
  mutate(category = ifelse(reg %in% c("IRF7","ARID5B", "MITF", "PPARG", "CEBPA"), "target", "neg")) %>%
  group_by(category, niche) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))%>%
  ungroup()

write_xlsx(all_eff, file.path(res_path, "tab_TF_assoc.xslx"))
write_xlsx(res, file.path(res_path, "tab_TF_assoc_meta.xslx"))

# PLOTS
#-------------------------------------------------------------------------------
for(i in c("scHPF_20", "scHPF_26")){
  
  gene_rows <- all_eff %>%
    filter(factor == i) %>%
    mutate(coef = ifelse(coef > 0, -log10(p.adj), log10(p.adj)), 
           reg = paste0(niche, "_", reg)) %>%
    dplyr::select(reg, gene, coef)
  
  meta_rows <- res %>%
    filter(factor == i) %>%
    mutate(coef = ifelse(meta.beta > 0, -log10(meta.p.adj), log10(meta.p.adj)), 
           reg = paste0(niche, "_", reg), gene = "RE Model") %>%
    dplyr::select(reg, gene, coef) 
  
  mat <- rbind(gene_rows, meta_rows) %>%
    spread(gene, coef) %>%
    column_to_rownames("reg") %>%
    as.matrix()
  
  row_split <- gsub("_", "", str_extract(rownames(mat), "[0-9a-zA-Z\\.]+_"))
  row_split <- factor(row_split , levels = c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "all"), labels = c("WM", "GM (L2-6)", "WM + GM (L1-6)"))
  rownames(mat) <- gsub("[0-9a-zA-Z\\.]+_", "", rownames(mat))
  
  genes <- names(sort(schpf@feature.loadings[, i], decreasing = T))[1:20]
  genes <- genes[genes %in% colnames(mat)]
  mat <- mat[,c(genes, "RE Model")]
  col_split = factor(c(rep("1", length(genes)), "2"))
  
  pdf(file.path(res_path, paste0(fig, ifelse(i == "scHPF_20", "panel_C", "_panel_D"), ".pdf")), 
      height = ifelse(nrow(mat) > 3, 5, 3), width = ifelse(nrow(mat) > 3, 5.5, 4.5))
  p <- Heatmap(mat, "mat",
               col = colorRamp2(breaks = c(ifelse(min(mat, na.rm=T) < 0, min(mat, na.rm=T), 0), 0, max(mat, na.rm=T)), 
                                colors = c(ifelse(min(mat, na.rm=T) < 0, proj_cols("blue"), "white"), "white", proj_cols("red"))),
               cluster_columns = F, cluster_rows = F, row_split = factor(row_split), 
               column_split = col_split,
               row_title_gp = gpar(fontsize = 10),
               column_title_gp = gpar(fontsize = 0),
               column_names_gp = gpar(fontface = c(rep("plain", length(genes)), "bold"), fontsize = 11),
               column_names_rot = 45,
               column_gap = unit(0.1, 'in'),
               row_gap = unit(0.1, 'in'),
               na_col = "#F0F0F0",
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if(!is.na(mat[i, j]) & abs(mat[i, j]) > -log10(0.05)) {
                   grid.points(pch = 8, x, y, size = unit(0.075, "in"))
                 } else if (is.na(mat[i, j])){
                   grid.points(pch = 4, x, y, size = unit(0.075, "in"))
                 }
               },
               left_annotation = rowAnnotation(foo = anno_empty(border = FALSE, width = unit(0.1, "in"))),
               heatmap_legend_param = list(
                 title = expression('Coef-signed -log10(FDR p-value)'),
                 legend_width = unit(1.5, "in"),
                 direction = "horizontal",
                 title_position = "topcenter",
                 border = "black"
               ))
  
  draw(p,
       align_heatmap_legend = "heatmap_center",
       heatmap_legend_side = "top",
       background = "transparent", 
       padding = unit(c(0.2, 0.2, 0.2, 0.2), "in"),
       legend_title_gp = gpar(fontsize = 12, hjust = 1))
  
  for(j in 1:3){
    decorate_heatmap_body("mat", slice = j, { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1.25)) })
    decorate_annotation("foo", slice = j, { grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = "#282828", col = NA), just = "left") })
  }
  dev.off()
}

#===============================================================================
#
# PANEL E - ZOOM PLOTS
#
#===============================================================================
fov_list = c("run.26.1", "run.9.2")

gscore <- readRDS(file.path(home_path, "merscope/spatial_stats/combined_getis_ord_20.rds"))

plots <- list()
for(i in 1:2){
  
  fov <- readRDS(file.path(home_path, "merscope", run_metadata$run_name[run_metadata$batch_code==fov_list[i]], "data_proc.rds"))
  fov@meta.data <- fov@meta.data %>% 
    left_join(cell_type_annotations, by = c("fov.y"="orig.ident","cell_names")) %>%
    left_join(niches, by = c("fov.y"="orig.ident","cell_names")) %>%
    left_join(scores, by = c("cell_names"))
  rownames(fov@meta.data) <- fov@meta.data$cell_names
  
  fov <- subset(fov, subset = cell_type !="Ambiguous")
  fov <- NormalizeData(fov, assay = "Vizgen")
  
  detected_transcripts <- fread(file.path(data_path, run_metadata$run_name[run_metadata$batch_code==fov_list[i]], "cp_2D_allLayers_nuclei", "detected_transcripts.csv"))
  
  fov@meta.data <- fov@meta.data %>%
    mutate(col = ifelse(cell_type!= "Microglia", "other", scHPF_20_class),
           col = factor(col, levels = c("Mic.20.high", "Mic.20.low", "other"),
                        labels = c("Mic.20.high", "Mic.20.low", "Non-microglia")))
  
  cols <- c(proj_cols("pink", "teal"), "#F0F0F0")
  names(cols) <- c("Mic.20.high", "Mic.20.low", "Non-microglia")
  
  target_cell <- gscore %>% 
    dplyr::select(orig.ident, localG, X.log10p_adj.Sim,scHPF_20, cell_names) %>%
    left_join(adjacent_counts, by = c("cell_names"="cell_name")) %>%
    filter(factor == 20 & orig.ident == fov_list[i] & stat == "High" & X.log10p_adj.Sim > -log10(0.05) & localG > 0) %>%
    arrange(desc(localG))
  
  zoom <- fov@meta.data %>%
    left_join(GetTissueCoordinates(fov), by = c("cell_names"="cell")) %>%
    filter(cell_names %in% target_cell$cell_names[1]) %>%
    dplyr::select(x, y)
  zoom <- unname(unlist(zoom))
  
  p <- fov@meta.data %>%
    left_join(GetTissueCoordinates(fov), by = c("cell_names"="cell")) %>%
    arrange(desc(col)) %>%
    ggplot(aes(x = x, y = y, color = col)) + 
    geom_point_rast(size = 0.1, raster.dpi = 600, scale = 0.5)+
    scale_color_manual(values = cols)+
    guides(color = guide_legend(override.aes = aes(size = 5)))+
    annotate("rect", xmin = zoom[1]-100, xmax = zoom[1]+100,
             ymin = zoom[2]-100, ymax = zoom[2]+100,
             size = 0.25, fill = NA, color = "black") +
    theme_umap()+
    theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          legend.position = "none") +
    NoAxes()
  
  plots <- c(plots, list(p))
  
  for(j in c("scHPF_20", "IRF7", "MX1", "OAS3")){
    p <- FOV_FeaturePlot(data = fov, 
                        fov = fov_list[i], 
                        zoom = c(zoom[1]-100, zoom[1]+100,zoom[2]-100, zoom[2]+100),
                        boundary = 200,
                        feature = j, 
                        bar.text.size = 0)+
      labs(title = j)+
      theme(plot.title = element_text(size = 8, face = "plain", margin = margin(b = 0.05, unit = "in")),
            legend.margin = margin(l = -0.2, unit = "in"),
            legend.text = element_text(size = 8, margin = margin(l = 0.05, unit = "in")),
            legend.key.width = unit(0.1, "in"),
            legend.key.height = unit(0.2, "in")) +
      NoAxes()
    plots <- c(plots, list(p))
  }
}

p <- wrap_plots(plots) + plot_layout(nrow = 2) & theme(plot.margin = margin(l = 0, r = 0, unit = "in"))
pdf(file.path(res_path, paste0(fig, "panel_E.pdf")), height = 3.5, width = 7)
print(p)
dev.off()

#===============================================================================
#
# PANEL F - IRF7/MX1/OAS3 in MICROGLIA
#
#===============================================================================
plot_data <- cbind(data@meta.data, t(data@assays$RNA@counts[c("IRF7", "MX1", "OAS3"),])) %>%
  filter(cell_type == "Microglia") %>%
  mutate(IRF7 = ifelse(IRF7 > 0, "IRF7+", "IRF7-"),
         MX1 = ifelse(MX1 > 0, "MX1+", "MX1-"),
         OAS3 = ifelse(OAS3 > 0, "OAS3+", "OAS3-")) %>%
  rowwise() %>%
  mutate(stat = ifelse(IRF7 == "IRF7+", paste0(IRF7,"/", MX1, "/", OAS3), IRF7),
         scHPF_20_class = recode(scHPF_20_class, "Mic.20.low"="IR<sup>Low</sup>", "Mic.20.high"="IR<sup>High</sup>"),
         scHPF_20_class = factor(scHPF_20_class, levels = c("IR<sup>Low</sup>", "IR<sup>High</sup>")),
         dx = ifelse(donor %in% c(26, 33), "AD", "non-AD"))

plots <- list()
for(j in c("AD", "non-AD")){
  tmp <- plot_data  %>% filter(dx == j) 
  chisq <- chisq.test(tmp$scHPF_20_class=="IR<sup>High</sup>", tmp$IRF7=="IRF7+")
  fisher <- fisher.test(tmp$scHPF_20_class=="IR<sup>High</sup>", tmp$IRF7=="IRF7+")
  
  or <- unname(formatC(fisher$estimate, digits = 3, format = "g"))
  upper <- unname(format(fisher$conf.int[2], digits = 3, format = "g"))
  lower <- unname(format(fisher$conf.int[1], digits = 3, format = "g"))
  lab <- paste0("OR (95%CI) = ", or, " (", lower, ", ", upper, ")", ifelse(chisq$p.value < 0.001, "***", ""))
  
  p <- plot_data %>%
    filter(dx == j) %>%
    group_by(scHPF_20_class, stat) %>%
    summarise(n = n()) %>%
    group_by(scHPF_20_class) %>%
    summarise(dx="AD", stat, pct = n/sum(n)) %>%
    ggplot(aes(x = scHPF_20_class, y = pct, fill = stat)) +
    geom_col() +
    geom_text(aes(label = ifelse(pct > 0.05, paste0(round(pct*100, 1), "%"), NA)), 
              color = 'white', 
              position = position_stack(vjust = 0.5))+
    scale_fill_manual(values = unname(proj_cols("blue", "purple","teal","yellow", "pink")))+
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25), labels = formatC(seq(0, 1, 0.25), format = "g"))+
    labs(y = "Proportion", title = j, subtitle = lab)+
    theme_proj()+
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_markdown(size = 12, margin = margin(t = 0, unit = "in")),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_blank(), 
          plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          plot.title = element_text(size = 14, margin = margin(b = 0.05, t = 0.05, unit = "in")),
          plot.subtitle = element_markdown(size = 10, hjust = 0.5))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(nrow = 2, guides = "collect") &
  theme(legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
        legend.margin = margin(l = -0.2, unit = "in"),
        legend.key.spacing.x = unit(0.05, "in"))

pdf(file.path(res_path, paste0(fig, "panel_F.pdf")), height = 6, width = 5)
print(p)
dev.off()

#===============================================================================
#
# PANEL I - DEGS IN OTHER CELL TYPES
#
#===============================================================================
cols <- c(unname(proj_cols("blue", "pink")), unname(proj_cols_grey("med grey")))
names(cols) <- c("1", "2", "3")
col_labs <- c("Downregulated", "Upregulated", "NS (Adj. p > 0.05)")

layer_labs <- c("WM", "GM (L1)", "GM (L2-6)")
names(layer_labs) <- c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")

plot_data <- read_xlsx(file.path(res_path, "spatial_stats/tab_adjacent_DEGs20.xlsx")) %>% 
  mutate(log10.p.adj = -log10(p_val_adj),
         col = ifelse(avg_log2FC > 0 & p_val_adj < 0.05, "2", "3"),
         col = ifelse(avg_log2FC < 0 & p_val_adj < 0.05, "1", col),
         lab = ifelse(p_val_adj < 0.05, genes, NA),
         niche = recode(niche, !!!layer_labs)) %>%
  filter(niche != "WM" & cell_type == "Astrocytes")

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
  geom_text_repel(data = plot_data %>% drop_na(), aes(label = lab), 
                  family = "Helvetica", seed = 1234, size = 3, show.legend = F) +
  labs(x = "Log FC (IR<sup>High</sup> vs. IR<sup>Low</sup>)", 
       y = "log<sub>10</sub>(Adj. <i>p</i>-value)")+
  facet_wrap(~niche, nrow = 2) +
  theme_proj()+
  theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
        legend.spacing.y = grid::unit(0.05, "in"),,
        legend.margin = margin(l = -0.2, unit = "in"),
        panel.border = element_rect(linewidth = 1, color = "black"),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12),
        axis.line.y = element_line(linewidth = .25, color = "black"),
        axis.line.x = element_line(linewidth = .25, color = "black"),
        strip.text = element_text(size = 10, margin = margin(t = 0.05, b = 0.05, unit = "in")))

pdf(file.path(res_path, paste0(fig, "panel_G.pdf")), height = 5, width = 4)
print(p)
dev.off()

#===============================================================================
#
# PANEL H
#
#===============================================================================
pred <- data.frame()
for(i in 1:12){
  for(j in  c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    file <- file.path(res_path, "nearest_neighbors", paste0(run_metadata$batch_code[i], "_", j, "_other_cells_seurat.rds"))
    if(file.exists(file)){
      tmp <- readRDS(file)
      pred <- rbind(pred, cbind(batch_code = run_metadata$batch_code[i],
                               donor = run_metadata$donor[i],
                               tmp@misc$localG_perm_scHPF_20))
      rm(tmp)
    }
  }
}

fov_list <- pred %>% 
  filter(!batch_code %in% c("run.1.1", "run.9.1")) %>%
  group_by(donor, batch_code) %>% 
  summarise(mean=mean(localG.pred)) %>%
  group_by(donor) %>%
  arrange(desc(mean)) %>%
  slice_head(n= 1) %>%
  pull(batch_code)
fov_list <- factor(fov_list, levels = run_metadata$batch_code)
fov_list <- sort(fov_list) 

plots <- list()
for(i in 1:4){
  
  fov <- readRDS(file.path(home_path, "merscope", run_metadata$run_name[run_metadata$batch_code==fov_list[i]], "data_proc.rds"))
  
  fov@meta.data = fov@meta.data %>%
    left_join(data@meta.data %>% dplyr::select(cell_names, cell_type, layer_niches, scHPF_20_class))
  
  pred <- data.frame()
  for(j in  c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    file <- file.path(res_path, "nearest_neighbors", paste0(fov_list[i], "_", j, "_other_cells_seurat.rds"))
    if(file.exists(file)){
      tmp <- readRDS(file)
      pred <- rbind(pred, tmp@misc$localG_perm_scHPF_20)
      rm(tmp)
    }
  }
  rownames(pred)<-NULL
  fov@meta.data <- fov@meta.data %>%
    left_join(pred, by = c("cell_names")) 
  rownames(fov@meta.data) = fov@meta.data$cell_names
  
  fov <- subset(fov, subset = cell_type !="Ambiguous")
  
  detected_transcripts <- fread(file.path(data_path, run_metadata$run_name[run_metadata$batch_code==fov_list[i]], "cp_2D_allLayers_nuclei", "detected_transcripts.csv"))
  
  fov@meta.data <- fov@meta.data %>%
    mutate(col = ifelse(cell_type == "Microglia", scHPF_20_class, as.character(cell_type)),
           col = recode(col, !!! cell_type_abb),
           col = factor(col, levels = c("Mic.20.high", "Mic.20.low", cell_type_abb)))
  
  cols <- c(proj_cols("pink"), proj_cols("teal"), alpha(cell_type_cols, 0.5))
  names(cols)<-c("Mic.20.high", "Mic.20.low", dplyr::recode(names(cell_type_cols), !!! cell_type_abb))
  cols["Oli"] <- "#EFE4EA"
  
  target_cell <- fov@meta.data %>%
    filter(cell_type=="Astrocytes") %>%
    arrange(desc(localG.pred)) %>%
    dplyr::select(cell_names, layer_niches, localG.pred) %>%
    left_join(GetTissueCoordinates(fov), by = c("cell_names"="cell"))
  
  target_cell1 <- target_cell[1, ]
  target_cell2 <- target_cell %>% 
    filter(!cell_names %in% target_cell1$cell_names &
             (x > target_cell$x[1]+100 | x < target_cell$x[1]-100) &
             (y > target_cell$y[1]+100 | y < target_cell$y[1]-100))
  target_cell2 = target_cell2[1, ]
  
  layers <- c(target_cell1$layer_niches, target_cell2$layer_niches)
  
  zoom_list <- list(c(target_cell1$x[1]-100, target_cell1$x[1]+100, target_cell1$y[1]-100, target_cell1$y[1]+100),
                   c(target_cell2$x[1]-100, target_cell2$x[1]+100, target_cell2$y[1]-100, target_cell2$y[1]+100))
  
  p <- fov@meta.data %>%
    left_join(GetTissueCoordinates(fov), by = c("cell_names"="cell")) %>%
    arrange(desc(col)) %>%
    ggplot(aes(x = x, y = y, color = col)) +
    geom_point_rast(size = 0.1, scale = 0.25)+
    scale_color_manual(values = cols)+
    guides(color = guide_legend(override.aes = aes(size = 5)))+
    labs(title = run_metadata$donor_name[run_metadata$batch_code==fov_list[i]])+
    annotate("rect", xmin = zoom_list[[1]][1], xmax =  zoom_list[[1]][2], ymin =  zoom_list[[1]][3], ymax = zoom_list[[1]][4], size = 0.5, fill = NA, color = "black")+
    annotate("text", label = "1", vjust = 0, hjust = 1.2, x = zoom_list[[1]][1], y =  zoom_list[[1]][3], size = 2.5, color = "black", family = "Helvetica", fontface = "bold")+
    annotate("rect", xmin = zoom_list[[2]][1], xmax =  zoom_list[[2]][2], ymin =  zoom_list[[2]][3], ymax = zoom_list[[2]][4], size = 0.5, fill = NA, color = "black")+
    annotate("text", label = "2", vjust = 0, hjust = 1.2, x = zoom_list[[2]][1], y =  zoom_list[[2]][3], size = 2.5, color = "black", family = "Helvetica", fontface = "bold")+
    theme_umap()+
    theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          legend.position = "none", 
          plot.title = element_text(size = 10, face = "bold", family = "Helvetica", margin = margin(t = 0, b = 0, unit = "in"))) +
    NoAxes()
  
  plots <- c(plots, list(p))
  
  for(k in 1:2){
    
    p <- FOV_DimPlot(data= fov, 
                    fov = as.character(fov_list[i]), 
                    zoom = zoom_list[[k]],
                    group.by = "col", 
                    boundary = 200,
                    cols = cols, 
                    mol_cols = unname(c("cyan", proj_cols("blue"), "#282828", "red", "purple", "forestgreen")),
                    fill.alpha = 0.7,
                    molecules = c("IRF7", "MX1", "OAS3", "IFI6", "IFIT3", "IFI44L"),
                    pt.size = 0.1,
                    pt.scale = 0.5,
                    detected_transcripts = detected_transcripts, 
                    bar.text.size = 0) +
      scale_y_continuous(limits = zoom_list[[k]][3:4], breaks = seq(zoom_list[[k]][3], zoom_list[[k]][4], 50))+
      scale_x_continuous(limits = zoom_list[[k]][1:2], breaks = seq(zoom_list[[k]][1], zoom_list[[k]][2], 50))+
      labs(title = paste0(k, ": ", recode(layers[k], !!!layer_labs)))+
      scale_fill_manual(values = cols, drop = FALSE) +
      theme(legend.key.spacing.y = grid::unit(0, unit = "in"),
            legend.key.width = grid::unit(0.1, unit = "in"),
            legend.key.height = grid::unit(0.05, unit = "in"),
            legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
            plot.margin = margin(0, 0, 0, 0, unit = "in"),
            plot.title = element_text(size = 10, face = "plain", family = "Helvetica", margin = margin(t = 0, b = 0, unit = "in"))) +
      NoAxes()
    
    plots <- c(plots, list(p))
  }
}

p <- wrap_plots(plots) + plot_layout(nrow = 2, guides = "collect")
pdf(file.path(res_path, paste0(fig, "panel_I.pdf")), height = 3.5, width = 10)
print(p)
dev.off()
