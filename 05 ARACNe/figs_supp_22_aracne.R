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
library(clusterProfiler)
library(ComplexHeatmap)
library(ggrastr)
library(ggVennDiagram)
library(grid)
library(patchwork)
library(readxl)
library(scCustomize)
library(tidyverse)
library(viridis)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sn")
fig <- "fig_s22_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# scHPF
#-------------------------------------------------------------------------------
schpf <- read_rds(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

# REGULON SUMAMRIES
#-------------------------------------------------------------------------------
sc <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
sn <- read_csv(file.path(home_path, "aracne_sn/aracne_w_tfmode.csv"))

sc %>% group_by(Regulator) %>% summarise(n = n()) %>% summarise(mean_n=mean(n), median(n), max(n), min(n))
sn %>% group_by(Regulator) %>% summarise(n = n()) %>% summarise(mean_n=mean(n), median(n), max(n), min(n))
universe <- rownames(schpf@feature.loadings)

#===============================================================================
#
# PANEL A
#
#===============================================================================
p <- ggVennDiagram(list(sc = unique(sc$Regulator), sn= unique(sn$Regulator)), 
                   label_alpha = 0, 
                   label_color = "white") +
  scale_fill_gradient(high = proj_cols('purple'), low = proj_cols("blue"))+
  theme(legend.position = "none")+
  coord_flip()

pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 2, width = 3)
print(p)
dev.off()
#===============================================================================
#
# PANEL A2
#
#===============================================================================

# AVERAGE REGULON SIZE
#-------------------------------------------------------------------------------
plot_data <- rbind(sc %>% mutate(source = "sc"), sn %>% mutate(source = "sn"))

p1 <- plot_data %>% 
  group_by(source, Regulator, dir) %>% 
  summarise(n = n()) %>%
  group_by(source, dir) %>% 
  summarise(mean = mean(n)) %>%
  ggplot(aes(x = factor(dir, labels = c("Down", "Up")), y = mean, fill = factor(source, levels = c("sc", "sn"), labels = c("scRNA-seq", "snucRNA-seq")))) + 
  geom_col(position = position_dodge()) + 
  scale_fill_manual(values = unname(proj_cols("blue", "yellow")))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25), limits = c(0, 100))+
  labs(y = "Mean N Targets")+
  theme_proj()+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, margin = margin(t = 0, unit = "in"), hjust = 1, angle = 45),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.1, "in"),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
        legend.key.spacing.x = unit(0.05, "in"),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top',
        legend.direction = "horizontal") 

# REGULON COMPLEXITY
#-------------------------------------------------------------------------------
all_regs <- unique(c(sc$Regulator, sn$Regulator))

p2 <- plot_data %>% 
  group_by(source, Regulator, dir) %>% 
  summarise(n = sum(Target %in% all_regs)) %>%
  group_by(source, dir) %>%
  summarise(mean = mean(n)) %>%
  ggplot(aes(x = factor(dir, labels = c("Down", "Up")), y = mean, fill = factor(source, levels = c("sc", "sn"), labels = c("scRNA-seq", "snucRNA-seq")))) + 
  geom_col(position = position_dodge())+
  scale_fill_manual(values = unname(proj_cols("blue", "yellow")))+
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10, 2), limits = c(0, 10))+
  labs(y = "Mean N TFs")+
  theme_proj()+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, margin = margin(t = 0, unit = "in"), hjust = 1, angle = 45),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.1, "in"),
        legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
        legend.key.spacing.x = unit(0.05, "in"),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top',
        legend.direction = "horizontal") 

# OVERLAP BETWEEN TARGET GENES
#-------------------------------------------------------------------------------
overlap_regs <- intersect(sc$Regulator, sn$Regulator)

res <- data.frame()
for(i in overlap_regs){
  for(j in c("up", "down")) {
    unique_sc_targets <- sc %>% filter(Regulator == i & dir == j) %>% pull(Target)
    unique_sn_targets <- sn %>% filter(Regulator == i & dir == j) %>% pull(Target)
    
    n_overlap <- length(intersect(unique_sc_targets, unique_sn_targets))
    n_unique <- length(unique(c(unique_sc_targets, unique_sn_targets)))
    
    res <- rbind(res, cbind(reg = i, dir = j, pct = n_overlap/n_unique))
  }
}


p3 <- ggplot(res, aes(x = factor(dir, labels = c("Down", "Up")), y = as.numeric(pct))) + 
  geom_boxplot(fill = "lightgrey") +
  labs(y = "% overlapping\ntargets")+
  theme_proj()+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, margin = margin(t = 0, unit = "in"), hjust = 1, angle = 45),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank()) 

p <- wrap_plots(p1, p2, p3) + plot_layout(guides = "collect") & theme(legend.position = 'top')

pdf(file.path(res_path, paste0(fig, "panel_A2.pdf")), height = 2, width = 4)
print(p)
dev.off()
#===============================================================================
#
# PANEL B
#
#===============================================================================
aracne_sc <- read_xlsx(file.path(home_path, "aracne_sc", "tab_aracne_scHPF_enrichment.xlsx"))
aracne_sn <- read_xlsx(file.path(home_path, "aracne_sn", "tab_aracne_scHPF_enrichment.xlsx"))

top_regs <- aracne_sc %>% 
  group_by(ID) %>% 
  top_n(-3, wt = p.adj.bonf) %>% 
  pull(signature)
top_regs <- unique(top_regs)

mat <- aracne_sn %>% 
  filter(signature %in% top_regs) %>%
  dplyr::select(ID, p.adj.bonf, signature) %>%
  mutate(p.adj.bonf = -log10(p.adj.bonf),
         ID = factor(ID, levels = factors, labels = names(factors))) %>%
  spread(signature, p.adj.bonf) %>%
  column_to_rownames("ID") %>%
  as.matrix()

cairo_pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 6, width = 16, bg = "transparent")
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(0, max(mat, na.rm = T)),  colors = c("white", proj_cols("red"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & (mat[i, j] >-log10(0.05) | mat[i, j] <log10(0.05))) {
                 grid.points(pch = 8, x, y, size = unit(0.075, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             cluster_rows = F,
             row_title = "Expression Programs",
             row_title_side = "left",
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 16),
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 14),
             row_names_side = "left",
             row_dend_side = "right",
             column_title = "sc-ARACNE TF Regulons",
             column_title_side = "bottom",
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 16),
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 14),
             column_names_rot = 45,
             cluster_columns = F,
             column_names_centered = F,
             heatmap_legend_param = list(
               title = expression('-log'[10]*'(Adj. '*italic(p)*'-value)'), 
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
decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1.5)) })
dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
regulators <- c("ARID5B", "CEBPA", "MITF", "PPARG", "IRF7")

terms <- sc  %>% 
  filter(dir == "up") %>% dplyr::select(gsid = Regulator, Gene = Target) %>% 
  filter(gsid %in% regulators)
# this is a temporary solution to: https://github.com/YuLab-SMU/clusterProfiler/issues/283
# to avoid truncation of the universe
# terms = rbind(terms, cbind(gsid = "universe", Gene = universe))

res <- data.frame()
for(i in regulators){
  
  genes <- sn %>% filter(Regulator %in% i) %>% filter(dir == "up") %>% pull(Target)
  
  set.seed(1234)
  out <- enricher(genes, 
                  universe = universe, 
                  TERM2GENE = terms,
                  pAdjustMethod = "none",
                  minGSSize = 3,
                  maxGSSize = 2500)
  if(!is.null(out)) {
    out <- out@result %>% mutate(sn_reg = i)
    res <- rbind(res, out)
  }
}
res <- res %>% 
  mutate(padj = p.adjust(pvalue, method = "BH"),
         richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))

mat <- res %>%
  dplyr::select(sn_reg, ID, richFactor) %>%
  spread(ID, richFactor) %>%
  column_to_rownames("sn_reg") %>%
  as.matrix()

pmat <- res %>%
  mutate(lab = ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", ""))),
         lab = paste0(round(richFactor, 2), lab)) %>%
  dplyr::select(sn_reg, ID, lab) %>%
  spread(ID, lab) %>%
  column_to_rownames("sn_reg") %>%
  as.matrix()

pdf(file.path(res_path, paste0(fig, "panel_C.pdf")), height = 3.5, width = 3.5)
p <- Heatmap(mat, name = "mat", 
             col = colorRamp2(breaks = c(0, max(mat, na.rm=T)),  colors = c("white", proj_cols("light blue"))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(pmat[i, j])) {
                 grid.text(pmat[i, j], x, y, gp = gpar(fontsize = 8))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             row_dend_side = "right",
             row_title = "Validation ARACNe (sn)",
             row_names_side = "left",
             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 10),
             column_title = "Discovery ARACNe (sc)",
             column_title_side = "bottom",
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize= 10),
             column_names_centered = F,
             heatmap_legend_param = list(
               title = "Fold Enrichment", 
               legend_width = unit(1.5, "in"),
               direction = "horizontal",
               title_position = "topcenter",
               border = "black"
             ))

draw(p, heatmap_legend_side = "top",padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
dev.off()

#===============================================================================
#
# PANEL D
#
#===============================================================================
source(file.path(home_path, "src/00 Helpers/load_ref_data.R"))
data <- NormalizeData(data) %>% ScaleData(features = c("ARID5B", "CEBPA", "MITF", "PPARG"))

snuc <- readRDS("~/projects/microglia_atlas/proj_snuc/projection/snuc.rds")
snuc <- NormalizeData(snuc) %>% ScaleData(features = c("ARID5B", "CEBPA", "MITF", "PPARG"))

labs <- c("scHPF_26", "ARID5B", "CEBPA", "MITF", "PPARG")
names(labs) <- c("<i>GPNMB</i><sup>high</sup>", "<i>ARID5B</i>", "<i>CEBPA</i>", "<i>MITF</i>", "<i>PPARG</i>")

plots <- list()
for(i in c("data", "snuc")){
  for(j in labs) {
    
    p <- Plot_Density_Custom(eval(parse(text = i)), reduction = 'schpf.umap', features = j)
    
    p <- ggplot(p$data, aes(x = schpfumap_1, y = schpfumap_2, color = feature)) +
      geom_point_rast(size = 0.1, scale = 0.5)+
      scale_color_viridis(option = "magma")+
      labs(title = names(labs)[labs %in% j]) +
      theme_umap()+
      theme(legend.position = "right",
            legend.key.height = unit(0.25, "in"),
            legend.key.width = unit(0.2, "in"),
            legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
            legend.margin = margin(l = -0.2, unit = "in"),
            legend.title = element_blank(),
            plot.title = element_markdown(size = 16, face = "bold", family = "Helvetica", margin = margin(b = 0)),
            panel.border = element_rect(linewidth = 1.5, color = "black"),
            axis.line.x = element_line(linewidth = .25, color = "black"),
            axis.line.y = element_line(linewidth = .25, color = "black"),
            plot.margin = margin(0.01, 0.01, 0.01, 0.01, 'in'))
    plots <- c(plots, list(p))
  }
}

plots[[1]] <- plots[[1]] + labs(y = "Discovery<br>(scRNA-seq)") +
  theme(axis.title.y = element_markdown(size = 16, face = "bold", family = "Helvetica", angle = 90))
plots[[6]] <- plots[[6]] + labs(y = "CUIMC1<br>(snRNA-seq)")+
  theme(axis.title.y = element_markdown(size = 16, face = "bold", family = "Helvetica", angle = 90))

p <- wrap_plots(plots) + plot_layout(nrow = 2)

pdf(file.path(res_path, paste0(fig, "panel_D.pdf")), height = 4.5, width = 12)
print(p)
dev.off()
