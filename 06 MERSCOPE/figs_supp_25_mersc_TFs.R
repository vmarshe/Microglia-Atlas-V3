#!/usr/bin/Rscript

#===============================================================================
# 
# MERSCOPE
# Figure S25. TF and ARACNE-predicted target expression across GPNMBHigh and GPNMBLow microglia. 
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib_spat")

library(ggrepel)
library(ggtext)
library(lmerTest)
library(patchwork)
library(readxl)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope")

figs <- "fig_s25"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))
source(file.path(home_path, "src/00 Helpers/funcs_seurat.R"))

# DATA
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

#===============================================================================
#
# PANEL A
#
#===============================================================================
plot_data <- cbind(data@meta.data, t(data@assays$RNA@counts[c("CEBPA", "ARID5B", "MITF", "PPARG"), ])) %>%
  filter(cell_type == "Microglia") %>%
  mutate(CEBPA = ifelse(CEBPA > 0, 1, 0),
         ARID5B = ifelse(ARID5B > 0, 1, 0),
         MITF = ifelse(MITF > 0, 1, 0),
         PPARG = ifelse(PPARG > 0, 1, 0)) %>%
  rowwise() %>%
  mutate(counts = sum(CEBPA, ARID5B, MITF, PPARG),
         scHPF_26_class = factor(scHPF_26_class, levels = c("Mic.26.low", "Mic.26.high")),
         dx = ifelse(donor %in% c(26, 33), "AD", "non-AD"))

plots <- list()
for(i in c("AD", "non-AD")){
  tmp <- plot_data %>% filter(dx == i) %>% mutate(scHPF_26_class = ifelse(scHPF_26_class == "Mic.26.high", 1, 0))
  annot <- coefficients(summary(lmer(counts ~ scHPF_26_class + scale(nCount_RNA)+ (1|orig.ident), data = tmp)))
  annot <- paste0("*β* (se) = ", 
                  round(annot[2, 1], 2), " (", round(annot[2, 2], 2), "), *p* = ",
                  formatC(annot[2, 5], digits = 2, format = "e"))
  
  p <- plot_data %>%
    filter(dx == i) %>%
    group_by(scHPF_26_class, counts) %>%
    summarise(n = n()) %>%
    group_by(scHPF_26_class) %>%
    reframe(counts, pct = n/sum(n)) %>%
    ggplot(aes(x = factor(scHPF_26_class, 
                          levels = c("Mic.26.low", "Mic.26.high"), 
                          labels = c("GPNMB<sup>Low</sup>", "GPNMB<sup>High</sup>")), 
               y = pct, 
               fill = factor(counts))) +
    geom_col() +
    geom_text(aes(label = ifelse(pct > 0.05, paste0(round(pct*100, 1), "%"), NA)), 
              color = 'white', 
              position = position_stack(vjust = 0.5))+
    scale_fill_manual(name = "N TFs", values = unname(proj_cols("blue", "purple","teal","yellow", "pink")))+
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25), labels = formatC(seq(0, 1, 0.25), format = "g"))+
    labs(y = "Proportion", title = i, subtitle = annot)+
    theme_proj()+
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_markdown(size = 12, margin = margin(t = 0, unit = "in")),
          axis.title.y = element_text(size = 14, margin = margin(r = -0.05, unit = "in")),
          axis.title.x = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = 0, unit = "in"),
          legend.key.spacing.y = unit(0.01, "in"),
          legend.title = element_text(size = 12),
          panel.grid.major = element_blank(),
          plot.subtitle = element_markdown(size = 11, hjust = 1, face = "plain"),
          plot.title = element_markdown(size = 14))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(plots) + plot_layout(nrow = 1, guides = "collect")

cairo_pdf(file.path(res_path, paste0(figs, "_TF+.pdf")), height = 4, width = 9)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
plots <- list()

for(i in c("AD", "non-AD")){
  
  plot_data <- cbind(data@meta.data, t(data@assays$RNA@counts[c("CEBPA", "GPNMB", "GLDN"), ])) %>%
    mutate(dx = ifelse(donor %in% c(26, 33), "AD", "non-AD")) %>%
    filter(dx == i & cell_type == "Microglia") %>%
    mutate(CEBPA = ifelse(CEBPA > 0, "CEBPA+", "CEBPA-"),
           GPNMB = ifelse(GPNMB > 0, "GPNMB+", "GPNMB-"),
           GLDN = ifelse(GLDN > 0, "GLDN+", "GLDN-")) %>%
    rowwise() %>%
    mutate(stat = ifelse(CEBPA == "CEBPA+", paste0(CEBPA,"/", GPNMB, "/", GLDN), CEBPA))
  
  tmp <- plot_data %>% 
    filter(dx == i) %>% 
    mutate(scHPF_26_class = ifelse(scHPF_26_class == "Mic.26.high", 1, 0), 
           stat = factor(ifelse(stat=="CEBPA+/GPNMB+/GLDN+", 1, 0)))
  
  out <- glmer(stat ~ scHPF_26_class + scale(nCount_RNA)+ (1|orig.ident), data = tmp, family = "binomial")
  
  cc <- confint(out,parm="scHPF_26_class")
  coef <- exp(cbind(est=fixef(out)[2],cc))
  
  annot <- paste0("OR (95% CI) = ", format(coef[1], digits = 3, format = "g"), " (", 
                  format(coef[2], digits = 3, format = "g"), ", ", 
                  format(coef[3], digits = 3, format = "g"),
                  "), *p* = ",
                  formatC(summary(out)$coefficients[2,4], digits = 2, format = "e"))
  
  p <- plot_data %>%
    group_by(scHPF_26_class, stat) %>%
    summarise(n = n()) %>%
    group_by(scHPF_26_class) %>%
    summarise(stat, pct = n/sum(n)) %>%
    ggplot(aes(x = factor(scHPF_26_class, levels = c("Mic.26.low", "Mic.26.high"), labels = c("GPNMB<sup>Low</sup>", "GPNMB<sup>High</sup>")), y = pct, fill = stat)) +
    geom_col() +
    geom_text(aes(label = ifelse(pct > 0.1, paste0(round(pct*100, 1), "%"), NA)), 
              color = 'white', 
              position = position_stack(vjust = 0.5))+
    scale_fill_manual(values = unname(proj_cols("blue", "purple","teal","yellow", "pink")))+
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25), labels = formatC(seq(0, 1, 0.25), format = "g"))+
    guides(fill = guide_legend(ncol = 1)) +
    labs(y = "Proportion", title = i, subtitle = annot)+
    theme_proj()+
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_markdown(size = 12, margin = margin(t = 0, unit = "in")),
          axis.title.y = element_text(size = 14, margin = margin(r = -0.05, unit = "in")),
          axis.title.x = element_blank(),
          legend.key.width = unit(0.15, "in"),
          legend.key.height = unit(0.1, "in"),
          legend.text = element_text(size = 12, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = 0, unit = "in"),
          legend.key.spacing.y = unit(0.01, "in"),
          legend.key.spacing.x = unit(0.01, "in"),
          panel.grid.major = element_blank(),
          legend.title = element_blank(), 
          plot.subtitle = element_markdown(size = 11, hjust = 1, face = "plain"),
          plot.title = element_markdown(size = 14)) 
  
  plots <- c(plots, list(p))
  
}

p <- wrap_plots(plots) + plot_layout(nrow = 1, guides = "collect", axis_titles = "collect")

cairo_pdf(file.path(res_path, paste0(figs, "_CEBPA+GPNMB+GLDN+.pdf")), height = 4, width = 9)
print(p)
dev.off()
