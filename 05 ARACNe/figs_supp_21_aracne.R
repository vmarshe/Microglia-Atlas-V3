#!/usr/bin/Rscript

#===============================================================================
# 
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(ggtext)
library(patchwork)
library(readxl)
library(Seurat)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sc")
fig <- "fig_s21_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
aracne <- read_csv(file.path(home_path, "aracne_sc/aracne_w_tfmode.csv"))
regulators <- unique(aracne$Regulator)
schpf <- readRDS(file.path(home_path, "data/scHPF_v3_reduc_FULL.rds"))

#===============================================================================
#
# PANEL A - TFs & FACTORS
#
#===============================================================================
results <- data.frame()

for(i in factors){
  
  genes <- names(sort(schpf@feature.loadings[,i], decreasing = T)[1:100])
  ind <- which(genes %in% regulators)
  regs_in_top <- genes[genes %in% regulators]
  
  if(length(regs_in_top) > 0){ results <- rbind(results, cbind(factor = i, regs = regs_in_top, ind)) }
}

p <- results %>%
  group_by(factor) %>% 
  arrange(ind) %>%
  summarise(n = n(), labs = paste0(na.omit(unique(regs)[1:3]), collapse = "/")) %>%
  mutate(factor = factor(factor, levels = factors, labels = names(factors))) %>%
  arrange(n) %>%
  mutate(factor = factor(factor, levels = factor)) %>%
  ggplot(aes(x = n, y = factor)) + 
  geom_col(fill = proj_cols("yellow")) +
  geom_text(aes(label = labs, 
                x = ifelse(n < 11, n+0.5, n-0.5), 
                hjust = ifelse(n < 11, 0, 1)), 
            color = "black", size = 2.75, family = "Helvetica") +
  scale_x_continuous(expand = c(0, 0))+
  labs(x = "*N* TFs")+
  theme_proj()+
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "#F0F0F0"),
        strip.text = element_text(size = 11, family = "Helvetica"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(r = 0, unit = "in")),
        axis.title.x = element_markdown(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        plot.margin = margin(t =0.1, r=0.1, b=0.1, l=0.1, unit = "in"),
        panel.grid.major = element_blank())

cairo_pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 4, width = 4.25)
print(p)
dev.off()
#===============================================================================
#
# PANEL C - SHARED ENRICHED TFs
#
#===============================================================================
plot_data <- read_xlsx(file.path(res_path, "tab_aracne_scHPF_enrichment.xlsx")) %>% 
  filter(qvalue < 0.05 & p.adj.bonf < 0.05) %>% 
  group_by(signature) %>% 
  arrange(p.adj.bonf) %>%
  summarise(n = n(), factors = paste0(ID, collapse = ", ")) %>%
  arrange(desc(n)) %>%
  top_n(20, wt = n) %>%
  mutate(signature = factor(signature, levels = rev(signature))) 

p <- ggplot(plot_data[1:20, ],aes(x = n, y = signature)) +
  geom_col(fill = proj_cols("teal"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 10, 2), labels = formatC(seq(0, 10, 2), format="g"))+
  labs(x = "*N* shared expression programs", y = "Shared TFs")+
  theme_proj() +
  theme(legend.position = "none", 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(r = 0, unit = "in")),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_markdown(size = 10),
        plot.margin = margin(t =0.1, r=0.2, b=0.1, l=0.1, unit = "in"),
        panel.grid.major = element_blank())

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 4, width = 3)
print(p)
dev.off()
