#!/usr/bin/Rscript
#===============================================================================
# 
# MERSCOPE: 
# Figure S29.	Measures of local clustering for IRHigh factor (scHPF_20).
# Figure S33.	Measures of local clustering for GPNMBHigh factor (scHPF_26).
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); renv::load("~/projects/microglia_atlas/src/rlib")

library(colorspace)
library(ggrastr)
library(ggridges)
library(ggsignif)
library(ggtext)
library(patchwork)
library(Seurat)
library(tidyverse)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "merscope/2025-02-14_cellular_niches")
figs <-  c(20, 26)
names(figs) <- c("figs_s28_", "figs_s32_")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
run_metadata <- read_csv(file.path(home_path, "merscope/batch_names.csv"))
adjacent_cells <- readRDS(file.path(res_path, "adjacent_cells.rds"))
adjacent_counts <- readRDS(file.path(res_path, "adjacent_counts.rds"))

source(file.path(home_path, "src/00 Helpers/load_merscope.R"))

data@meta.data <- data@meta.data %>% 
  mutate(dx = ifelse(donor %in% c(26, 33), "AD", "non-AD"))

data@meta.data %>% group_by(dx, layer_niches) %>% summarise(mean(scHPF_20_class=="Mic.20.high", na.rm=T))
data@meta.data %>% group_by(layer_niches) %>% summarise(mean(scHPF_20_class=="Mic.20.high", na.rm=T))

data@meta.data %>% group_by(dx, layer_niches) %>% summarise(mean(scHPF_26_class=="Mic.26.high", na.rm=T))
data@meta.data %>% group_by(layer_niches) %>% summarise(mean(scHPF_26_class=="Mic.26.high", na.rm=T))


layer_labs <- c("WM", "GM.SLC17A7.GAS7.HSP90AB1", "GM.AHNAK.FOS.ACTB")
names(layer_labs) <- c("WM", "GM (L2-6)", "GM (L1)")
#===============================================================================
# 
# PANEL A
#
#===============================================================================
for(factor in c(20, 26)){
  
  cols <- proj_cols("teal", "pink")
  names(cols) <- c(paste0("Mic.", factor, ".low"), paste0("Mic.", factor, ".high"))
  
  plots <- list()
  for(niche in c("WM", "GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    
    if(!(factor==26 & niche=="GM.AHNAK.FOS.ACTB")){
      all_res <- data.frame()
      for(i in 1:12){
        file = file.path(res_path, paste0(run_metadata$batch_code[i], "_", niche, "_mglia_sfe.rds"))
        if(file.exists(file)){
          sfe <- readRDS(file)
          tmp <- cbind(sfe@colData, sfe@int_colData$localResults[["localG_perm"]]@listData[[paste0("scHPF_", factor)]]) %>% as_tibble() 
          all_res <- rbind(all_res, tmp)
        }
      }
      
      out <- cor.test(all_res[[paste0("scHPF_", factor)]], all_res$localG)
      p <- ifelse(out$p.value != 0, paste0(" <", formatC(out$p.value, digits = 2, format = "e")), " <0.001")
      lab <- paste0("r = ", formatC(out$estimate, digits = 2), ", p ", p)
      
      p <- ggplot(all_res, aes(x = !!sym(paste0("scHPF_", factor)), y = localG)) +
        geom_point_rast(aes(color = !!sym(paste0("scHPF_", factor, "_class"))), alpha=0.9, scale = 0.5) +
        geom_smooth(formula = y ~ x, method = "lm") +
        geom_hline(yintercept = mean(tmp[["localG"]]), lty = 2, color = proj_cols_grey("dark grey")) +
        geom_vline(xintercept = mean(tmp[[paste0("scHPF_", factor)]]), lty = 2, color = "gray") +
        geom_density2d() +
        labs(y = "Getis-Ord Local G<sub>i</sub>*", 
             x = paste0(names(factors)[factors %in% paste0("scHPF_", factor)], " score"),
             title = names(layer_labs)[layer_labs %in% niche])+
        scale_color_manual(values = cols, 
                           labels = if(factor == 20){c("IR<sup>High</sup>", "IR<sup>Low</sup>")} else {c("<i>GPNMB</i><sup>High</sup>", "<i>GPNMB</i><sup>Low</sup>")})+ 
        guides(color = guide_legend(override.aes = aes(size = 4)))+
        annotate("text", x = Inf, y = Inf, vjust = 1.2, hjust = 1, label= lab, family = "Helvetica")+
        theme_proj()+
        theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, "in"), 
              legend.title = element_blank(),
              legend.text = element_markdown(size = 12, margin = margin(l = 0.05, unit = "in")),
              legend.spacing.y = grid::unit(0.05, "in"),
              legend.direction = "horizontal",
              legend.position = "top",
              legend.margin = margin(t = 0, unit = "in"),
              panel.border = element_rect(linewidth = 1, color = "black"),
              panel.grid.major = element_blank(),
              plot.title = element_text(size = 14, face = "bold", family = "Helvetica"),
              axis.text = element_text(size = 12),
              axis.title.x = element_markdown(size = 12),
              axis.title.y = element_markdown(size = 12),
              axis.line.y = element_line(linewidth = .25, color = "black"),
              axis.line.x = element_line(linewidth = .25, color = "black"))
      
      plots = c(plots, list(p))
    }
  }
  
  p <- wrap_plots(plots) + plot_layout(nrow = 1, guides = "collect") &
    theme(legend.direction = "horizontal", legend.position = "bottom")
  
  pdf(file.path(res_path, paste0(names(figs)[figs %in% factor], "panelA.pdf")), height = 3, width = ifelse(factor == 26, 6, 9))
  print(p)
  dev.off()
  
}
#===============================================================================
# 
# PANEL B
#
#===============================================================================
cols <- proj_cols("pink", "teal")
names(cols) <- c("High", "Low")

for(factor_name in c(20, 26)){
  plots <- list()
  for(niche_name in c("WM","GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    if(!(factor_name==26 & niche_name=="GM.AHNAK.FOS.ACTB")){
      
      p <- adjacent_counts %>%
        filter(niche == niche_name & factor== factor_name) %>%
        group_by(stat, counts) %>%
        summarise(n=n()) %>%
        group_by(stat) %>%
        reframe(counts, n, pct = n/sum(n)) %>%
        mutate(niche = niche_name) %>%
        ggplot(aes(x = factor(counts), y = pct, fill = stat))+
        geom_col(color = "black", linewidth = 0.5, position = position_dodge(preserve = "single")) +
        labs(x = paste0("Number of adjacent ", ifelse(factor_name == 20, "IR", "<i>GPNMB</i>"), "<sup>High</sup> microglia"), 
             y = "Proportion of<br>layer-specific<br>microglia", 
             title = names(layer_labs)[layer_labs %in% niche_name])+
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25), limits = c(0, 1))+
        scale_fill_manual(values = cols, 
                          labels = if(factor_name == 20){c("IR<sup>High</sup>", "IR<sup>Low</sup>")} else {c("<i>GPNMB</i><sup>High</sup>", "<i>GPNMB</i><sup>Low</sup>")})+
        theme_proj() +
        theme(panel.grid.major = element_blank(),
              axis.title.x = element_markdown(size = 14),
              axis.title.y = element_markdown(size = 14),
              legend.title = element_blank(),
              legend.text = element_markdown(size = 12, margin = margin(l = 0.05, unit = "in")),
              legend.spacing.y = grid::unit(0.05, "in"),
              legend.key.width = grid::unit(0.2, "in"),
              legend.key.height = grid::unit(0.2, "in"),
              legend.margin = margin(t = 0, unit = "in"),
              plot.title = element_text(size = 14, face = "bold", family = "Helvetica"))
      
      plots <- c(plots, list(p))
    }
  }
  
  p <- wrap_plots(plots)+plot_layout(nrow=1, guides = "collect", axis_titles = "collect", axes = "collect") &
    theme(legend.direction = "horizontal",
          legend.position = "bottom")
  pdf(file.path(res_path, paste0(names(figs)[figs %in% factor_name], "panelB.pdf")), height = 3, width = ifelse(factor_name == 26, 6, 9))
  print(p)
  dev.off()
}

#===============================================================================
# 
# PANEL C
#
#===============================================================================
cols <- proj_cols("teal", "pink")
names(cols) <- c("Low", "High")

for(factor_name in c(20, 26)){
  plots <- list()
  for(niche_name in c("WM","GM.AHNAK.FOS.ACTB", "GM.SLC17A7.GAS7.HSP90AB1")){
    if(!(factor_name==26 & niche_name=="GM.AHNAK.FOS.ACTB")){
      
      
      tmp <- adjacent_counts %>% filter(niche == niche_name & factor == factor_name) 
      out <- fisher.test(ifelse(tmp$stat == "Low", 0, 1), ifelse(as.numeric(as.character(tmp$counts)) > 0, 1, 0))
      or <- unname(formatC(out$estimate, digits = 2, format = "g"))
      upper <- unname(format(out$conf.int[2], digits = 2, format = "g"))
      lower <- unname(format(out$conf.int[1], digits = 2, format = "g"))
      
      lab <- paste0("OR (95%CI) = ", or, " (", lower, ", ", upper, ")\np = ", formatC(out$p.value, digits = 2, format = "e"))
      mean <- adjacent_counts %>%
        filter(niche == niche_name & factor== factor_name) %>%
        group_by(stat) %>%
        summarise(pct = mean(counts>0)) %>%
        mutate(lab = paste0(round(pct*100, 2), "%"))
      
      p <- adjacent_counts %>%
        filter(niche == niche_name & factor== factor_name) %>%
        group_by(stat) %>%
        summarise(pct = mean(counts>0)) %>%
        ggplot(aes(x = stat, y = pct, fill = stat))+
        geom_col(color = "black", linewidth = 0.5, position = position_dodge(preserve = "single"))+
        geom_text(data = mean, aes(label = lab), vjust = 0)+
        labs(y = paste0("Proportion of microglia<br>sharing an ", ifelse(factor_name == 20, "IR", "<i>GPNMB</i>"), "<sup>High</sup><br>neighborhood"), 
             title = names(layer_labs)[layer_labs %in% niche_name])+
        annotate(geom="text", label = lab, hjust = 1, x = Inf, y = Inf, vjust = 1)+
        scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.125, 0.025), limits = c(0, 0.125))+
        scale_x_discrete(labels = if(factor_name == 20){c("IR<sup>High</sup>", "IR<sup>Low</sup>")} else {c("<i>GPNMB</i><sup>High</sup>", "<i>GPNMB</i><sup>Low</sup>")})+
        scale_fill_manual(values = cols)+
        theme_proj() +
        theme(panel.grid.major = element_blank(),
              legend.position = "none",
              axis.text.x =element_markdown(),
              plot.title = element_text(size = 14, face = "bold", family = "Helvetica"),
              axis.title.x = element_blank(),
              axis.title.y = element_markdown(size = 14))
      
      plots <- c(plots, list(p)) 
    }
  }
  
  p <- wrap_plots(plots)+plot_layout(nrow=1, axis_titles = "collect", axes = "collect")
  pdf(file.path(res_path, paste0(names(figs)[figs %in% factor_name], "panelC.pdf")), height = 3, width = ifelse(factor_name == 26, 6.5, 9))
  print(p)
  dev.off()
}