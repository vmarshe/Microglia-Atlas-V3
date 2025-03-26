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
library(Seurat)
library(tidyverse)
library(vegan)
library(VennDiagram)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sn")
fig <- "fig_s23_"

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))


#===============================================================================
#
# DATA
#
#===============================================================================

# ROSMAP METADATA
#-------------------------------------------------------------------------------
rosmap_meta <- "/path/to/dataset_707_basic_11-01-2023.txt"
rosmap_meta <- read_delim(rosmap_meta, delim="\t") %>% 
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  arrange(projid) %>%
  dplyr::select(projid, msex, age_death, educ, pmi, study) %>%
  distinct()

multiome_ids <- read_delim("dataset_basic_n3638.txt") %>%
  dplyr::select(projid, wgs_SampleID) %>%
  mutate(projid = str_pad(projid, 8, pad = "0")) %>%
  distinct()

non_overlapping_donors <- readRDS(file.path(home_path, "data/snuc_non_overlap_donors.rds"))

# CUIMC1
#-------------------------------------------------------------------------------
snuc <- readRDS(file.path(home_path, "proj_snuc/projection/snuc.rds"))
snuc@meta.data$cell_names <- rownames(snuc@meta.data)
snuc@meta.data <- snuc@meta.data %>% 
  left_join(snuc@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(snuc@meta.data) <- snuc@meta.data$cell_names

snuc@meta.data <- snuc@meta.data %>% mutate(projid = as.character(str_pad(projid, 8, pad = "0")))

for(i in c("projid", "msex", "pathoAD", "pmi", "age_death", "study")){
  
  tmp <- snuc@meta.data %>% dplyr::select(projid, cell_names) %>%
    left_join(rosmap_meta %>% dplyr::select(projid, all_of(i)))
  var <- tmp[,i]
  names(var) <- tmp$cell_names
  
  if(any(grepl(i, names(snuc@meta.data)))){ snuc@meta.data[,i] = NULL }
  
  snuc <- AddMetaData(snuc, var, i)
}
snuc_agg <- readRDS(file.path(home_path, "proj_snuc/projection/snuc_aggregated.rds"))%>%
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  left_join(rosmap_meta) 

# MIT
#-------------------------------------------------------------------------------
kellis <- readRDS(file.path(home_path, "proj_kellis/projection/kellis.rds"))
kellis@meta.data$cell_names <- rownames(kellis@meta.data)
kellis@meta.data <- kellis@meta.data %>% 
  left_join(kellis@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(kellis@meta.data) <- kellis@meta.data$cell_names

kellis@meta.data <- kellis@meta.data %>% mutate(projid = as.character(str_pad(projid, 8, pad = "0")))
rownames(kellis@meta.data) <- kellis@meta.data$cell_names

for(i in c("projid", "msex", "pathoAD", "pmi", "age_death", "study")){
  
  tmp <- kellis@meta.data %>% dplyr::select(projid, cell_names) %>%
    left_join(rosmap_meta %>% dplyr::select(projid, all_of(i)))
  var <- tmp[,i]
  names(var) <- tmp$cell_names
  
  if(any(grepl(i, names(kellis@meta.data)))){ kellis@meta.data[,i] = NULL }
  
  kellis <- AddMetaData(kellis, var, i)
}

keep_ids = kellis@meta.data %>% group_by(projid) %>% summarise(n=n()) %>% filter(n >= 10) %>% pull(projid)
kellis_agg <- kellis@meta.data %>% 
  filter(projid %in% keep_ids) %>%
  group_by(projid) %>% 
  summarise_at(all_of(unname(factors)), list(mean = mean)) %>%
  rename_at(paste0(factors, "_mean"), ~gsub("_mean","",.)) %>%
  left_join(rosmap_meta) 

kellis_ids <- non_overlapping_donors %>% filter(reference == "Kellis-Tsai") %>% pull(projid) #132
kellis <- subset(kellis, subset = orig.ident %in% kellis_ids) # 130

kellis_agg_overlap <- kellis_agg %>% filter(projid %in% snuc_agg$projid)
kellis_agg <- kellis_agg %>% filter(projid %in% kellis_ids)

# CUIMC2
#-------------------------------------------------------------------------------
multiome <- readRDS(file.path(home_path, "proj_multiome/projection/multiome.rds"))
multiome@meta.data$cell_names <- rownames(multiome@meta.data)
multiome@meta.data <- multiome@meta.data %>% 
  left_join(multiome@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(multiome@meta.data) <- multiome@meta.data$cell_names

multiome@meta.data <- multiome@meta.data %>% 
  left_join(multiome_ids, by = c("orig.ident"="wgs_SampleID")) %>%
  mutate(projid = as.character(str_pad(projid, 8, pad = "0")))

rownames(multiome@meta.data) <- multiome@meta.data$cell_names

for(i in c("projid", "msex", "pathoAD", "pmi", "age_death", "study")){
  
  tmp <- multiome@meta.data %>% dplyr::select(projid, cell_names) %>%
    left_join(rosmap_meta %>% dplyr::select(projid, all_of(i)))
  var <- tmp[,i]
  names(var) <- tmp$cell_names
  
  if(any(grepl(i, names(multiome@meta.data)))){ multiome@meta.data[,i] = NULL }
  
  multiome <- AddMetaData(multiome, var, i)
}
multiome_ids <- non_overlapping_donors %>% filter(reference == "Multiome") %>% pull(projid) #219
multiome <- subset(multiome, subset = projid %in% multiome_ids) #212

keep_ids <- multiome@meta.data %>% group_by(projid) %>% summarise(n=n()) %>% filter(n >= 10) %>% pull(projid)
multiome_agg <- multiome@meta.data %>% 
  filter(projid %in% keep_ids) %>%
  group_by(projid) %>% 
  summarise_at(all_of(unname(factors)), list(mean = mean)) %>%
  rename_at(paste0(factors, "_mean"), ~gsub("_mean", "", .)) %>%
  left_join(rosmap_meta) 

# Keep unique IDs
#-------------------------------------------------------------------------------
kellis_ids <- non_overlapping_donors %>% filter(reference == "Kellis-Tsai") %>% pull(projid) #132
multiome_ids <- non_overlapping_donors %>% filter(reference == "Multiome") %>% pull(projid) #219

kellis <- subset(kellis, subset = orig.ident %in% kellis_ids) # 130
multiome <- subset(multiome, subset = projid %in% multiome_ids) #212

# NORMALIZE
#-------------------------------------------------------------------------------
regs <- c("ARID5B", "CEBPA", "MITF", "PPARG")
snuc <- NormalizeData(snuc) %>% ScaleData(features = regs)
multiome <- NormalizeData(multiome) %>% ScaleData(features = regs)
kellis <- NormalizeData(kellis) %>% ScaleData(features = regs)

#===============================================================================
#
# PANEL A
#
#===============================================================================
pdf(file.path(res_path, paste0(fig, "panel_A.pdf")), height = 6, width = 6, onefile = T)
for(j in c("snuc", "multiome", "kellis")) { 
  
  data <- eval(parse(text = j))
  
  counts <- data@assays$RNA@scale.data[regs, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell_names")
  
  counts <- data@meta.data %>% 
    left_join(counts) %>%
    group_by(projid) %>% 
    summarise(factor = mean(scHPF_26), 
              ARID5B = mean(ARID5B), 
              umi = mean(nCount_RNA), 
              MITF = mean(MITF), 
              CEBPA = mean(CEBPA), 
              PPARG = mean(PPARG),
              pmi = unique(pmi),
              age = unique(age_death), 
              sex = unique(msex),
              study = unique(study))%>%
    drop_na()
  
  plot(varpart(counts$factor, counts$ARID5B, counts$CEBPA, counts$MITF, counts$PPARG),
       Xnames = c("ARID5B", "CEBPA", "MITF", "PPARG"),
       bg = proj_cols("blue", "yellow", "teal", "pink"), 
       alpha = 80,
       digits = 2,
       cex = 1,
       cutoff = 0.005)
  
}
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================
cols <- unname(proj_cols("blue", "yellow", "teal", "pink"))
names(cols) <- c("ARID5B", "CEBPA", "MITF","PPARG")

dataset_labs = c("snuc", "multiome", "kellis")
names(dataset_labs) = c("**CUIMC1** (*n*=424)", "**CUIMC2** (*n*=212)", "**MIT** (*n*=130)")

plots <- list()
for(j in c("snuc", "multiome", "kellis")) { 
  
  data <- eval(parse(text = j))
  
  counts <- data@assays$RNA@scale.data[regs, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell_names")
  
  counts <- data@meta.data %>% 
    left_join(counts) %>%
    group_by(projid) %>% 
    summarise(factor = mean(scHPF_26), 
              umi = mean(nCount_RNA), 
              ARID5B = mean(ARID5B), 
              MITF = mean(MITF), 
              CEBPA = mean(CEBPA), 
              PPARG = mean(PPARG),
              pmi = unique(pmi),
              age = unique(age_death), 
              sex = unique(msex),
              study = unique(study))%>%
    drop_na()
  
  res <- data.frame()
  for(i in regs){
    vars <- counts %>% dplyr::select(umi, pmi, age, sex, study)
    part <- varpart(counts$factor, counts[,i], vars)
    sum_part <- summary(part)
    
    sig <- anova.cca(rda(counts$factor, counts[,i], vars), permutations = 10000)
    res <- rbind(res, cbind(reg = i, contrib = sum_part$uniqpart[1], p = sig$`Pr(>F)`[1]))
  }
  
  res <- res %>% 
    mutate(contrib = as.numeric(as.character(contrib)), 
           p.adj = p.adjust(p, method = "BH"), 
           lab = ifelse(p.adj < 0.001, "***", ifelse(p.adj < 0.01, "**", ifelse(p.adj < 0.05, "*", "*"))),
           lab = paste0("Adj. *R*<sup>2</sup> = ", round(contrib, 2), lab))
  
  p <- counts %>%
    gather("reg", "count", ARID5B, CEBPA, MITF, PPARG) %>%
    mutate(dataset = names(dataset_labs)[dataset_labs %in% j]) %>%
    ggplot(aes(x = count, y = factor, color = reg)) + 
    geom_point(alpha = .7) +
    geom_smooth(method = "lm", color = "black") +
    geom_richtext(data = res, aes(x = Inf, y = Inf, label = lab), vjust = 1, hjust = 1, color = "black", fill = alpha("lightgrey", 0.1),
                  label.padding = unit(0.01, "in"), label.color = NA)+
    facet_grid(dataset~reg, scales = "free_x") + 
    scale_color_manual(values = cols)+
    labs(y = "Mean GPNMB<sup>high</sup> (26) program score", x = "Mean scaled TF expression level")+
    theme_proj()+
    theme(legend.position = "none",
          panel.spacing.x = unit(0.3, 'in'),
          panel.grid.major = element_blank(),
          axis.title.y = element_markdown(size = 14),
          strip.text = element_markdown(face = "bold", size = 12, margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in")))
  
  plots <- c(plots, list(p))
}

p <- wrap_plots(list(plots[[1]], 
                    plots[[2]]+ theme(strip.text.x = element_blank(), strip.background.x = element_blank()),
                    plots[[3]] + theme(strip.text.x = element_blank(), strip.background.x = element_blank()))) + 
  plot_layout(nrow = 3, axis_titles = "collect")

pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 9, width = 12)
print(p)
dev.off()
