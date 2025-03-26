#!/usr/bin/Rscript

#===============================================================================
#
# MICROGLIA ATLAS V3
#
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
library(renv); load("~/projects/microglia_atlas/src/rlib")

library(circlize)
library(ComplexHeatmap)
library(forestplot)
library(ggrastr)
library(patchwork)
library(rms)
library(schex)
library(Seurat)
library(tidyverse)
library(writexl)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "proj_snuc/analysis")
fig <- "fig_3_"

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
rosmap_meta <- read_delim(file.path(home_path, "data/dataset_707_basic_11-01-2023.txt"), delim="\t") %>% 
  mutate(projid = as.character(str_pad(projid, 8, pad = "0"))) %>%
  arrange(projid) %>%
  dplyr::select(projid, msex, age_death, educ, pmi, study, amyloid, tangles, cogng_demog_slope, niareagansc) %>%
  distinct() %>%
  mutate(pathoAD = factor(ifelse(niareagansc %in% c(3, 4), 0, 1), levels = c(0, 1)),
         amyloid_sqrt = sqrt(amyloid),
         tangles_sqrt = sqrt(tangles)) 

multiome_ids <- read_delim(file.path(home_path, "data/dataset_basic_n3638.txt")) %>%
  dplyr::select(projid, wgs_SampleID) %>%
  mutate(projid = str_pad(projid, 8, pad = "0")) %>%
  distinct()

non_overlapping_donors <- readRDS(file.path(home_path, "data/snuc_non_overlap_donors.rds"))

# SNUC
#-------------------------------------------------------------------------------
snuc <- readRDS(file.path(home_path, "proj_snuc/projection/snuc.rds"))
snuc@meta.data$cell_names <- rownames(snuc@meta.data)
snuc@meta.data <- snuc@meta.data %>% 
  left_join(snuc@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(snuc@meta.data) <- snuc@meta.data$cell_names

snuc@meta.data <- snuc@meta.data %>% mutate(projid = as.character(str_pad(projid, 8, pad = "0")))

for(i in c("projid", "msex", "pathoAD", "amyloid", "tangles", "cogng_demog_slope")){
  
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

# KELLIS
#-------------------------------------------------------------------------------
kellis <- readRDS(file.path(home_path, "proj_kellis/projection/kellis.rds"))
kellis@meta.data$cell_names <- rownames(kellis@meta.data)
kellis@meta.data <- kellis@meta.data %>% 
  left_join(kellis@reductions$scHPF@cell.embeddings %>% as.data.frame() %>% rownames_to_column("cell_names"))
rownames(kellis@meta.data) <- kellis@meta.data$cell_names

kellis@meta.data <- kellis@meta.data %>% mutate(projid = as.character(str_pad(projid, 8, pad = "0")))
rownames(kellis@meta.data) <- kellis@meta.data$cell_names

for(i in c("projid", "msex", "pathoAD", "amyloid", "tangles", "cogng_demog_slope")){
  
  tmp <- kellis@meta.data %>% dplyr::select(projid, cell_names) %>%
    left_join(rosmap_meta %>% dplyr::select(projid, all_of(i)))
  var <- tmp[,i]
  names(var) <- tmp$cell_names
  
  if(any(grepl(i, names(kellis@meta.data)))){ kellis@meta.data[,i] = NULL }
  
  kellis <- AddMetaData(kellis, var, i)
}

keep_ids <- kellis@meta.data %>% group_by(projid) %>% summarise(n=n()) %>% filter(n >= 10) %>% pull(projid)
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

# MULTIOME
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

for(i in c("projid", "msex", "pathoAD", "amyloid", "tangles", "cogng_demog_slope")){
  
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

#===============================================================================
#
# PANEL A
#
#===============================================================================
ref_hex <- readRDS(file.path(home_path, "data/ref_hex.rds"))

hex <- make_hexbin(snuc, nbins = 100, dimension_reduction = "schpf.umap")

p <- ggplot() +
  geom_hex(data = ref_hex, aes(x = x, y = y), fill = "#F0F0F0", stat = "identity") +
  geom_hex(data = data.frame(hex@misc$hexbin$hexbin.matrix), aes(x = x, y = y, fill = number_of_cells), stat = "identity") +
  scale_fill_viridis_c(name = "N nuclei\nper hexbin")+
  theme_umap()+
  theme(legend.key.height = unit(0.25, "in"),
        legend.key.width = unit(0.25, "in"),
        legend.text = element_text(size = 12, margin = margin(0, 0, 0, 0.05, 'in'), family = "Helvetica"),
        legend.title = element_text(size = 14, margin = margin(b = 0.1, unit = "in")),
        legend.margin = margin(l = -0.2, unit = 'in'),
        legend.ticks.length = unit(0.05, "in"),
        legend.ticks = element_line(linewidth = 0.5, color = "white"),
        plot.margin = margin(0, 0, 0, 0, "in"))

pdf(file.path(home_path, "proj_snuc", "analysis", paste0(fig, "_panelA_schex.pdf")), height = 3.5, width = 4)
print(p)
dev.off()

#===============================================================================
#
# PANEL B
#
#===============================================================================

# MODELS - CUIMC1
#-------------------------------------------------------------------------------
snuc_cont <- data.frame()

for(trait in c("amyloid_sqrt", "tangles_sqrt", "cogng_demog_slope")){
  for(i in factors){
    
    tmp <- snuc_agg %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      drop_na()
    
    d <- datadist(tmp); options(datadist = "d")
    if(trait == "cogng_demog_slope"){
      mod <- ols(trait ~ study + age_death + msex + educ + factor, data = tmp)
    } else {
      mod <- ols(trait ~ study + age_death + msex + pmi + factor, data = tmp)
    }
    
    snuc_cont <- rbind(snuc_cont, 
                       cbind(dataset = "snuc", 
                             factor = i, 
                             trait, 
                             n = nrow(tmp),
                             mod_p = unname(1-pchisq(mod$stats[2], mod$stats[3])),
                             coef = coefficients(summary.lm(mod))["factor","Estimate"], 
                             se = coefficients(summary.lm(mod))["factor","Std. Error"],
                             lower = confint(mod, 'factor', level=0.95)[1],
                             upper = confint(mod, 'factor', level=0.95)[2],
                             p = summary.lm(mod)$coefficients["factor", "Pr(>|t|)"]))
    
    rm(list = c("tmp", "d", "mod"))
  }
}

snuc_cat <- data.frame()
for(trait in c("pathoAD")){
  for(i in factors){
    
    tmp <- snuc_agg %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      mutate(trait = factor(trait)) %>%
      drop_na()
    
    d <- datadist(tmp); options(datadist = "d")
    mod <- lrm(trait ~ study + age_death + msex + pmi + factor, data = tmp)
    
    coefs <- summary(mod) %>% as.data.frame() %>% rownames_to_column("var") %>% filter(var=="factor")
    
    snuc_cat <- rbind(snuc_cat, 
                      cbind(dataset = "snuc", 
                            factor = i, 
                            trait, 
                            n = nrow(tmp),
                            mod_p = mod$stats["P"],
                            coef = coefs$Effect, 
                            se = coefs$S.E.,
                            lower = coefs$`Lower 0.95`,
                            upper = coefs$`Upper 0.95`,
                            p = anova(mod)["factor", "P"]))
    
    rm(list = c("tmp", "d", "mod", "coefs"))
  }
}

# MODELS - CUIMC2
#-------------------------------------------------------------------------------
multiome_cont <- data.frame()
for(trait in c("amyloid_sqrt", "tangles_sqrt", "cogng_demog_slope")){
  for(i in factors){
    
    tmp <- multiome_agg %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      drop_na()
    
    d <- datadist(tmp); options(datadist = "d")
    mod <- ols(trait ~ study + age_death + msex + pmi + factor, data = tmp)
    
    multiome_cont <- rbind(multiome_cont, 
                           cbind(dataset = "multiome", 
                                 factor = i, 
                                 trait, 
                                 n = nrow(tmp),
                                 mod_p = unname(1-pchisq(mod$stats[2], mod$stats[3])),
                                 coef = coefficients(summary.lm(mod))["factor","Estimate"], 
                                 se = coefficients(summary.lm(mod))["factor","Std. Error"],
                                 lower = confint(mod, 'factor', level=0.95)[1],
                                 upper = confint(mod, 'factor', level=0.95)[2],
                                 p = summary.lm(mod)$coefficients["factor", "Pr(>|t|)"]))
    
    rm(list = c("tmp", "d", "mod"))
  }
}

multiome_cat <- data.frame()
for(trait in c("pathoAD")){
  for(i in factors){
    
    tmp <- multiome_agg %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      mutate(trait = factor(trait)) %>%
      drop_na()
    
    d <- datadist(tmp); options(datadist = "d")
    mod <- lrm(trait ~ study + age_death + msex + pmi + factor, data = tmp)
    
    coefs <- summary(mod) %>% as.data.frame() %>% rownames_to_column("var") %>% filter(var=="factor")
    
    multiome_cat <- rbind(multiome_cat, 
                          cbind(dataset = "multiome", 
                                factor = i, 
                                trait, 
                                n = nrow(tmp),
                                mod_p = mod$stats["P"],
                                coef = coefs$Effect, 
                                se = coefs$S.E.,
                                lower = coefs$`Lower 0.95`,
                                upper = coefs$`Upper 0.95`,
                                p = anova(mod)["factor", "P"]))
    
    rm(list = c("tmp", "d", "mod", "coefs"))
  }
}

# MODELS - MIT
#-------------------------------------------------------------------------------
kellis_cont <- data.frame()
for(trait in c("amyloid_sqrt", "tangles_sqrt", "cogng_demog_slope")){
  for(i in factors){
    
    tmp <- kellis_agg %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      drop_na()
    
    d <- datadist(tmp); options(datadist = "d")
    mod <- ols(trait ~ study + age_death + msex + pmi + factor, data = tmp)
    
    kellis_cont <- rbind(kellis_cont, 
                         cbind(dataset = "kellis", 
                               factor = i, 
                               trait, 
                               n = nrow(tmp),
                               mod_p = unname(1-pchisq(mod$stats[2], mod$stats[3])),
                               coef = coefficients(summary.lm(mod))["factor","Estimate"], 
                               se = coefficients(summary.lm(mod))["factor","Std. Error"],
                               lower = confint(mod, 'factor', level=0.95)[1],
                               upper = confint(mod, 'factor', level=0.95)[2],
                               p = summary.lm(mod)$coefficients["factor", "Pr(>|t|)"]))
    
    rm(list = c("tmp", "d", "mod"))
  }
}

kellis_cat <- data.frame()
for(trait in c("pathoAD")){
  for(i in factors){
    
    tmp <- kellis_agg %>%
      filter(!projid %in% intersect(snuc$projid, kellis$projid)) %>%
      rename_at(i, ~"factor") %>%
      left_join(rosmap_meta %>% dplyr::select(projid, trait = all_of(trait), age_death, msex, pmi, study, educ)) %>%
      mutate(trait = factor(trait)) %>%
      drop_na()
    
    d <- datadist(tmp); options(datadist = "d")
    mod <- lrm(trait ~ study + age_death + msex + pmi + factor, data = tmp)
    
    coefs <- summary(mod) %>% as.data.frame() %>% rownames_to_column("var") %>% filter(var=="factor")
    
    kellis_cat <- rbind(kellis_cat, 
                        cbind(dataset = "kellis", 
                              factor = i, 
                              trait, 
                              n = nrow(tmp),
                              mod_p = mod$stats["P"],
                              coef = coefs$Effect, 
                              se = coefs$S.E.,
                              lower = coefs$`Lower 0.95`,
                              upper = coefs$`Upper 0.95`,
                              p = anova(mod)["factor", "P"]))
    
    rm(list = c("tmp", "d", "mod", "coefs"))
  }
}

rbind(snuc_cont, snuc_cat, multiome_cont, multiome_cat, kellis_cont, kellis_cat) %>%
  mutate_at(all_of(c("mod_p", "coef", "se", "p", "n", "lower", "upper")), ~as.numeric(as.character(.))) %>%
  group_by(dataset, trait) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  write_xlsx(file.path(home_path, "proj_snuc", "analysis", "individual_mods.xlsx"))

# META-ANALYSIS (KELLIS & MULTIOME)
#-------------------------------------------------------------------------------
cont <- rbind(kellis_cont, multiome_cont) %>%
  mutate_at(all_of(c("mod_p", "coef", "se", "p", "n", "lower", "upper")), ~as.numeric(as.character(.)))

cont_meta_two <- data.frame()
for(i in factors){
  for(j in c("amyloid_sqrt", "tangles_sqrt", "cogng_demog_slope")){
    tmp <- cont %>% filter(factor == i & trait == j)
    meta <- rma(yi = tmp$coef, sei = tmp$se, weights = tmp$n, data = tmp, method = "REML")
    
    cont_meta_two <- rbind(cont_meta_two, cbind(factor = i, 
                                                trait = j,
                                                meta.beta = c(meta$beta), 
                                                meta.se = meta$se, 
                                                meta.ci.lb = c(meta$ci.lb),
                                                meta.ci.ub = c(meta$ci.ub),
                                                meta.z = meta$zval,
                                                meta.p = c(meta$pval)))
  }
}

cont_meta_two <- cont_meta_two %>%
  mutate_at(all_of(grep("meta", names(cont_meta_two), v=T)), ~as.numeric(as.character(.)))

cat <- rbind(kellis_cat, multiome_cat) %>%
  mutate_at(all_of(c("mod_p", "coef", "se", "p", "n", "lower", "upper")), ~as.numeric(as.character(.)))

cat_meta_two <- data.frame()
for(i in factors){
  for(j in c("pathoAD")){
    tmp <- cat %>% filter(factor == i & trait == j)
    meta <- rma(yi = tmp$coef, sei = tmp$se, weights = tmp$n, data = tmp, method = "REML")
    
    cat_meta_two <- rbind(cat_meta_two, cbind(factor = i, 
                                              trait = j,
                                              meta.beta = c(meta$beta), 
                                              meta.se = meta$se, 
                                              meta.ci.lb = c(meta$ci.lb),
                                              meta.ci.ub = c(meta$ci.ub),
                                              meta.z = meta$zval,
                                              meta.p = c(meta$pval)))
  }
}

cat_meta_two <- cat_meta_two %>%
  mutate_at(all_of(grep("meta", names(cont_meta_two), v=T)), ~as.numeric(as.character(.)))

# META-ANALYSIS
#-------------------------------------------------------------------------------
cont <- rbind(snuc_cont, kellis_cont, multiome_cont) %>%
  mutate_at(all_of(c("mod_p", "coef", "se", "p", "n", "lower", "upper")), ~as.numeric(as.character(.)))

cont_meta_three <- data.frame()
for(i in factors){
  for(j in c("amyloid_sqrt", "tangles_sqrt", "cogng_demog_slope")){
    tmp <- cont %>% filter(factor == i & trait == j)
    meta <- rma(yi = tmp$coef, sei = tmp$se, weights = tmp$n, data = tmp, method = "REML")
    
    cont_meta_three <- rbind(cont_meta_three, cbind(factor = i, 
                                                    trait = j,
                                                    meta.beta = c(meta$beta), 
                                                    meta.se = meta$se, 
                                                    meta.ci.lb = c(meta$ci.lb),
                                                    meta.ci.ub = c(meta$ci.ub),
                                                    meta.z = meta$zval,
                                                    meta.p = c(meta$pval)))
  }
}

cont_meta_three <- cont_meta_three %>%
  mutate_at(all_of(grep("meta", names(cont_meta_two), v=T)), ~as.numeric(as.character(.)))

cat <- rbind(snuc_cat, kellis_cat, multiome_cat) %>%
  mutate_at(all_of(c("mod_p", "coef", "se", "p", "n", "lower", "upper")), ~as.numeric(as.character(.)))

cat_meta_three <- data.frame()
for(i in factors){
  for(j in c("pathoAD")){
    tmp <- cat %>% filter(factor == i & trait == j)
    meta <- rma(yi = tmp$coef, sei = tmp$se, weights = tmp$n, data = tmp, method = "REML")
    
    cat_meta_three <- rbind(cat_meta_three, cbind(factor = i, 
                                                  trait = j,
                                                  meta.beta = c(meta$beta), 
                                                  meta.se = meta$se, 
                                                  meta.ci.lb = c(meta$ci.lb),
                                                  meta.ci.ub = c(meta$ci.ub),
                                                  meta.z = meta$zval,
                                                  meta.p = c(meta$pval)))
  }
}

cat_meta_three <- cat_meta_three %>%
  mutate_at(all_of(grep("meta", names(cont_meta_two), v=T)), ~as.numeric(as.character(.)))

rbind(cont_meta_two, cat_meta_two) %>%
  group_by(trait) %>%
  mutate(meta.p.adj = p.adjust(meta.p, method = "BH")) %>%
  write_xlsx(file.path(home_path, "proj_snuc", "analysis", "meta_mods_two.xlsx"))

rbind(cont_meta_three, cat_meta_three) %>%
  group_by(trait) %>%
  mutate(meta.p.adj = p.adjust(meta.p, method = "BH")) %>%
  write_xlsx(file.path(home_path, "proj_snuc", "analysis", "meta_mods_three.xlsx"))

# PLOTS
#-------------------------------------------------------------------------------
mat <- data.frame()
for(i in c("snuc","kellis", "multiome")){
  
  tmp <- rbind(eval(parse(text = paste0(i, "_cont"))), eval(parse(text = paste0(i, "_cat")))) %>%
    mutate_at(all_of(c("coef", "se", "p")), ~as.numeric(as.character(.))) %>%
    group_by(trait) %>%
    mutate(p.adj = p.adjust(p, method = "BH"),
           log10.padj = ifelse(coef > 0, -log10(p.adj), log10(p.adj))) %>%
    dplyr::select(factor, trait, log10.padj) %>%
    mutate(dataset = i)
  
  mat <- rbind(mat, tmp)
}

tmp_two <- rbind(cat_meta_two, cont_meta_two) %>%
  group_by(trait) %>%
  mutate(p.adj = p.adjust(meta.p, method = "BH"),
         log10.padj = ifelse(meta.beta > 0, -log10(p.adj), log10(p.adj))) %>%
  dplyr::select(factor, trait, log10.padj) %>%
  mutate(dataset = "meta_two")

tmp_three <- rbind(cat_meta_three, cont_meta_three) %>%
  group_by(trait) %>%
  mutate(p.adj = p.adjust(meta.p, method = "BH"),
         log10.padj = ifelse(meta.beta > 0, -log10(p.adj), log10(p.adj))) %>%
  dplyr::select(factor, trait, log10.padj) %>%
  mutate(dataset = "meta_three")

labs <- c("Pathological AD", "Amyloid", "Tangles", "Cognitive decline")

mat <- rbind(mat, tmp_two) %>%
  rbind(tmp_three) %>%
  mutate(factor = factor(factor, levels = factors, labels = names(factors)),
         trait = factor(trait, 
                        levels = c("pathoAD", "amyloid_sqrt", "tangles_sqrt",  "cogng_demog_slope"),
                        labels = labs)) %>%
  mutate(trait = paste0(dataset, "_", trait)) %>%
  dplyr::select(factor, trait, log10.padj) %>%
  spread(factor, log10.padj) %>%
  column_to_rownames("trait") %>%
  as.matrix()

mat <- mat[c(paste0("snuc_", labs), paste0("multiome_", labs), paste0("kellis_", labs), paste0("meta_two_", labs), paste0("meta_three_", labs)), ]
row_split <- str_extract(rownames(mat), "snuc|multiome|kellis|meta_two|meta_three")
row_split <- factor(row_split, levels = c("snuc","multiome", "kellis", "meta_two", "meta_three"), 
                    labels = c("ROS-MAP\nCUIMC1", "ROS-MAP\nCUIMC2", 
                               "ROS-MAP\nMIT", 
                               "Meta-analysis\n(CUIMC2+MIT)", "Meta-analysis\n(all three)"))
rownames(mat) <- gsub("snuc_|multiome_|kellis_|meta_two_|meta_three_", "", rownames(mat))

cols <- c(proj_cols("yellow", "blue", "pink"),proj_cols_grey("dark grey"), "black")


cairo_pdf(file.path(res_path, paste0(fig, "panel_B.pdf")), height = 7.5, width = 9)
set.seed(1234)
p <- Heatmap(mat, name = "mat",
             col = colorRamp2(breaks = c(min(mat, na.rm=T), 0, max(mat, na.rm=T)),
                              colors = c(unname(proj_cols("blue")),"white", unname(proj_cols("red")))),
             cell_fun = function(j, i, x, y, w, h, fill) {
               if(!is.na(mat[i, j]) & abs(mat[i, j]) > -log10(0.05)) {
                 grid.points(pch = 8, x, y, size = unit(0.065, "in"))
               }
             },
             rect_gp = gpar(col = "white", lwd = 1.5),
             row_names_side = "left",
             row_names_gp = gpar(fontfamily = "Helvetica", fontsize = 12),
             row_title_side = "right",
             row_title_rot = 270,
             row_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 9.5),
             row_gap = unit(0.1, "in"),
             right_annotation = rowAnnotation(foo = anno_empty(border = FALSE, width = unit(0.3, "in"))),
             row_split = row_split,
             cluster_rows = F,
             cluster_columns = T,
             column_names_rot = 45,
             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 12),
             column_title_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = 0),
             heatmap_legend_param = list(
               title = expression('Coefficient-signed -log'[10]*'(FDR '*italic(p)*'-value)'),
               legend_width = unit(1.5, "in"),
               direction = "horizontal",
               title_position = "lefttop",
               border = "black",
               title_gp = gpar(fontfamily = "Helvetica", fontsize = 12),
               labels_gp = gpar(fontfamily = "Helvetica", fontsize = 10)
             ))

draw(p,
     align_heatmap_legend = "heatmap_center",
     heatmap_legend_side = "top",
     background = "transparent", 
     padding = unit(c(0.2, 0.2, 0.2, 0.2), "in"))
for(i in 1:5){
  decorate_heatmap_body("mat", slice = i, { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2)) })
  decorate_annotation("foo", slice = i, { grid.rect(x = 0, width = unit(0.1, "in"), gp = gpar(fill = cols[i], col = NA), just = "left") })
}

dev.off()

#===============================================================================
#
# PANEL C
#
#===============================================================================
all_meta <- rbind(cat_meta_three, cont_meta_three) %>% 
  group_by(trait) %>%
  mutate(meta.p.adj = p.adjust(meta.p, method = "BH")) %>%
  mutate(p.adj = p.adjust(meta.p, method = "BH"))

pdf(file.path(res_path, paste0(fig, "panel_C")), height = 2, width = 3.5, onefile = T)
for(j in "pathoAD"){
  for(i in paste0("scHPF_", c(26))){
    
    tmp <- all_meta %>% filter(trait == j & factor == i)
    tmp2 <- cat %>% filter(trait == j & factor ==i)
    
    base_data <- tibble::tibble(mean = exp(tmp2$coef), lower = exp(tmp2$lower),  upper = exp(tmp2$upper), study = tmp2$dataset) %>%
      mutate(study = recode(study, "snuc"="CUIMC1", "kellis"="MIT", "multiome"="CUIMC2"),
             study = factor(study, levels = c("CUIMC1", "CUIMC2", "MIT"))) %>%
      arrange(study)
    
    sum_col <- ifelse(tmp$p.adj < 0.05 & tmp$meta.beta > 0, proj_cols("red"), ifelse(tmp$p.adj < 0.05 & tmp$meta.beta < 0, proj_cols("blue"), proj_cols_grey("med grey")))
    p <- base_data |>
      forestplot(labeltext = c(study), xlab = "OR (95% CI)") |>
      fp_set_style(box = "#282828", 
                   line = "#282828", 
                   summary = sum_col, 
                   txt_gp = fpTxtGp(ticks = gpar(fontfamily = "Helvetica", cex = 1),
                                    label = gpar(fontfamily = "Helvetica", cex = 0.9),
                                    xlab = gpar(fontfamily = "Helvetica", cex = 0.9))) |> 
      fp_add_header(study = "") |>
      fp_append_row(mean  = exp(tmp$meta.beta), lower = exp(tmp$meta.ci.lb), upper = exp(tmp$meta.ci.ub), 
                    study = "Meta-analysis", is.summary = TRUE) |> 
      fp_decorate_graph(grid = structure(1, gp = gpar(lty = 1, col = "#282828")))
    
    print(p)
  }
}

for(j in "amyloid_sqrt"){
  for(i in paste0("scHPF_", c(26))){
    
    tmp <- all_meta %>% filter(trait == j & factor == i)
    tmp2 <- cont %>% filter(trait == j & factor ==i)
    
    base_data <- tibble::tibble(mean  = tmp2$coef, lower = tmp2$lower, upper = tmp2$upper, study = tmp2$dataset) %>%
      mutate(study = recode(study, "snuc"="CUIMC1", "kellis"="MIT", "multiome"="CUIMC2"),
             study = factor(study, levels = c("CUIMC1", "CUIMC2", "MIT"))) %>%
      arrange(study)
    
    sum_col <- ifelse(tmp$p.adj < 0.05 & tmp$meta.beta > 0, proj_cols("red"), ifelse(tmp$p.adj < 0.05 & tmp$meta.beta < 0, proj_cols("blue"), proj_cols_grey("med grey")))
    p <- base_data |>
      forestplot(labeltext = c(study), xlab = "Beta (95% CI)") |>
      fp_set_style(box = "#282828", line = "#282828", summary = sum_col, 
                   txt_gp = fpTxtGp(ticks = gpar(fontfamily = "Helvetica", cex = 1),
                                    label = gpar(fontfamily = "Helvetica", cex = 0.9),
                                    xlab = gpar(fontfamily = "Helvetica", cex = 0.9))) |> 
      fp_add_header(study = "") |>
      fp_append_row(mean  = tmp$meta.beta, lower = tmp$meta.ci.lb, upper = tmp$meta.ci.ub, study = "Meta-analysis", is.summary = TRUE) |> 
      fp_decorate_graph(grid = structure(0, gp = gpar(lty = 1, col = "#282828")))
    
    print(p)
  }
}
dev.off()