#!/usr/bin/Rscript

#===============================================================================
# 
# DISCOVERY ARACNE MODEL
# This script gathers the output from ARACNe and calculates the TFMoA index (-1; 1).
# 
#===============================================================================

# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(readxl)
library(tidyverse)
library(viper)

# PATHS
#-------------------------------------------------------------------------------
home_path <- "~/projects/microglia_atlas"
res_path <- file.path(home_path, "aracne_sc")

# FUNCS
#-------------------------------------------------------------------------------
source(file.path(home_path, "src/00 Helpers/factors.R"))
source(file.path(home_path, "src/00 Helpers/theme.R"))

# DATA
#-------------------------------------------------------------------------------
aracne <- read_table(file.path(res_path, "network.txt.gz"))
regulon <- readRDS(file.path(res_path, "regulon.rds"))

#===============================================================================
#
# GET TFMODE
#
#===============================================================================
regulators <- unique(aracne$Regulator)

updated <- data.frame()

pb <- txtProgressBar(min = 1, max = length(regulators), style = 3)

for(i in 1:length(regulators)){
  
  setTxtProgressBar(pb, i)
  
  tfmode <- regulon[[regulators[i]]]$tfmode
  tfmode <- data.frame(Regulator = regulators[i], 
                      Target = names(tfmode),
                      tfmode = tfmode)
  
  likelihood = regulon[[regulators[i]]]$likelihood
  likelihood = data.frame(Regulator = regulators[i], 
                          Target = names(regulon[[regulators[i]]]$tfmode), 
                          likelihood = likelihood)
  
  tmp <- aracne %>% 
    filter(Regulator == regulators[i]) %>%
    left_join(tfmode) %>%
    left_join(likelihood)
  
  updated <- rbind(updated, tmp)
  
  rm(list = c("tfmode", "likelihood", "tmp"))
  
}
close(pb)

#===============================================================================
#
# CORRELATIONS
#
#===============================================================================
meta_mat <- read_delim(file.path(res_path, "metacell_matrix_norm_mean_nopc_names.norm.txt.gz"))
meta_mat_genes <- meta_mat$gene
meta_mat <- as.matrix(meta_mat[,2:ncol(meta_mat)])
rownames(meta_mat) <- meta_mat_genes

corr <- data.frame()

pb <- txtProgressBar(min = 1, max = length(regulators), style = 3)

for(i in 1:length(regulators)){
  
  setTxtProgressBar(pb, i)
  
  targets <- aracne %>% filter(Regulator == regulators[i]) %>% pull(Target)
  
  rho <- lapply(targets, function(x) cor.test(meta_mat[regulators[i],], meta_mat[x,], method = "spearman")$estimate)
  rho <- unlist(rho)
  
  tmp <- data.frame(regulator = regulators[i],
                   target = targets,
                   r = as.numeric(as.character(rho)))
  
  corr <- rbind(corr, tmp)
  
  rm(list = c("rho", "targets", "tmp"))
  
}
close(pb)


updated <- updated %>%
  mutate(up = ifelse(tfmode > 0, 1, 0),
         down = ifelse(tfmode < 0, 1, 0),
         dir = ifelse(up == 1, "up", ifelse(down == 1, "down", NA))) %>%
  left_join(corr, by = c("Regulator"="regulator", "Target"="target"))

write.csv(updated, file.path(res_path, "aracne_w_tfmode.csv"), row.names = F)