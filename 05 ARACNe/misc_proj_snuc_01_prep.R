
# LIBRARIES
#-------------------------------------------------------------------------------
renv::load("~/projects/microglia_atlas/src/rlib")

library(tidyverse)
library(Seurat)
library(SeuratDisk)

# PATHS
#-------------------------------------------------------------------------------
path <- "~/projects/microglia_atlas"
file_names <- c(snuc = "ROSMAP_singlenuc", 
                multiome = "multiome", 
                kellis = "kellis_mic")

# SCHPF GENES
#-------------------------------------------------------------------------------
genes <- read_delim(file.path(path, "data/scHPF_v3_genes.txt"), col_names = F, delim = "\t")$X1

# PROCESS DATA
#-------------------------------------------------------------------------------
for(i in c("kellis")){ #"snuc", "kellis", "multiome"
  
  res_path <- file.path(path, paste0("proj_", i, "_validation"))
  
  if(!dir.exists(file.path(res_path, "input"))) {dir.create(file.path(res_path, "input"))}
  
  # DATA
  #-----------------------------------------------------------------------------
  file_name <- file.path(path, paste0("proj_", i), "input", paste0(file_names[names(file_names) %in% i], ".loom"))
  data <- Connect(file_name, mode = "r", force = T)
  data <- as.Seurat(data)
  
  # IDs
  #-----------------------------------------------------------------------------
  if(i == "multiome"){
    multiome_ids <- read_delim(file.path(path, "datadataset_basic_n3638.txt")) %>%
      dplyr::select(projid, wgs_SampleID) %>%
      mutate(projid = str_pad(projid, 8, pad = "0")) %>%
      distinct()
    data@meta.data <- data@meta.data %>% 
      left_join(multiome_ids, by = c("orig.ident"="wgs_SampleID")) %>%
      mutate(projid = as.character(str_pad(projid, 8, pad = "0")))
    
    data@meta.data$orig.ident <- data@meta.data$projid
    
  } else {
    data@meta.data$orig.ident <- as.character(str_pad(data@meta.data$projid, 8, pad = "0"))
  }
  
  # MISSING GENES
  #-----------------------------------------------------------------------------
  length(genes[!genes %in% rownames(data)])
  
  # REMOVE GENES
  #-----------------------------------------------------------------------------
  counts <- data[["RNA"]]@counts[genes, ]
  counts["ARID5B", ] <- 0
  counts["CEBPA", ] <- 0
  counts["PPARG", ] <- 0
  counts["MITF", ] <- 0
  counts["IRF7", ] <- 0
  
  # SAVE
  #-----------------------------------------------------------------------------
  tmp <- CreateSeuratObject(counts[genes, ])
  tmp@meta.data$orig.ident <- data@meta.data$orig.ident
  
  if(i == "kellis"){
    tmp <- subset(tmp, subset = orig.ident != "10102206") # donor has only one cell 
  }
  
  SaveLoom(tmp, filename = file.path(res_path, "input/data.loom"))
  write_delim(data.frame(Cells(tmp)), file.path(res_path, "input/cells.txt"), col_names = F)
  
  rm(list = c("res_path", "data", "counts", "tmp")); gc()
  
}