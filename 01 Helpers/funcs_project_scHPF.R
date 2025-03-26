project_query = function(data, proj_data, ref_model){

  set.seed(1234)
  proj.umap <- umap_transform(X = proj_data, model = ref_model)
  colnames(proj.umap) <- paste0("schpfumap_", 1:2)

  data@reductions$schpf.umap <- CreateDimReducObject(embedding = proj.umap,
                                                    assay = "RNA",
                                                    key = "schpfumap_")

  return(data)
}

project_scHPF <- function(loom_file = NULL,
                         name = NULL, 
                         factors = NULL, 
                         cell_names = NULL, 
                         gene_names, 
                         key = "scHPF", 
                         assay = "RNA",
                         umap = NULL, 
                         cell_score_file = NULL, 
                         gene_score_file = NULL, 
                         res_path = NULL){

  require(logger)
  require(uwot)

  loom <- Connect(filename = loom_file,  mode = "r",  force = T)
  data <- as.Seurat(x = loom)
  loom$close_all()

  data@meta.data$cell_names <- Cells(data)

  cell_scores <- as.matrix(read_delim(cell_score_file, col_names = F))
  colnames(cell_scores) <- paste0("scHPF_", 1:ncol(cell_scores))
  rownames(cell_scores) <- cell_names
  cell_scores <- cell_scores[,factors]

  gene_scores <- as.matrix(read_delim(gene_score_file, col_names = F))
  colnames(gene_scores) <- paste0("scHPF_", 1:ncol(gene_scores))
  rownames(gene_scores) <- gene_names
  gene_scores <- gene_scores[,factors]

  data[[key]] <- CreateDimReducObject(embedding = cell_scores,
                                     loadings = gene_scores,
                                     assay = assay,
                                     key = paste0(key, "_"))

  data <- project_query(data = data,
                       proj_data = data[[key]]@cell.embeddings,
                       ref_model = umap)

  if(!dir.exists(file.path(res_path, name))) dir.create(file.path(res_path, name))
  saveRDS(data, file.path(res_path, name, paste0(name, ".rds")))

}
