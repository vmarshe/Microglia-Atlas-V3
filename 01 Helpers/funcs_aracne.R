#===============================================================================
# 
# Courtesy of John Tuddenham with modifications.
#
#===============================================================================
make_metacells <- function(data = NULL,
                           dims = 1:40,
                           reduction = NULL,
                           k.param = 50,
                           features = NULL,
                           res_path = getwd(),
                           seed = 1234){
  
  require(logger)
  start <- Sys.time()
  
  data <- FindNeighbors(object = data,
                        dims = dims,
                        reduction = reduction,
                        k.param = k.param,
                        return.neighbor = T)
  
  neighbors <- data@neighbors$RNA.nn
  
  neighbor_locs <- neighbors@nn.idx
  rownames(neighbor_locs) <- colnames(data)
  
  samp_matrix <- data[["RNA"]]@counts
  samp_matrix <- samp_matrix[rownames(samp_matrix) %in% features,]
  
  set.seed(seed)
  metacell_num <- floor(dim(samp_matrix)[2]/k.param)
  
  metacells <- sample(1:dim(samp_matrix)[2], size = metacell_num, replace = FALSE)
  
  meta_mat <- matrix(0, nrow = dim(samp_matrix)[1], ncol = length(metacells))
  colnames(meta_mat) <- colnames(samp_matrix)[metacells]
  rownames(meta_mat) <- rownames(samp_matrix)

  pb <- txtProgressBar(min = 0,
                       max = length(metacells),
                       initial = 0,
                       width = 100,
                       style = 3)

  for(i in 1:length(metacells)){
    cell_name <- colnames(samp_matrix)[metacells[i]]
    neighbor_loc <- neighbor_locs[metacells[i],]
    neighbor_profiles <- samp_matrix[,neighbor_loc]
    neighbor_profiles <- apply(neighbor_profiles, MARGIN = 2, FUN = function(x){x/sum(x)})
    meta_mat[,which(colnames(meta_mat) == cell_name)] <- apply(neighbor_profiles, MARGIN = 1, mean)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  

  saveRDS(meta_mat, file = file.path(res_path, "metacell_matrix.rds"))
  saveRDS(data, file = file.path(res_path, "metacell_prep_data.rds"))
  
  return(meta_mat)
}
