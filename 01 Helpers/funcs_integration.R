calculate_lisi <- function(data = NULL,
                           split_by = NULL,
                           reduction = "umap"){
  require(lisi)
  
  labels <- data@meta.data %>%
    dplyr::select(dplyr::all_of(split_by))
  
  which_na <- unlist(lapply(labels, function(x) any(is.na(x))))
  split_by <- names(which_na)[which_na == F]
  
  coords <- data[[reduction]]@cell.embeddings
  
  metrics <- data.frame()
  scores <- data.frame()
  
  for(i in 2:ncol(coords)){
    
    res <- lisi::compute_lisi(coords[,c(i-1, i)], labels, split_by)
    
    calc_res <- res %>%
      gather() %>%
      group_by(key) %>%
      summarise(mean = mean(value, na.rm = T),
                sd = sd(value, na.rm = T),
                median = median(value, na.rm = T),
                min = min(value, na.rm = T),
                max = max(value, na.rm = T)) %>%
      mutate(coord1 = dimnames(coords)[[2]][i-1],
             coord2 = dimnames(coords)[[2]][i])
    
    plot_label <- function(x){
      paste0("Mean = ", round(mean(x, na.rm=T), 2), " (", round(sd(x, na.rm=T), 2), ")\n",
             "Med = ", round(median(x, na.rm=T), 2), " [", round(min(x, na.rm=T), 2), "-", round(max(x, na.rm=T), 2), "]")
    }
    
    if(i == 2){
      metrics <- calc_res
      scores <- res %>% rownames_to_column("orig.ident")
    } else {
      metrics <- full_join(metrics, calc_res)
      scores <-  full_join(scores, res %>% rownames_to_column("orig.ident"), by = "orig.ident")
    }
  }
  
  for(i in 2:ncol(scores)){
    
    tmp <- scores[,i]
    names(tmp) <- scores$orig.ident
    
    data <- AddMetaData(object = data,
                       metadata = tmp,
                       col.name = paste0("LISI_", names(scores)[i]))
  }
  
  return(data)
}

