FOV_DimPlot <- function(data = NULL, 
                        cell_name = NULL, 
                        boundary = 100, 
                        fov = NULL, 
                        assay = "Vizgen", 
                        molecules = NULL, 
                        zoom = NULL,
                        group.by = NULL, 
                        cols = NULL, 
                        bar_col = "white",
                        dectected_transcript = NULL, 
                        pt.size = 0.5, 
                        fill.alpha = 1, 
                        mol_cols = NULL,
                        pt.scale = 1, 
                        bar.text.size = 3, ...){
  
  if(!is.null(cell_name)){
    ind <- which(data@images[[fov]]$centroids@cells %in% cell_name)
    x <- data@images[[fov]]$centroids@coords[ind,"x"]
    y <- data@images[[fov]]$centroids@coords[ind,"y"]
  }
  
  cell_coords <- data.frame(cell_names = data@images[[fov]]@boundaries$centroids@cells,
                           x = data@images[[fov]]@boundaries$centroids@coords[,1],
                           y = data@images[[fov]]@boundaries$centroids@coords[,2])
  
  if(!is.null(zoom)){
    min_x <- ifelse(zoom[1] <= max(cell_coords$x) & zoom[1] >= min(cell_coords$x), zoom[1], min(cell_coords$x))
    max_x <- ifelse(zoom[2] >= min(cell_coords$x) & zoom[2] <= max(cell_coords$x), zoom[2], max(cell_coords$x))
    min_y <- ifelse(zoom[3] <= max(cell_coords$y) & zoom[3] >= min(cell_coords$y), zoom[3], min(cell_coords$y))
    max_y <- ifelse(zoom[4] >= min(cell_coords$y) & zoom[4] <= max(cell_coords$y), zoom[4], max(cell_coords$y))
    
  } else {
    min_y <- ifelse(y-boundary <= max(cell_coords$y) & y-boundary >= min(cell_coords$y), y - boundary, max(cell_coords$y)-boundary)
    max_y <- ifelse(y+boundary >= min(cell_coords$y) & y+boundary <= max(cell_coords$y), y + boundary, max(cell_coords$y))
    min_x <- ifelse(x-boundary <= max(cell_coords$x) & x-boundary >=  min(cell_coords$x), x - boundary, max(cell_coords$x)-boundary)
    max_x <- ifelse(x+boundary >= min(cell_coords$x) & x+boundary <= max(cell_coords$x), x + boundary,max(cell_coords$x))
  }

  # FIND CENTROIDS IN THE ZOOM AREA
  #-----------------------------------------------------------------------------
  cells <- data@images[[fov]]$centroids@cells
  centroids <- data.frame(cells, data@images[[fov]]$centroids@coords) %>%
    filter(x >= min_x & x <= max_x & y >= min_y & y <= max_y)
  
  # GET POLYGONS & COORDS
  #-----------------------------------------------------------------------------
  polygons <- lapply(centroids$cells, function(x) sf::st_polygon(list(data@images[[fov]]@boundaries$segmentation@polygons[[x]]@Polygons[[1]]@coords)))
  
  # extract individual polygon box
  bbox_list <- lapply(polygons, function(x) st_bbox(st_geometry(x)))
  poly_coords <- data.frame(matrix(unlist(bbox_list), nrow = length(bbox_list), byrow=TRUE))
  names(poly_coords) <- names(bbox_list[[1]])
  
  x_min <- min(min(poly_coords$xmin), zoom[1])
  x_max <- max(max(poly_coords$xmax), zoom[2])
  y_min <- min(min(poly_coords$ymin), zoom[3])
  y_max <- max(max(poly_coords$ymax), zoom[4])
  
  metadata <- data.frame(cell_names = centroids$cells) %>% 
    left_join(data@meta.data %>% dplyr::select(cell_names, all_of(group.by))%>% rename_at(all_of(group.by), ~"col"))
  
  metadata$geometry <- NULL
  
  for(j in 1:nrow(metadata)){
    metadata$geometry[j] <- st_geometry(polygons[[j]])
  }
  
  polygons <- st_as_sf(metadata, sf_column_name="geometry")
  
  # FINAL PLOT
  #-----------------------------------------------------------------------------
  p <- ggplot() +
    theme_proj() +
    geom_blank(data = data.frame(x = c(x_min, x_max), y = c(y_min, y_max)), aes(x, y)) +
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          plot.title = element_text(family = "Helvetica", size = 14, hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, -0.1, "in"),
          legend.spacing.x = grid::unit(-0.05, "in"),
          legend.key.width = grid::unit(0.2, unit = "in"),
          legend.key.height = grid::unit(0.2, unit = "in"),
          legend.text = element_text(size =10, family = "Helvetica", margin = margin(l = 0.1, unit = "in")),
          panel.border = element_rect(fill = NA, color = "black", size = 1),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.1),
          axis.line.y = element_line(color = "black", linewidth = 0.1),
          axis.ticks.x =  element_line(color = "black", linewidth = 0.5),
          axis.ticks.y =  element_line(color = "black", linewidth = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size =10, family = "Helvetica"),
          axis.text.y = element_text(size =10, family = "Helvetica")) +
    geom_sf(data = polygons, aes(fill = col), color = "white", show.legend = T, linewidth = 0, alpha = fill.alpha) +
    scale_fill_manual(values = cols)
  
  if(!is.null(molecules)){
    
    mol_data <- detected_transcripts %>%
      mutate(cell_id = as.character(cell_id)) %>%
      filter(gene %in% molecules & cell_id %in% centroids$cells) %>%
      filter(global_x > min_x & global_x < max_x & global_y > min_y & global_y < max_y) %>%
      mutate(gene = factor(gene, levels = molecules))
    
    if(is.null(mol_cols)){
      mol_cols <- DiscretePalette_scCustomize(num_colors = length(molecules), palette = "polychrome")
      names(mol_cols) <- molecules
    }
    
    p <- p + 
      geom_point_rast(data = mol_data, aes(x = global_x, y = global_y, color = gene), 
                      size = pt.size, 
                      scale = pt.scale, alpha = 1)+
      scale_color_manual(values = mol_cols)+
      guides(color = guide_legend(override.aes=list(fill=NA, size = 3)))
    
  }
  
  p <- p +
    geom_segment(aes(x = min_x + 25, xend = min_x + 25 + (boundary/4), 
                     y = min_y + 25, yend =  min_y + 25), 
                 color = "black", linewidth = 1) +
    annotate("text", x = ((min_x + 25 + (boundary/4))+(min_x + 25))/2, 
             y = min_y + 35, label=paste0((boundary/4), "μm"), 
             hjust = 0.5, vjust = 0, size = bar.text.size)
  return(p)
  
}

FOV_FeaturePlot <- function(data = NULL, 
                           cell_name = NULL, 
                           boundary = 100, 
                           fov = NULL, 
                           cell_subset = NULL,
                           assay = "Vizgen", feature = NULL, zoom = NULL,
                           bar_col = "white", bar.text.size = 3, ...){
  
  if(!is.null(cell_name)){
    ind <- which(data@images[[fov]]$centroids@cells %in% cell_name)
    x <- data@images[[fov]]$centroids@coords[ind,"x"]
    y <- data@images[[fov]]$centroids@coords[ind,"y"]
  }
  
  cell_coords <- data.frame(cell_names = data@images[[fov]]@boundaries$centroids@cells,
                           x = data@images[[fov]]@boundaries$centroids@coords[,1],
                           y = data@images[[fov]]@boundaries$centroids@coords[,2])
  
  if(!is.null(zoom)){
    min_x <- ifelse(zoom[1] <= max(cell_coords$x) & zoom[1] >= min(cell_coords$x), zoom[1], min(cell_coords$x))
    max_x <- ifelse(zoom[2] >= min(cell_coords$x) & zoom[2] <= max(cell_coords$x), zoom[2], max(cell_coords$x))
    min_y <- ifelse(zoom[3] <= max(cell_coords$y) & zoom[3] >= min(cell_coords$y), zoom[3], min(cell_coords$y))
    max_y <- ifelse(zoom[4] >= min(cell_coords$y) & zoom[4] <= max(cell_coords$y), zoom[4], max(cell_coords$y))
    
  } else {
    min_y <- ifelse(y-boundary <= max(cell_coords$y) & y-boundary >= min(cell_coords$y), y - boundary, max(cell_coords$y)-boundary)
    max_y <- ifelse(y+boundary >= min(cell_coords$y) & y+boundary <= max(cell_coords$y), y + boundary, max(cell_coords$y))
    min_x <- ifelse(x-boundary <= max(cell_coords$x) & x-boundary >=  min(cell_coords$x), x - boundary, max(cell_coords$x)-boundary)
    max_x <- ifelse(x+boundary >= min(cell_coords$x) & x+boundary <= max(cell_coords$x), x + boundary,max(cell_coords$x))
  }
  
  # FIND CENTROIDS IN THE ZOOM AREA
  #-----------------------------------------------------------------------------
  cells <- data@images[[fov]]$centroids@cells
  centroids <- data.frame(cells, data@images[[fov]]$centroids@coords) %>%
    filter(x >= min_x & x <= max_x & y >= min_y & y <= max_y)
  
  # GET POLYGONS & COORDS
  #-----------------------------------------------------------------------------
  polygons <- lapply(centroids$cells, function(x) sf::st_polygon(list(data@images[[fov]]@boundaries$segmentation@polygons[[x]]@Polygons[[1]]@coords)))
  
  # extract individual polygon box
  bbox_list <- lapply(polygons, function(x) st_bbox(st_geometry(x)))
  poly_coords <- data.frame(matrix(unlist(bbox_list), nrow = length(bbox_list), byrow=TRUE))
  names(poly_coords) <- names(bbox_list[[1]])
  
  x_min <- min(min(poly_coords$xmin), zoom[1])
  x_max <- max(max(poly_coords$xmax), zoom[2])
  y_min <- min(min(poly_coords$ymin), zoom[3])
  y_max <- max(max(poly_coords$ymax), zoom[4])
  
  if(feature %in% colnames(data@meta.data)){
    feature_data <- data.frame(cell_names = rownames(data@meta.data), feature = data@meta.data[,feature])
  } else {
    feature_data <- data.frame(cell_names = colnames(data[[assay]]@data), feature = data[[assay]]@data[feature,])
  }
  
  if(!is.null(cell_subset)){
    feature_data <- feature_data %>% mutate(feature = ifelse(cell_names %in% cell_subset, feature, NA))
  } 
  
  metadata <- data.frame(cell_names = centroids$cells) %>% 
    left_join(feature_data)
  
  metadata$geometry <- NULL
  
  for(j in 1:nrow(metadata)){
    metadata$geometry[j] <- st_geometry(polygons[[j]])
  }
  
  polygons <- st_as_sf(metadata, sf_column_name="geometry")
  
  # FINAL PLOT
  #-----------------------------------------------------------------------------
  p <- ggplot() +
    theme_proj() +
    geom_blank(data = data.frame(x = c(x_min, x_max), y = c(y_min, y_max)), aes(x, y)) +
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.margin = margin(0.05, 0.05, 0.05, 0.05, unit = "in"),
          plot.title = element_text(family = "Helvetica", size = 14, hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, -0.1, "in"),
          legend.spacing.x = unit(-0.05, "in"),
          legend.key.width = unit(0.2, unit = "in"),
          legend.key.height = unit(0.2, unit = "in"),
          legend.text = element_text(size =10, family = "Helvetica", margin = margin(l = 0.1, unit = "in")),
          panel.border = element_rect(fill = NA, color = "black", size = 1),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(color = "black", linewidth = 0.1),
          axis.line.y = element_line(color = "black", linewidth = 0.1),
          axis.ticks.x =  element_line(color = "black", linewidth = 0.5),
          axis.ticks.y =  element_line(color = "black", linewidth = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size =10, family = "Helvetica"),
          axis.text.y = element_text(size =10, family = "Helvetica"))+
    scale_fill_gradient(low = "lightgrey", high = "blue", na.value = "#F0F0F0")
  
  if(is.null(cell_subset)){
    p <- p + geom_sf(data = polygons, aes(fill = feature), color = "white", show.legend = T, linewidth = 0) 
  } else {
    p <- p + 
      geom_sf(data = subset(polygons, subset = cell_names %in% cell_subset), 
              aes(fill = feature), color = "white", show.legend = T, linewidth = 0) +
      geom_sf(data = subset(polygons, subset = !cell_names %in% cell_subset), 
              fill = "white", color = "lightgrey", show.legend = T, linewidth = 0.2) 
  }
  
  p <- p +
    geom_segment(aes(x = min_x + 25, xend = min_x + 25 + (boundary/4), 
                     y = min_y + 25, yend =  min_y + 25), 
                 color = "black", linewidth = 1) +
    annotate("text", x = ((min_x + 25 + (boundary/4))+(min_x + 25))/2, 
             y = min_y + 35, label=paste0((boundary/4), "μm"), 
             hjust = 0.5, vjust = 0, size = bar.text.size)

  return(p)
  
}

Clustered_DotPlot_Single_Group_v2 <- function(
    seurat_object,
    features,
    identity_lab = "",
    exp_title = "Avg. Exp",
    pct_title = "% Exp",
    colors_use_exp = viridis_plasma_dark_high,
    exp_color_min = -2,
    exp_color_middle = NULL,
    exp_color_max = 2,
    max_col_size = 2,
    print_exp_quantiles = FALSE,
    colors_use_idents = NULL,
    x_lab_rotate = TRUE,
    plot_padding = NULL,
    flip = FALSE,
    k = 1,
    feature_km_repeats = 1000,
    ident_km_repeats = 1000,
    row_label_size = 8,
    row_label_fontface = "plain",
    grid_color = NULL,
    cluster_feature = TRUE,
    cluster_ident = TRUE,
    column_label_size = 8,
    legend_label_size = 10,
    legend_title_size = 10,
    raster = FALSE,
    plot_km_elbow = TRUE,
    elbow_kmax = NULL,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    show_parent_dend_line = TRUE,
    ggplot_default_colors = FALSE,
    color_seed = 123,
    seed = 123
) {
  
  
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)
  
  # set padding
  if (!is.null(x = plot_padding)) {
    if (isTRUE(x = plot_padding)) {
      # Default extra padding
      # 2 bottom: typically mirrors unpadded plot
      # 15 left: usually enough to make rotated labels fit in plot window
      padding <- unit(c(2, 15, 0, 0), "mm")
    } else {
      if (length(x = plot_padding) != 4) {
        cli_abort(message = c("{.code plot_padding} must be numeric vector of length 4 or TRUE",
                              "i" = "Numeric vector will correspond to amount of padding to be added to bottom, left, top, right).",
                              "i" = "Seeting {.field TRUE} will set padding to {.code c(2, 10, 0, 0)}",
                              "i" = "Default is {.val NULL} for no extra padding."))
      }
      padding <- unit(plot_padding, "mm")
    }
  }
  
  # Check acceptable fontface
  if (!row_label_fontface %in% c("plain", "bold", "italic", "oblique", "bold.italic")) {
    cli_abort(message = c("{.code row_label_face} {.val {row_label_face}} not recognized.",
                          "i" = "Must be one of {.val plain}, {.val bold}, {.val italic}, {.val olique}, or {.val bold.italic}."))
  }
  
  # Check unique features
  features_unique <- unique(x = features)
  
  if (length(x = features_unique) != length(x = features)) {
    cli_warn("Feature list contains duplicates, making unique.")
  }
  
  
  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    cli_abort(message = c("Expression color min/max values are not compatible.",
                          "i" = "The value for {.code exp_color_min}: {.field {exp_color_min}} must be less than the value for {.code exp_color_max}: {.field {exp_color_max}}.")
    )
  }
  
  # Get DotPlot data
  seurat_plot <- DotPlot(object = seurat_object, features = features, assay = assay, group.by = group.by, scale = TRUE, idents = idents, col.min = NULL, col.max = NULL)
  
  data <- seurat_plot$data
  
  # Get expression data
  exp_mat <- data %>%
    select(-any_of(c("pct.exp", "avg.exp"))) %>%
    pivot_wider(names_from = any_of("id"), values_from = any_of("avg.exp.scaled")) %>%
    as.data.frame()
  
  row.names(x = exp_mat) <- exp_mat$features.plot
  
  # Check NAs if idents
  if (!is.null(x = idents)) {
    # Find NA features and print warning
    excluded_features <- exp_mat[rowSums(is.na(x = exp_mat)) > 0,] %>%
      rownames()
    cli_warn(message = c("Some scaled data missing.",
                         "*" = "The following features were removed as there is no scaled expression present in subset (`idents`) of object provided:",
                         "i" = "{.field {glue_collapse_scCustom(input_string = excluded_features, and = TRUE)}}.")
    )
    
    # Extract good features
    good_features <- rownames(x = exp_mat)
    
    # Remove rows with NAs
    exp_mat <- exp_mat %>%
      filter(.data[["features.plot"]] %in% good_features)
  }
  
  exp_mat <- exp_mat[,-1] %>%
    as.matrix()
  
  # Get percent expressed data
  percent_mat <- data %>%
    select(-any_of(c("avg.exp", "avg.exp.scaled"))) %>%
    pivot_wider(names_from = any_of("id"), values_from = any_of("pct.exp")) %>%
    as.data.frame()
  
  row.names(x = percent_mat) <- percent_mat$features.plot
  
  # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
  if (!is.null(x = idents)) {
    percent_mat <- percent_mat %>%
      filter(.data[["features.plot"]] %in% good_features)
  }
  
  percent_mat <- percent_mat[,-1] %>%
    as.matrix()
  
  # print quantiles
  if (isTRUE(x = print_exp_quantiles)) {
    cli_inform(message = "Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }
  
  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use_idents) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }
  if (is.null(x = colors_use_idents)) {
    # set default plot colors
    colors_use_idents <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }
  
  # Reduce color length list due to naming requirement
  colors_use_idents <- colors_use_idents[1:group_by_length]
  
  # Modify if class = "colors"
  if (inherits(x = colors_use_idents, what = "colors")) {
    colors_use_idents <- as.vector(x = colors_use_idents)
  }
  
  # Pull Annotation and change colors to ComplexHeatmap compatible format
  Identity <- colnames(x = exp_mat)
  
  identity_colors <- colors_use_idents
  names(x = identity_colors) <- Identity
  identity_colors_list <- list(Identity = identity_colors)
  
  # check grid color
  if (is.null(x = grid_color)) {
    grid_color <- NA
  } else {
    if (length(x = grid_color) > 1) {
      cli_abort(message = "{.code grid_color} can only be a single value.")
    }
    if (isTRUE(x = Is_Color(colors = grid_color))) {
      grid_color <- grid_color
    } else {
      cli_abort(message = "Value provided to {.code grid_color} ({.field {grid_color}}) is not valid value for color in R.")
    }
  }
  
  # Create identity annotation
  if (isTRUE(x = flip)) {
    column_ha <- ComplexHeatmap::rowAnnotation(Identity = Identity,
                                               annotation_label= identity_lab,
                                               annotation_name_gp = gpar(fontfamily = "Helvetica", fontsize = column_label_size),
                                               col =  identity_colors_list,
                                               gp = gpar(col = "black"),
                                               na_col = "grey",
                                               name = identity_lab,
                                               show_legend = FALSE
    )
  } else {
    column_ha <- ComplexHeatmap::HeatmapAnnotation(Identity = Identity,
                                                   annotation_label= identity_lab,
                                                   annotation_name_gp = gpar(fontfamily = "Helvetica", fontsize = column_label_size),
                                                   col =  identity_colors_list,
                                                   na_col = "grey",
                                                   gp = gpar(col = "black"),
                                                   name = identity_lab,
                                                   show_legend = FALSE
    )
  }
  
  # Set middle of color scale if not specified
  if (is.null(x = exp_color_middle)) {
    exp_color_middle <- Middle_Number(min = exp_color_min, max = exp_color_max)
  }
  
  palette_length <- length(x = colors_use_exp)
  palette_middle <- Middle_Number(min = 0, max = palette_length)
  pal_mid_bottom <- (0+palette_middle)/2
  pal_mid_top <- (0+palette_length)/2
  
  # Create palette
  col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])
  # col_fun = colorRamp2(c(exp_color_min, 
  #                        (exp_color_min+exp_color_middle)/2,
  #                        exp_color_middle, 
  #                        (exp_color_middle + exp_color_max)/2,
  #                        exp_color_max), colors_use_exp[c(1, pal_mid_bottom, palette_middle, pal_mid_top, palette_length)])
  
  
  # prep heatmap
  if (isTRUE(x = flip)) {
    if (isTRUE(x = raster)) {
      layer_fun_flip = function(i, j, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(max_col_size, "mm"),
                    gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = "black"))
      }
    } else {
      cell_fun_flip = function(i, j, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(max_col_size, "mm"),
                    gp = gpar(fill = col_fun(exp_mat[i, j]), col = "black"))
      }
    }
  } else {
    if (isTRUE(x = raster)) {
      layer_fun = function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(max_col_size, "mm"),
                    gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = "black"))
      }
    } else {
      cell_fun = function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(max_col_size, "mm"),
                    gp = gpar(fill = col_fun(exp_mat[i, j]), col = "black"))
      }
    }
  }
  
  # Create legend for point size
  lgd_list = list(
    #ComplexHeatmap::Legend(at = Identity, title = identity_lab, legend_gp = gpar(fill = identity_colors_list[[1]]), labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
    ComplexHeatmap::Legend(labels = c(10,25,50,75,100), title = pct_title,
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.50) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                              gp = gpar(fill = "black"))),
                           labels_gp = gpar(fontsize = legend_label_size),
                           title_gp = gpar(fontsize = legend_title_size, fontface = "plain")
    )
  )
  
  # Set x label roration
  if (is.numeric(x = x_lab_rotate)) {
    x_lab_rotate <- x_lab_rotate
  } else if (isTRUE(x = x_lab_rotate)) {
    x_lab_rotate <- 45
  } else {
    x_lab_rotate <- 0
  }
  
  # Create Plot
  set.seed(seed = seed)
  if (isTRUE(x = raster)) {
    if (isTRUE(x = flip)) {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(t(exp_mat),
                                                  heatmap_legend_param=list(title=exp_title, 
                                                                            labels_gp = gpar(fontsize = legend_label_size), 
                                                                            title_gp = gpar(fontsize = legend_title_size, fontface = "plain")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  layer_fun = layer_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  column_km = k,
                                                  row_km_repeats = ident_km_repeats,
                                                  border = "black",
                                                  left_annotation = column_ha,
                                                  column_km_repeats = feature_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_ident,
                                                  cluster_columns = cluster_feature)
    } else {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                  heatmap_legend_param=list(title=exp_title, 
                                                                            labels_gp = gpar(fontsize = legend_label_size), 
                                                                            title_gp = gpar(fontsize = legend_title_size, fontface = "plain")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  layer_fun = layer_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  row_km = k,
                                                  row_km_repeats = feature_km_repeats,
                                                  border = "black",
                                                  top_annotation = column_ha,
                                                  column_km_repeats = ident_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_feature,
                                                  cluster_columns = cluster_ident)
    }
  } else {
    if (isTRUE(x = flip)) {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(t(exp_mat),
                                                  heatmap_legend_param=list(title=exp_title, labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "plain")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  cell_fun = cell_fun_flip,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  column_km = k,
                                                  row_km_repeats = ident_km_repeats,
                                                  border = "black",
                                                  left_annotation = column_ha,
                                                  column_km_repeats = feature_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_ident,
                                                  cluster_columns = cluster_feature)
    } else {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                  heatmap_legend_param=list(title=exp_title, labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "plain")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  cell_fun = cell_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  row_km = k,
                                                  row_km_repeats = feature_km_repeats,
                                                  border = "black",
                                                  top_annotation = column_ha,
                                                  column_km_repeats = ident_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_feature,
                                                  cluster_columns = cluster_ident)
    }
  }
  
  if (!is.null(x = plot_padding)) {
    return(ComplexHeatmap::draw(cluster_dot_plot, 
                                align_heatmap_legend = "heatmap_center",
                                annotation_legend_list = lgd_list, padding = padding, legend_grouping = "original"))
  } else {
    return(ComplexHeatmap::draw(cluster_dot_plot, align_heatmap_legend = "heatmap_center",
                                annotation_legend_list = lgd_list, legend_grouping = "original"))
  }
}

Middle_Number <- function(
    min,
    max
) {
  min_max <- c(min, max)
  middle <- min_max[-length(min_max)] + diff(min_max) / 2
  return(middle)
}


