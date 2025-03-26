plot_milo_nhood_size <- function(milo, bins = 50, res_path = NULL, save = T){
  p <- plotNhoodSizeHist(milo) +
    geom_histogram(bins = bins, fill = proj_cols("yellow"), color = proj_cols_grey("dark grey")) +
    labs(x = "Nhood Size", y = "Frequency")+
    theme_proj()+
    theme(legend.position = "none",
          plot.margin = margin(t = 0.1, r = 0.2, b = 0.1, l = 0.1, 'in'))
  p <- p + geom_vline(xintercept = mean(p$data$nh_size), color = proj_cols("red"))
  p <- p + geom_vline(xintercept = median(p$data$nh_size), color = proj_cols("light blue"))
  
  if(save){
    ggsave(file.path(res_path, "milo_nhood_size_dist.png"), height = 3, width = 4, units = "in", dpi = 600)
  }
}


get_nhood_comp <- function(milo = NULL, res_path = NULL){
  nhood_comp <- list()
  for(j in 1:length(colnames(milo@nhoods))){
    keep_cells <- colnames(milo)[milo@nhoods[,j]==1]
    nhood_comp <- c(nhood_comp, list(keep_cells))
  }
  saveRDS(nhood_comp, file.path(res_path, "milo_nhood_membership.rds"))
}

get_nhood_schpf_scores <- function(milo = NULL,
                                   da_results = NULL,
                                   cell_scores = NULL,
                                   sum_cols = NULL,
                                   res_path = NULL){
  
  all_agg <- data.frame()
  
  for(j in 1:length(colnames(milo@nhoods))){
    keep_cells <- colnames(milo)[milo@nhoods[,j]==1]
    
    aggregate <- cell_scores %>%
      filter(cell_names %in% keep_cells) %>%
      summarise_at(.vars = all_of(unname(sum_cols)), list(median = median))
    names(aggregate) <- gsub("_median", "", names(aggregate))
    all_agg <- rbind(all_agg, cbind(Nhood = j, aggregate))
  }
  all_agg <- left_join(all_agg, da_results)
  saveRDS(all_agg, file.path(res_path, "milo_agg_factors.rds"))
  return(all_agg)
}

plot_milo_qc <- function(res = NULL, bins = 50, res_path = NULL){
  
  p1 <- ggplot(res, aes(PValue)) +
    geom_histogram(bins = bins, fill = proj_cols("light blue"), color = proj_cols_grey("dark grey")) +
    labs(x = "P-value", y = "Frequency")+
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75","1"))+
    theme_proj()
  
  p2 = res %>%
    mutate(col = ifelse(SpatialFDR < 0.05 & logFC < 0, "1", "2"),
           col = ifelse(SpatialFDR < 0.05 & logFC > 0, "3", col)) %>%
    ggplot(aes(x = logFC, y = -log10(SpatialFDR), color = col)) +
    geom_point() +
    geom_hline(yintercept = 1, lty = 2, color = proj_cols_grey("med grey")) +
    labs(x = "logFC", y = expression('log'[10]*'(Spatial FDR)'))+
    scale_color_manual(values = unname(c(proj_cols("blue"), proj_cols_grey("med grey"), proj_cols("pink"))))+
    theme_proj()+
    theme(legend.position = "none")
  
  if("trait_fraction" %in% names(res)){
    p3 <- ggplot(res, aes(trait_fraction)) +
      geom_histogram(bins = bin, fill = proj_cols("teal"), color = proj_cols_grey("dark grey")) +
      geom_vline(xintercept = 0.7, lty = 2, color = proj_cols("red")) +
      labs(x = "Fraction Assigned", y = "Frequency")+
      theme_proj()
  }
  
  plots <- eval(parse(text = paste0("list(", paste0(grep("p1|p2|p3", ls(), v=T), collapse= ","), ")")))
  p <- wrap_plots(plots)
  
  ggsave(file.path(res_path, "milo_qc_plots.png"),
         height = 3, width = 10, units = "in", dpi = 600)
}

plot_milo_res <- function(milo = NULL, da_results = NULL, res_path = NULL, label_size = 5, spatial.fdr = 0.05){
  
  n_groups <- length(unique(da_results$NhoodGroup))
  
  require(ggraph)
  
  set.seed(1234)
  group_cols <- colorspace::lighten(sample(pal()(n_groups)), 0.2, space = "HLS")
  
  p1 <- plotNhoodGroups(milo, da_results, layout = "SCHPF.UMAP")
  
  # ensure the order of the legend is in numerical order
  order <- sort(as.numeric(unique(p1$data$colour_by)))
  p1$data$colour_by <- factor(p1$data$colour_by, levels = order)
  
  p1 <- p1 +
    scale_fill_manual(name = "Nhood Group", values = group_cols)+
    scale_size(name = "Nhood Size", range = c(0.5, 3), guide = "none")+
    scale_edge_width(name = "Overlap Size", range = c(0.2, 3), guide ="none") +
    guides(fill = guide_legend(override.aes = aes(size = 8, alpha = 1, stroke = 0),
                               title = "Meta-neighborhood identity", 
                               position = "bottom", 
                               direction = "horizontal", 
                               nrow = 1, 
                               keywidth = 0.15,
                               keyheight = 0.2,
                               default.unit = "in", 
                               theme = theme(legend.title.position = "top",
                                             legend.title = element_text(size = 10, 
                                                                         family = "Helvetica", 
                                                                         margin = margin(t = 0, r = 0, b = 0.05, l = 0, 'in')),
                                             legend.margin = margin(0, 0, 0, 0),
                                             legend.text = element_text(size = 10, 
                                                                        family = "Helvetica", 
                                                                        hjust = 0.5,
                                                                        margin = margin(l = -0.2, 0, 0, 0, 'in'), 
                                                                        color = "white"),
                                             legend.key.spacing.x = unit(0.05, "in"))))
  
  annot <- p1$data %>% 
    group_by(colour_by) %>% 
    summarise(x = median(x), y = median(y))
  
  p1 <- p1 + 
    geom_point(data = annot, aes(x = x, y = y), 
               color = "black", 
               size = 5.5, 
               alpha = 0.4) +
    geom_text(data = annot, aes(x = x, y = y, label = colour_by), 
              color = "white",
              size = label_size, 
              fontface = "bold",
              family = "Helvetica")
  
  pdf(file.path(res_path, "milo_meta_neighborhoods.pdf"), height = 4.5, width = 4)
  print(p1)
  dev.off()
  
  p2 <- plotNhoodGraphDA(milo, da_results, alpha = spatial.fdr, layout ="SCHPF.UMAP") +
    scale_fill_gradient2(name = "logFC",
                         low = unname(proj_cols("blue")),
                         high = unname(proj_cols("red")),
                         mid = "white")+
    scale_size(name = "Nhood\nSize", range = c(0.5, 3))+
    scale_edge_width(name = "Overlap\nSize", range = c(0.2, 3)) +
    guides(size = guide_legend(order = 2, 
                               byrow = T, 
                               theme = theme(legend.title = element_text(margin = margin(b = 0, unit = "in")), 
                                             legend.key.spacing.y = unit(-0.05, "in"),
                                             legend.text = element_text(margin = margin(l = 0.05, unit = "in")))), 
           fill = guide_colorbar(order = 1, 
                                 barheight = unit(1, "in"),  
                                 theme = theme(legend.title = element_text(margin = margin(b = 0.1, unit = "in")),
                                               legend.text = element_text(margin = margin(l = 0.05, unit = "in")))),
           edge_width = guide_legend(order = 3, 
                                     byrow = T, 
                                     theme = theme(legend.title = element_text(margin = margin(b = 0, unit = "in")), 
                                                   legend.key.spacing.y = unit(-0.05, "in"),
                                                   legend.text = element_text(margin = margin(l = 0.05, unit = "in")))))+
    theme(legend.margin = margin(t=0, r=0, b=0, l=-0.1, "in"),
          legend.text = element_text(size = 10, family = "Helvetica"),
          legend.title = element_text(size = 10, family = "Helvetica"))
  
  pdf(file.path(res_path, "milo_nhood_logFC.pdf"), height = 4, width = 4.5)
  print(p2)
  dev.off()
  
  res <- data.frame()
  for(i in as.numeric(as.character(unique(da_results$NhoodGroup)))) {
    group1 <- da_results$logFC[da_results$NhoodGroup==i]
    group2 <- da_results$logFC[da_results$NhoodGroup!=i]
    
    if(length(group1) > 5){
      out <- t.test(group1, group2)
      res <- rbind(res, cbind(nhood = i, 
                              size = length(group1),
                              t = out$statistic, 
                              p = out$p.value, 
                              diff = unname(out$estimate[1]-out$estimate[2]), 
                              lower = out$conf.int[1], 
                              upper = out$conf.int[2]))
    } else {
      res <- rbind(res, cbind(nhood = i, 
                              size = length(group1),
                              t = NA, 
                              p = NA, 
                              diff = NA, 
                              lower = NA, 
                              upper = NA))
    }
  }
  
  res$p.adj <- p.adjust(res$p, method = "BH")
  res$nhood <- factor(res$nhood, levels = order)
  saveRDS(res, file.path(res_path, "comparison_test.rds"))
  
  p3 <- milo_beeswarm(da_results, res, group_by = "NhoodGroup", alpha = spatial.fdr)
  
  pdf(file.path(res_path, "milo_beeswarm.pdf"), height = 4, width = 4.5)
  print(p3)
  dev.off()
  
}


milo_beeswarm <- function(da_results, res, group_by = "NhoodGroup", alpha = 0.1){
  
  bee_data <- da_results %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood)), 
           NhoodGroup = factor(NhoodGroup, levels = as.numeric(as.character(unique(da_results$NhoodGroup)))))
  
  cols <- c(proj_cols("blue", "red"), proj_cols_grey("med grey"))
  names(cols) <- c("1", "2", "3")
  res$cols<-= ifelse(res$t > 0 & res$p.adj < 0.05, "2", ifelse(res$t < 0 & res$p.adj < 0.05, "1", "3"))
  
  set.seed(1234)
  ggplot() +
    ggbeeswarm::geom_quasirandom(data = bee_data, aes(!!sym(group_by), logFC, color=logFC_color), alpha=1) +
    scale_color_gradient2(name = "logFC",
                          low = unname(proj_cols("blue")),
                          high = unname(proj_cols("red")),
                          mid = "white",
                          na.value = proj_cols_grey("light grey"))+
    geom_hline(yintercept = 0, linewidth = 0.5)+
    geom_errorbar(data = res, aes(ymin = lower, ymax = upper, x = nhood), width = 0.3)+
    geom_point(data = res, aes(y = diff, x = nhood, fill = cols), size = 2, shape = 21, stroke = 0.5) +
    scale_fill_manual(values = cols, guide = "none")+
    scale_y_continuous(labels = function(x) formatC(x, format = "g"))+
    guides(color="none") +
    labs(y = "logFC", x = "Meta-neighborhood")+
    coord_flip()+
    theme_proj()+
    theme(panel.border = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(family = "Helvetica"),
          axis.text.y = element_text(family = "Helvetica"),
          axis.ticks.y = element_blank())
  
}

compared_nhoods <- function(run1_path = NULL, 
                            run2_path = NULL, 
                            name1 = NULL, 
                            name2 = NULL, 
                            file_name = "nhood_comparison.png", 
                            res_path = NULL,
                            fontsize = 12){
  
  require(grid)
  require(ComplexHeatmap)
  
  nhoods1 <- read_rds(file.path(run1_path, "milo_nhood_membership.rds"))
  meta_nhoods1 <- read_rds(file.path(run1_path, "milo_grouped_Nhoods.rds"))
  nhoods2 <- read_rds(file.path(run2_path, "milo_nhood_membership.rds"))
  meta_nhoods2 <- read_rds(file.path(run2_path, "milo_grouped_Nhoods.rds"))
  
  all_shared <- length(intersect(unlist(nhoods1), unlist(nhoods2)))
  all_unique <- length(unique(c(unlist(nhoods1), unlist(nhoods2))))
  
  sim <- data.frame()
  
  for(i in unique(meta_nhoods1$NhoodGroup)){
    for(j in unique(meta_nhoods2$NhoodGroup)){
      
      keep_nhoods1 <- rownames(meta_nhoods1 %>% filter(NhoodGroup == i))
      cells1 <- unique(unlist(nhoods1[as.numeric(keep_nhoods1)]))
      keep_nhoods2 <- rownames(meta_nhoods2 %>% filter(NhoodGroup == j))
      cells2 <- unique(unlist(nhoods2[as.numeric(keep_nhoods2)]))
      
      int <- length(intersect(cells1, cells2))/length(unique(c(cells1, cells2)))
      sim <- rbind(sim, cbind(i, j, int))
      
      rm(list = ls()[ls() %in% c("keep_nhoods1", "keep_nhoods2", "cell2", "cell1", "int")])
    }
  }
  sim <- sim %>%
    mutate_at(all_of(c("j", "j", "int")), ~as.numeric(as.character(.))) %>%
    mutate(int = ifelse(is.na(int), 0, int)) %>%
    mutate(i = factor(i, levels = sort(as.numeric(unique(i)))),
           j = factor(j, levels = sort(as.numeric(unique(j))))) %>%
    arrange(i, j)
  
  mat <-sim %>%
    spread(i, int) %>%
    column_to_rownames("j") %>%
    mutate_all(~as.numeric(.)) %>%
    as.matrix()
  
  
  pdf(file.path(res_path, file_name), height = 4, width = 3.5)
  set.seed(1234)
  p<-ComplexHeatmap::Heatmap(mat, name = "mat", na_col = "#F0F0F0",
                             col = circlize::colorRamp2(breaks = c(0, max(mat, na.rm=T)),
                                                        colors = c("white", proj_cols("light blue"))),
                             cell_fun = function(j, i, x, y, w, h, fill) {
                               if(!is.na(mat[i, j])) {
                                 grid.text(if(mat[i, j]>=0.05){sprintf("%.2f", mat[i, j])}, x, y, gp = gpar(fontsize = 6))
                               }
                             },
                             rect_gp = grid::gpar(col = "#F0F0F0", lwd = 1.5),
                             cluster_rows = F,
                             column_title = name1,
                             row_title = name2,
                             column_title_side = "bottom",
                             row_names_side = "left",
                             row_dend_side = "right",
                             row_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = fontsize),
                             column_names_rot = 0,
                             column_names_gp = gpar(fontface = "plain", fontfamily = "Helvetica", fontsize = fontsize),
                             cluster_columns = F,
                             column_names_centered = T,
                             heatmap_legend_param = list(
                               title = "% Overlap",
                               legend_width = unit(2, "in"),
                               direction = "horizontal",
                               title_position = "topcenter",
                               border = "black"
                             ))
  
  draw(p, heatmap_legend_side = "top", background = "transparent", padding = unit(c(0.1, 0.1, 0.1, 0.1), "in"))
  decorate_heatmap_body("mat", { grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 1.5)) })
  dev.off()
}

plot_nhood_dotplot <- function(milo_agg_factors, res_path, factors, height = 4.5, width = 8.5){
  
  groups <- unique(milo_agg_factors$NhoodGroup)
  order <- sort(as.numeric(groups))
  
  res <- data.frame()
  for(i in unname(factors)){
    for(j in groups){
      sel_group <- milo_agg_factors[,i][milo_agg_factors$NhoodGroup==j]
      other_group <- milo_agg_factors[,i][milo_agg_factors$NhoodGroup!=j]
      
      test <- wilcox.test(sel_group, other_group)
      
      res <- rbind(res, cbind(factor = i, 
                              NhoodGroup=j, 
                              mean_sel = mean(sel_group),
                              sd_sel = sd(sel_group),
                              mean_other = mean(other_group),
                              sd_other = sd(other_group),
                              log2fc = mean(log2(sel_group)) - mean(log2(other_group)), 
                              pvalue = test$p.value))
    }
  }
  
  res <- res %>%
    mutate_at(all_of(c("pvalue", "log2fc", "mean_sel", "sd_sel", "mean_other", "sd_other")), ~as.numeric(as.character(.))) %>%
    mutate(p.adj = p.adjust(pvalue, method = "BH"),
           factor = factor(factor, levels = factors, labels = names(factors)),
           NhoodGroup = factor(NhoodGroup, levels = rev(order)))
  
  saveRDS(res, file.path(res_path, "factor_comparison.rds"))
  
  p <- res %>%
    mutate(p.adj.cor = ifelse(p.adj == 0, min(p.adj[p.adj!=0]), p.adj),
           log10.p.adj = ifelse(log2fc > 0, -log10(p.adj.cor), log10(p.adj.cor)),
           log10.p.adj = ifelse(p.adj < 0.001, log10.p.adj, NA)) %>%
    ggplot(aes(x = factor, y = NhoodGroup, size = abs(log2fc), fill = log10.p.adj)) +
    geom_point(aes(color = ifelse(is.na(log10.p.adj), "1", "0")), shape = 21) +
    scale_fill_gradient2(low = proj_cols("blue"), high = proj_cols("red"), midpoint = 0, na.value = proj_cols_grey("med grey")) +
    scale_color_manual(values = c(`1`=unname(proj_cols_grey("med grey")), `0`="black"), breaks = "1", labels = "NS")+
    scale_size(range = c(0.5, 5.5))+
    labs(y = "Meta-neighborhood")+
    theme_proj()+
    theme(plot.margin = margin(l = 1, t = 0.1, r = 0.1, b = 0.1,  unit = "in"), 
          axis.text.x = element_text(angle = 45, hjust = 1, family = "Helvetica", size = 12),
          axis.text.y = element_text(family = "Helvetica", size = 12),
          axis.title.y = element_text(family = "Helvetica", size = 12),
          panel.grid.major = element_blank(),
          panel.border = element_rect(linewidth = 1, color = "black"),
          axis.line.y = element_line(linewidth = .25, color = "black"),
          axis.line.x = element_line(linewidth = .25, color = "black"),
          axis.title.x = element_blank()) +
    guides(color = guide_legend(order = 2, 
                                override.aes = aes(shape = 15, size = 4, stroke = 1.5), 
                                theme = theme(legend.title = element_blank(),
                                              legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
                                              legend.margin = margin(l = -0.1, t = -0.2, unit = "in"))),
           size = guide_legend(title = "Avg abs<br>log<sub>2</sub>FC", 
                               order = 3, 
                               byrow = T, 
                               theme = theme(legend.title = element_markdown(size = 10, margin = margin(b = 0.05, unit = "in")), 
                                             legend.key.spacing.y = unit(0.01, "in"),
                                             legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
                                             legend.margin = margin(l=-0.1, unit = "in"))), 
           fill = guide_colorbar(title = "-log<sub>10</sub>FDR", 
                                 order = 1, 
                                 barheight = unit(0.75, "in"), 
                                 barwidth = unit(0.15, "in"),
                                 theme = theme(legend.title = element_markdown(size = 10, margin = margin(b = 0.05, unit = "in")),
                                               legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
                                               legend.margin = margin(l=-0.1, unit = "in"),
                                               legend.ticks = element_line(linewidth = 0.5, color = "white"),
                                               legend.ticks.length = unit(0.025, "in")))) 
  
  cairo_pdf(file.path(res_path, "milo_factor_comparison.pdf"), height = height, width = width)
  print(p)
  dev.off()
  
  return(p)
  
}
