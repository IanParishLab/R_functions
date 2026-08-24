plot_UCell <- function(seu.obj, gn, reduction = "pca", reduction.name, cluster.ident, 
                       label = TRUE, alpha = 0.4, cluster_col = NULL,
                       save.loc = "plots", plot.name, plot.width = length(gn)*3, plot.height = 6,
                       vln = FALSE, ncol = length(gn), raster = FALSE, print = FALSE){
  suppressPackageStartupMessages({
    require(Seurat)
    require(UCell)
    require(patchwork)
    require(RColorBrewer)
  })
  
  # make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # get colour paletter if not provided
  if (is.null(cluster_col) && exists("cluster_col", envir = .GlobalEnv)) {
    cluster_col <- get("cluster_col", envir = .GlobalEnv)
  }
  
  # Run UCell
  Idents(seu.obj) <- cluster.ident
  seu.obj1 <- AddModuleScore_UCell(seu.obj, features = gn)
  signature.names <- paste0(names(gn), "_UCell")
  
  # `SmoothKNN` function performs smoothing of single-cell scores by weighted average of the k-nearest neighbors. It can be useful to 'impute' scores by neighboring cells and partially correct data sparsity.
  seu.obj1 <- SmoothKNN(seu.obj1, signature.names = signature.names, reduction = reduction)
  
  if(isTRUE(label)){
    fp <- FeaturePlot(seu.obj1, reduction = reduction.name, order = TRUE, repel = TRUE, label = TRUE, alpha = alpha,
                      features = c(signature.names, paste0(signature.names,"_kNN")), ncol = ncol, raster = raster)
  } else {
    fp <- FeaturePlot(seu.obj1, reduction = reduction.name, order = TRUE, repel = TRUE, alpha = alpha,
                      features = c(signature.names, paste0(signature.names,"_kNN")), ncol = ncol, raster = raster)
  }
  
  fp <- fp &
    NoAxes() &
    theme(aspect.ratio = 1, text = element_text(size = 7), 
          axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  plot_features <- c(signature.names, paste0(signature.names,"_kNN"))
  if(isTRUE(vln)){
    vln_p <- VlnPlot(seu.obj1, group.by = cluster.ident,
                     features = plot_features, 
                     cols = cluster_col, ncol = 2, pt.size = 0, raster = raster) &
      geom_boxplot(width = 0.1, color = "black", alpha = alpha) &
      theme(text = element_text(size = 7), 
            axis.text = element_text(size = 7),
            axis.title.x = element_text(size = 7), 
            axis.title.y = element_text(size = 7))
  }
  
  if(isTRUE(print)){
    print(fp)
  }
  
  # save plot
  png(file.path(save.loc, "plots", paste0(plot.name, ".plotUCell.png")), width = plot.width, height = plot.height, units = "in", res = 300)
  print(fp)
  dev.off()
  
  if(isTRUE(vln)){
    print(vln_p)
  }
  
  return(seu.obj1)
}
