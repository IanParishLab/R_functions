plot_UCell <- function(seu.obj, gn, reduction.name, cluster.ident, label = TRUE, 
                       save.loc = "./plots", plot.name, plot.width = 10, plot.height = 6,
                       vln = FALSE, ncol = 4, raster = FALSE, print = FALSE){
  require(Seurat)
  require(UCell)
  require(patchwork)
  require(RColorBrewer)
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  
  Idents(seu.obj) <- cluster.ident
  
  seu.obj1 <- AddModuleScore_UCell(seu.obj, features = gn)
  
  signature.names <- paste0(names(gn), "_UCell")
  
  # `SmoothKNN` function performs smoothing of single-cell scores by weighted average of the k-nearest neighbors. It can be useful to 'impute' scores by neighboring cells and partially correct data sparsity.
  seu.obj1 <- SmoothKNN(seu.obj1,
                        signature.names = signature.names,
                        reduction="pca",)
  
  if(isTRUE(label)){
    fp <- FeaturePlot(seu.obj1, reduction = reduction.name, order = TRUE, repel = TRUE, label = TRUE,
                      features = c(signature.names, paste0(signature.names,"_kNN")), ncol = ncol, raster = raster)
  } else {
    fp <- FeaturePlot(seu.obj1, reduction = reduction.name, order = TRUE, repel = TRUE,
                      features = c(signature.names, paste0(signature.names,"_kNN")), ncol = ncol, raster = raster)
  }
  
  fp <- fp &
    NoAxes() &
    theme(aspect.ratio = 1, text = element_text(size = 7), 
          axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  if(isTRUE(vln)){
    vln_p <- VlnPlot(seu.obj1, group.by = cluster.ident,
                     features = c(signature.names, paste0(signature.names,"_kNN")), 
                     cols = cluster_col, ncol = 2, pt.size = 0, raster = raster) &
      geom_boxplot(width = 0.1, color = "black", alpha = 0.3) &
      theme(text = element_text(size = 7), 
            axis.text = element_text(size = 7),
            axis.title.x = element_text(size = 7), 
            axis.title.y = element_text(size = 7))
  }
  
  if(isTRUE(print)){
    print(fp)
  }
    
  # save plot
  pdf(file.path(save.loc, paste0(plot.name, ".plotUCell.pdf")), width = plot.width, height = plot.height)
  print(fp)
  if(isTRUE(vln)){
    print(vln_p)
  }
  dev.off()
  
  return(seu.obj1)
}
