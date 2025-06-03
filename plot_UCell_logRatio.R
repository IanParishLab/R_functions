plot_UCell_logRatio <- function(so, group1, group2, reduction = 'umap.scvi', group.ident = NULL, 
                                save.loc = "plots", alpha = 0.4, raster = FALSE, nrow, ncol, width, height){
  suppressPackageStartupMessages({
    library(UCell)
    library(Seurat)
    library(tidyverse)
    library(RColorBrewer)
  })
  
  # so = tmp
  # group1 = 'Tolerance_UCell'
  # group2 = 'Exhaustion_UCell'
  # reduction = "umap"
  # group.ident = 'PatientID'
  # save.loc = save.loc
  
  # make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # get test name
  test_name = paste0(c(group1,group2), collapse = "_") %>% gsub("_UCell","",.)
  
  # calculate ratios
  so[[test_name]] <- so[[group1]]/so[[group2]]
  so <- SmoothKNN(so, signature.names = test_name, reduction = "pca")
  
  # calculate log ratios 
  knn_name = paste0(test_name,"_kNN")
  log_name = paste0("log_",test_name)
  log_knn_name = paste0("log_",test_name,"_kNN")
  
  so@meta.data[c(test_name, knn_name)] %>% log ->
    so@meta.data[c(log_name, log_knn_name)]
  
  # define cells to plot
  ind1 <- which(is.na(so[[log_name]]))
  ind2 <- which(is.infinite(so[[log_knn_name]] %>% unlist))
  rm_cells <- unique(c(ind1,ind2))
  plot_cells <- colnames(so)[-rm_cells]
  
  # make and save plot
  p <- FeaturePlot(so, cells = plot_cells, 
                   reduction = reduction, order = TRUE, repel = TRUE, alpha = alpha, raster = raster,
                   features = log_knn_name)
  q <- FeaturePlot(so, cells = plot_cells, split.by = group.ident, 
                   reduction = reduction, order = TRUE, repel = TRUE, alpha = alpha, raster = raster,
                   features = log_knn_name, ncol = ncol) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  plot <- (p + q) &
    NoAxes() &
    theme(aspect.ratio = 1, text = element_text(size = 7),
          axis.title.x = element_text(size = 7), 
          axis.title.y = element_text(size = 7)) 
  
  pdf(file.path(save.loc,"plots",paste0(test_name,'.FeaturePlot.pdf')), width = width, height = height)
  print(plot + plot_layout(width = c(0.1,0.4), nrow = nrow, ncol=ncol))
  dev.off()
  
  return(so)
}


