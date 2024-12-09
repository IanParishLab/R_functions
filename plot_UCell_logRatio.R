plot_UCell_logRatio <- function(so,group1,group2,group.ident = "PatientID", save.loc){
  suppressPackageStartupMessages({
    library(UCell)
    library(Seurat)
    library(tidyverse)
    library(RColorBrewer)
  })
  
  # get test name
  test_name = paste0(c(group1,group2), collapse = "_") %>% gsub("_UCell","",.)
  
  # calculate ratios
  so[[test_name]] <- so[[group1]]/so[[group2]]
  so <- SmoothKNN(so,
                  signature.names = test_name,
                  reduction = "pca")
  
  # calculate log ratios 
  knn_name = paste0(test_name,"_kNN")
  log_name = paste0("log_",test_name)
  log_knn_name = paste0("log_",test_name,"kNN")
  
  # define cells to plot
  so@meta.data[c(test_name, knn_name)] %>% log ->
    so@meta.data[c(log_name, log_knn_name)]
  
  plot_cells <- c(which(is.na(so[[log_name]])), 
                  which(is.infinite(so[[log_knn_name]]))) 
  plot_cells <- colnames(so)[-plot_cells]
  
  # make and save plot
  pdf(file.path(save.loc,paste0(test_name,'.FeaturePlot.pdf')), width = 15, height = 5)
  FeaturePlot(so, cells = plot_cells, split.by = group.ident, reduction = "umap.scvi", order = TRUE, repel = TRUE,
              features = log_knn_name, ncol = 6) &
    NoAxes() &
    theme(aspect.ratio = 1, text = element_text(size = 7),
          axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  dev.off()
  
  return(so)
}
