plot_Overview <- function(seu.obj, reduction.name = "umap.scvi", ident = "PatientID", alpha = 0.4,
                          hash.ident = "sampleSource", other.hash.ident = NULL,
                          res = "RNA_snn_res.1", save.loc = "./plots", plot.name = "",
                          HTO_col = HTO_col, capture_col = capture_col, cluster_col = cluster_col, 
                          plot.width = 12, plot.height = 10, ncol=3, nrow=3, raster = FALSE, print = FALSE){
  
  # seu.obj = scvi
  # reduction.name = "umap.scvi"
  # ident = "PatientID"
  # hash.ident = "sampleSource"
  # other.hash.ident = "sampleSourceSuperset"
  # res = "seuratClusters"
  # save.loc = save.loc
  # plot.name = int_obj_name
  # plot.width = 10
  # plot.height = 10
  # ncol=3
  # nrow=3
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(ggpubr)
    require(tidyverse)
    require(patchwork)
  })
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  # set theme
  theme = theme(aspect.ratio = 1, 
                text = element_text(size = 10), 
                plot.title = element_text(hjust = 0.5, size = 10)) &
    NoAxes() &
    NoLegend()
  
  # patient plot
  p1 <- DimPlot(seu.obj, reduction = reduction.name, group.by = ident, alpha = 0.4,
                ncol = ncol, order = TRUE, repel = TRUE, cols = capture_col, raster = raster) + 
    ggtitle("Pooled Participant Data")
  p2 <- DimPlot(seu.obj, reduction = reduction.name, group.by = ident, split.by = ident,
                ncol = ncol, order = TRUE, repel = TRUE, cols = capture_col, raster = raster) + 
    ggtitle("")
  
  p_plot <- p1 + p2 + plot_layout(widths = c(1,4)) & 
    theme
  
  # seurat cluster plot
  sc_plot <- DimPlot(seu.obj, reduction = reduction.name, group.by = res, alpha = 0.4, 
                     order = TRUE, label = TRUE, label.box = TRUE, label.size = 2, repel = TRUE, 
                     cols = cluster_col, raster = raster)  
  
  # hash plot
  if(!is.null(hash.ident)){
    hash <- intersect(levels(seu.obj@meta.data[[hash.ident]]),
                      seu.obj@meta.data[[hash.ident]])
    Idents(seu.obj) <- hash.ident
    q <- list()
    for(HTO in hash){
      q[[HTO]] <- DimPlot(seu.obj, reduction = reduction.name, group.by = hash.ident, 
                          order =TRUE, repel = TRUE, alpha = 0.4,
                          cells.highlight = list(WhichCells(seu.obj, idents = HTO)),
                          cols.highlight = list(HTO_col[HTO]),
                          sizes.highlight = 0.3,
                          cols='grey', alpha = 0.4,
                          raster = raster) +
        ggtitle(HTO) + 
        theme
    }
  }
  
  # other hash ident
  if(!is.null(other.hash.ident)){
    other_hash <- intersect(levels(seu.obj@meta.data[[other.hash.ident]]),
                            seu.obj@meta.data[[other.hash.ident]])
    other_hash <- other_hash[!(other_hash %in% hash)]
    print(paste0("[MSG] Also plotting other hashes ", paste0(other_hash,collapse = ",")))
    
    Idents(seu.obj) <- other.hash.ident
    for(HTO in other_hash){
      q[[HTO]] <- DimPlot(seu.obj, reduction = reduction.name, group.by = other.hash.ident, 
                          order =TRUE, repel = TRUE, alpha = 0.4,
                          cells.highlight = list(WhichCells(seu.obj, idents = HTO)),
                          cols.highlight = list(HTO_col[HTO]),
                          sizes.highlight = 0.3,
                          cols='grey', 
                          raster = raster) +
        ggtitle(HTO) + 
        theme
    }
    hash_plot <- ggarrange(plotlist = q[c(hash,other_hash)], ncol = ncol, nrow = nrow)
  } else {
    hash_plot <- ggarrange(plotlist = q, ncol = ncol, nrow = nrow)
  }
  
  if(!is.null(hash.ident)){
    cluster_plot <- sc_plot + hash_plot + plot_layout(widths = c(1,4)) & theme 
  } else {
    cluster_plot <- sc_plot & theme
  }
  
  if(isTRUE(print)){
    print(p_plot)
    print(cluster_plot)
  }
    
  # save plot
  pdf(file.path(save.loc, paste0(plot.name,".plotOverview.pdf")), width = plot.width, height = plot.height)
  print(p_plot)
  print(cluster_plot)
  dev.off()
  
}
