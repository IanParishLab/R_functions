plot_Overview <- function(seu.obj, reduction.name = "umap.scvi", ident = "PatientID", alpha = 0.4,
                          hash.ident = "sampleSource", other.hash.ident = NULL,
                          res = "RNA_snn_res.1", save.loc = "plots", plot.name = "",
                          HTO_col, capture_col, cluster_col, 
                          plot.width = 20, plot.height = 10, ncol = 3, nrow = 3, raster = FALSE, print = FALSE) {
  
  
  # seu.obj = tmp
  # reduction.name = "umap"
  # ident = "Capture"
  # alpha = 0.4
  # hash.ident = "sampleSource_simple"
  # other.hash.ident = "sampleSourceSuperset_simple"
  # res = "RNA_snn_res.1"
  # save.loc = save.loc
  # plot.name = obj_name
  # HTO_col <- HTO_col_simple
  # capture_col <- data_col
  # # cluster_col
  # plot.width = 8
  # plot.height = 8
  # ncol = 4
  # nrow = 3
  # raster = FALSE
  # print = FALSE
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(ggpubr)
    require(tidyverse)
    require(patchwork)
  })
  
  # Make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # get color palette if not provided
  # if (is.null(capture_col) && exists("capture_col", envir = .GlobalEnv)) {
  #   capture_col <- get("capture_col", envir = .GlobalEnv)
  # }
  if (is.null(cluster_col) && exists("cluster_col", envir = .GlobalEnv)) {
    cluster_col <- get("cluster_col", envir = .GlobalEnv)
  }
  # if (is.null(HTO_col) && exists("HTO_col", envir = .GlobalEnv)) {
  #   HTO_col <- get("HTO_col", envir = .GlobalEnv)
  # }
  
  # Shared theme
  theme <- theme(aspect.ratio = 1, 
                 text = element_text(size = 10), 
                 plot.title = element_text(hjust = 0.5, size = 10)) &
    NoAxes() & NoLegend()
  
  # Participant-level overview
  p1 <- DimPlot(seu.obj, reduction = reduction.name, group.by = ident, alpha = alpha,
                ncol = ncol, order = TRUE, repel = TRUE, cols = capture_col, raster = raster) + 
    ggtitle("Pooled Participant Data")
  p2 <- DimPlot(seu.obj, reduction = reduction.name, group.by = ident, split.by = ident, alpha = alpha,
                ncol = ncol, order = TRUE, repel = TRUE, cols = capture_col, raster = raster)
  p_plot <- p1 + p2 + plot_layout(widths = c(1, 4)) & theme
  
  # Cluster overview
  sc_plot <- DimPlot(seu.obj, reduction = reduction.name, group.by = res, alpha = alpha, 
                     order = TRUE, label = TRUE, label.box = TRUE, label.size = 2, repel = TRUE, 
                     cols = cluster_col, raster = raster)
  
  # Hash-level plots
  q <- list()
  hash <- if (!is.null(hash.ident)) {
    intersect(levels(seu.obj@meta.data[[hash.ident]]), seu.obj@meta.data[[hash.ident]])
  } else NULL
  
  if (!is.null(hash)) {
    Idents(seu.obj) <- hash.ident
    for (HTO in hash) {
      q[[HTO]] <- DimPlot(seu.obj, reduction = reduction.name, group.by = hash.ident, 
                          order = TRUE, repel = TRUE, alpha = alpha,
                          cells.highlight = list(WhichCells(seu.obj, idents = HTO)),
                          cols.highlight = list(HTO_col[[HTO]]),
                          sizes.highlight = 0.3,
                          cols = 'grey', raster = raster) +
        ggtitle(HTO) + theme
    }
  }
  
  if (!is.null(other.hash.ident)) {
    other_hash <- setdiff(
      intersect(levels(seu.obj@meta.data[[other.hash.ident]]), seu.obj@meta.data[[other.hash.ident]]),
      hash
    )
    message("[MSG] Also plotting other hashes: ", paste0(other_hash, collapse = ", "))
    
    Idents(seu.obj) <- other.hash.ident
    for (HTO in other_hash) {
      q[[HTO]] <- DimPlot(seu.obj, reduction = reduction.name, group.by = other.hash.ident, 
                          order = TRUE, repel = TRUE, alpha = alpha,
                          cells.highlight = list(WhichCells(seu.obj, idents = HTO)),
                          cols.highlight = list(HTO_col[[HTO]]),
                          sizes.highlight = 0.3,
                          cols = 'grey', raster = raster) +
        ggtitle(HTO) + theme
    }
    hash_plot <- ggarrange(plotlist = q[c(hash, other_hash)], ncol = ncol, nrow = nrow)
  } else if (!is.null(hash)) {
    hash_plot <- ggarrange(plotlist = q, ncol = ncol, nrow = nrow)
  }
  
  cluster_plot <- if (!is.null(hash)) {
    sc_plot + hash_plot + plot_layout(widths = c(1, 4)) & theme
  } else {
    sc_plot & theme
  }
  
  # Optionally display
  if (print) {
    print(p_plot)
    print(cluster_plot)
  }
  
  # Save to PDF
  png(file.path(save.loc, "plots", paste0(plot.name, ".plotOverview.png")), width = plot.width, height = plot.height, units = "in", res = 300)
  print(p_plot | cluster_plot)
  dev.off()
}
