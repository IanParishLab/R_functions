plot_Composition <- function(seu.obj, ident = "PatientID", res = "seuratClusters", 
                             hash.ident = "sampleSource", other.hash.ident = NULL,
                             save.loc = "plots", plot.name = "test",
                             cluster_col, HTO_col,
                             plot.width = 8, plot.height = 8, ncol = 3, nrow = 2,
                             print = FALSE,...) {
  
  suppressPackageStartupMessages({
    require(dittoSeq)
    require(Seurat)
    require(tidyverse)
    require(patchwork)
  })
  
  # make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # get color palette if not provided
  if (is.null(cluster_col) && exists("cluster_col", envir = .GlobalEnv)) {
    cluster_col <- get("cluster_col", envir = .GlobalEnv)
  }
  if (is.null(HTO_col) && exists("HTO_col", envir = .GlobalEnv)) {
    HTO_col <- get("HTO_col", envir = .GlobalEnv)
  }
  
  theme_custom <- theme(text = element_text(size = 7), 
                        axis.text.x = element_text(angle = 90, size = 7, hjust = 1))
  
  # Helper for composition (returns a patchwork object)
  make_composition <- function(by1, by2, col1, col2, split.by = ident, split.nrow = nrow, split.ncol = ncol) {
    a <- dittoBarPlot(seu.obj, var = by2, group.by = by1, xlab = "Sample source", 
                      retain.factor.levels = TRUE, color.panel = col1) +
      ggtitle("Cell cluster proportions per sample type") +
      theme_custom +
      guides(color = guide_legend(ncol = 2)) +
      NoLegend()
    a_split <- dittoBarPlot(seu.obj, var = by2, group.by = by1, split.by = split.by, 
                            xlab = "Sample source", retain.factor.levels = TRUE, color.panel = col1,
                            split.nrow = split.nrow, split.ncol = split.ncol) +
      ggtitle("Split by patient") + theme_custom + guides(color = guide_legend(ncol = 2))
    a <- a + a_split
    
    b <- dittoBarPlot(seu.obj, group.by = by2, var = by1, xlab = "Cluster", 
                      retain.factor.levels = TRUE, color.panel = col2) +
      ggtitle("Sample type proportions per cell cluster") +
      theme_custom +
      guides(color = guide_legend(ncol = 2)) +
      NoLegend()
    b_split <- dittoBarPlot(seu.obj, var = by1, group.by = by2, split.by = split.by, 
                            xlab = "Cluster", retain.factor.levels = TRUE, color.panel = col2,
                            split.nrow = split.nrow, split.ncol = split.ncol) +
      ggtitle("Split by patient") + guides(color = guide_legend(ncol = 2)) + theme_custom
    b <- b + b_split
    
    a / b
  }
  
  message("Generating cluster composition plot(s)...")
  cluster_composition <- make_composition(hash.ident, res, cluster_col, HTO_col)
  if (!is.null(other.hash.ident)) {
    cluster_composition_superset <- make_composition(other.hash.ident, res, cluster_col, HTO_col)
  }
  
  png(file.path(save.loc, "plots", paste0(plot.name, ".plotComposition.png")), width = plot.width, height = plot.height,
      units = "in", res = 300)
  print(cluster_composition)
  if (!is.null(other.hash.ident)) print(cluster_composition|cluster_composition_superset)
  dev.off()
  
  if (isTRUE(print)) {
    message("Print plot...")
    return(if (!is.null(other.hash.ident)) (cluster_composition | cluster_composition_superset) else cluster_composition)
  }
}
