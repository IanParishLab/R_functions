alt_plot_Composition <- function(seu.obj,
                             ident = "PatientID",
                             res = "seuratClusters",
                             hash.ident = "sampleSource",
                             other.hash.ident = "sampleSourceSuperset",
                             save.loc = "./plots",
                             plot.name = seu.obj@project.name,
                             cluster_col = NULL,
                             HTO_col = NULL,
                             plot.width = 8,
                             plot.height = 8,
                             ncol = 3,
                             nrow = 2,
                             print = FALSE) {
  suppressPackageStartupMessages({
    library(dittoSeq)
    library(Seurat)
    library(tidyverse)
  })
  
  # Validate inputs
  if (!dir.exists(save.loc)) dir.create(save.loc, recursive = TRUE)
  if (!all(c(res, ident) %in% colnames(seu.obj@meta.data))) {
    stop("[ERROR] Missing required metadata: res or ident")
  }
  
  # Set theme
  base_theme <- theme(text = element_text(size = 7),
                      axis.text.x = element_text(angle = 90, size = 7, hjust = 1))
  
  # Helper function to build bar plots
  make_bar_plots <- function(var_group, var_split, color.panel, group.label, split.label) {
    if (!all(c(var_group, var_split) %in% colnames(seu.obj@meta.data))) return(NULL)
    
    # levels_x <- levels(factor(seu.obj@meta.data[[var_split]]))
    
    p1 <- dittoBarPlot(seu.obj, var = var_group, group.by = var_split, xlab = group.label,
                       retain.factor.levels = TRUE, color.panel = color.panel) +
      ggtitle(paste0(var_group, " per ", var_split)) +
      base_theme +
      scale_x_discrete(limits = levels_x) +
      guides(color = guide_legend(ncol = 2)) +
      NoLegend()
    
    p2 <- dittoBarPlot(seu.obj, group.by = var_group, var = var_split, xlab = split.label,
                       retain.factor.levels = TRUE, color.panel = color.panel) +
      ggtitle(paste0(var_split, " per ", var_group)) +
      base_theme +
      guides(color = guide_legend(ncol = 2)) +
      NoLegend()
    
    # Add patient-split panels
    p1 <- p1 + dittoBarPlot(seu.obj, var = var_group, group.by = var_split, split.by = ident,
                            xlab = group.label, retain.factor.levels = TRUE, color.panel = color.panel,
                            split.nrow = nrow, split.ncol = ncol) +
      ggtitle("Split by patient") + base_theme +
      scale_x_discrete(limits = levels_x) +
      guides(color = guide_legend(ncol = 2))
    
    p2 <- p2 + dittoBarPlot(seu.obj, group.by = var_group, var = var_split, split.by = ident,
                            xlab = split.label, retain.factor.levels = TRUE, color.panel = color.panel,
                            split.nrow = nrow, split.ncol = ncol) +
      ggtitle("Split by patient") + base_theme +
      guides(color = guide_legend(ncol = 2))
    
    return(p1 / p2)
  }
  
  # ---- Main Composition Plots ----
  message("[INFO] Creating cluster composition plot...")
  cluster_composition <- make_bar_plots(res, hash.ident, cluster_col, "Sample source", "Cluster")
  
  # ---- Superset Composition Plots ----
  cluster_composition_superset <- NULL
  if (!is.null(other.hash.ident) && other.hash.ident %in% colnames(seu.obj@meta.data)) {
    message("[INFO] Creating superset composition plot...")
    cluster_composition_superset <- make_bar_plots(res, other.hash.ident, cluster_col, "Sample source", "Cluster")
  }
  
  # ---- Save to PDF ----
  message("[INFO] Saving to PDF...")
  pdf(file.path(save.loc, paste0(plot.name, ".plotComposition.pdf")),
      width = plot.width, height = plot.height)
  if (!is.null(cluster_composition)) print(cluster_composition)
  if (!is.null(cluster_composition_superset)) print(cluster_composition_superset)
  dev.off()
  
  # ---- Return if requested ----
  if (isTRUE(print)) {
    if (!is.null(cluster_composition_superset)) {
      return(cluster_composition_superset)
    } else {
      return(cluster_composition)
    }
  }
  
  invisible(NULL)
}
