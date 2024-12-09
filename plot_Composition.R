plot_Composition <- function(seu.obj, ident = "PatientID", res  = "seuratClusters", 
                             hash.ident = "sampleSource", other.hash.ident = "sampleSourceSuperset",
                             save.loc = "./plots", plot.name,
                             cluster_col = cluster_col, HTO_col = HTO_col,
                             plot.width = 8, plot.height = 8, ncol = 3, nrow = 2){
  
  require(dittoSeq)
  require(Seurat)
  require(tidyverse)
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  a <- dittoBarPlot(seu.obj, var = res, group.by = hash.ident, xlab = "Sample source", 
                    retain.factor.levels = TRUE, color.panel = cluster_col) + 
    ggtitle("Cell cluster proportions per sample type") +
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) + 
    scale_x_discrete(limits= levels(seu.obj@meta.data[[hash.ident]])) + 
    NoLegend()
  b <- dittoBarPlot(seu.obj, group.by = res, var = hash.ident, xlab = "Cluster", 
                    retain.factor.levels = TRUE, color.panel = HTO_col) + 
    ggtitle("Sample type proportions per cell cluster") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) + 
    NoLegend()
  
  a <- a + dittoBarPlot(seu.obj, var = res, group.by = hash.ident, split.by = ident, xlab = "Sample source", 
                        retain.factor.levels = TRUE, color.panel = cluster_col, 
                        split.nrow = nrow,split.ncol = ncol) + 
    ggtitle("Split by patient") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) +
    scale_x_discrete(limits= levels(seu.obj@meta.data[[hash.ident]]))
  
  b <- b + dittoBarPlot(seu.obj, var = hash.ident, group.by = res, split.by = ident, xlab = "Cluster", 
                        retain.factor.levels = TRUE, color.panel = HTO_col, 
                        split.nrow = nrow,split.ncol = ncol) + 
    ggtitle("Split by patient") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1))
  
  cluster_composition <- a/b
  
  # superset
  a <- dittoBarPlot(seu.obj, var = res, group.by = other.hash.ident, xlab = "Sample source",
                    retain.factor.levels = TRUE, color.panel = cluster_col) + 
    ggtitle("Cell cluster proportions per sample type") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) + 
    scale_x_discrete(limits= levels(seu.obj@meta.data[[other.hash.ident]])) + 
    NoLegend()
  b <- dittoBarPlot(seu.obj, group.by = res, var = other.hash.ident, xlab = "Cluster", 
                    retain.factor.levels = TRUE, color.panel = HTO_col) + 
    ggtitle("Sample type proportions per cell cluster") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) + 
    NoLegend()
  
  a <- a + dittoBarPlot(seu.obj, var = res, group.by = other.hash.ident, xlab = "Sample source",
                        retain.factor.levels = TRUE, color.panel = cluster_col, 
                        split.by = ident, split.nrow = nrow,split.ncol = ncol) +
    ggtitle("Split by patient") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) + 
    scale_x_discrete(limits= levels(seu.obj@meta.data[[other.hash.ident]]))
  
  b <- b + dittoBarPlot(seu.obj, group.by = res, var = other.hash.ident, xlab = "Cluster", 
                        retain.factor.levels = TRUE, color.panel = HTO_col,
                        split.by = ident, split.nrow = nrow,split.ncol = ncol) +
    ggtitle("Split by patient") + 
    theme(text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 7, hjust = 1)) 
  
  cluster_composition_superset <- a/b
  
  # save plot
  pdf(file.path(save.loc, paste0(plot.name,".plotComposition.pdf")), width = plot.width, height = plot.height)
  print(cluster_composition)
  print(cluster_composition_superset)
  dev.off()
  
}
