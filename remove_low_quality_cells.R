remove_low_quality_cells <- function(so, assay="RNA", genome, low_abundant_genes_limit = -5, libsize_drop_nmads = 2, mito_drop_nmads = 3){
  suppressWarnings(suppressPackageStartupMessages({
    library(ruvIIInb) 
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(scuttle)
    library(edgeR)
    library(celldex)
    library(hrbrthemes)
    library(tidyverse)
    library(ggplot2)
    library(scMerge)
    library(Seurat)
    library(randomcoloR)
    library(dittoSeq)
    library(pheatmap)
    library(gridExtra)
    library(igraph)
    library(DelayedArray)
  }))
  
  print("[MSG] Creating SingleCellExperiment...")
  # so is Seurat Object
  sce <- SingleCellExperiment(assays = list(counts = so[[assay]]$counts),
                              colData = so@meta.data)
  
  # 2.2.1 Filtering low quality cells
  # The cell-level quality control metrics such as the number of genes that have non-zero counts and the percentage of counts that comes from Mitochondrial genes for each cell can be computed and added to the SingleCellExperiment object as follows:
  if (genome == "mouse"){
    sce <- addPerCellQCMetrics(x = sce,subsets=list(Mito=grep("Mt-",rownames(sce))))
  } else if (genome == "human"){
    sce <- addPerCellQCMetrics(x = sce,subsets=list(Mito=grep("MT-",rownames(sce))))
  } else {
    stop("genome not included - choose only human or mouse.")
  }
  
  libsize_drop <- isOutlier(
    metric = sce$total,
    nmads = libsize_drop_nmads,
    type = "lower",
    log = TRUE)
  colData(sce)$libsize_drop<-libsize_drop
  
  mito_drop <- isOutlier(
    metric = colData(sce)$subsets_Mito_percent,
    nmads = mito_drop_nmads,
    type = "higher")
  colData(sce)$mito_drop<-mito_drop
  
  # 2.2.2 Filtering low abundant genes
  # Similarly, we recommend computing gene-level quality control metrics such as the mean count across all cells for each gene, and the percentage of cells with non-zero counts for each gene. The code below adds these measures to the SingleCellExperiment object and flags genes with a low mean cell count (see Figure 2.3).
  sce <- addPerFeatureQCMetrics(x = sce)
  
  #Remove genes with zero counts for each gene
  sce <- subset(sce, rowData(sce)$mean > 0 )
  
  # detect low abundant genes
  lowcount_drop <- log(rowData(sce)$mean) < low_abundant_genes_limit
  
  print("[MSG] Subsetting SingleCellExperiment to filter out low quality cells...")
  # subset 
  print("[MSG] Number of low counts genes:")
  print(table(lowcount_drop))
  print("[MSG] Number of cell-level total counts with number of median absolute deviations (nMADs = 2) away from the median value:")
  print(table(libsize_drop))
  print("[MSG] Number of cell-level mitochondrial percentage with number of median absolute deviations (nMADs = 3) away from the median value:")
  print(table(mito_drop))
  sce <- sce[!(lowcount_drop), !(libsize_drop | mito_drop)]

  # logNormCounts() function from scuttle will compute a log-transformed normalized expression matrix and store it as another assay.
  sce <- scuttle::logNormCounts(sce)

  seu.obj <- as.Seurat(sce, counts = "counts", data = "logcounts", project = "SeuratExpt")
  print("[MSG] returning Seurat Object... ")
  return(seu.obj)
  
}
