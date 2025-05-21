limma_cluster_DEGs <- function(
    seurat_obj,
    assay = "RNA",
    comparison_type = c("one_vs_rest", "pairwise"),
    cluster_column = NULL,
    sample_column = NULL,
    pseudobulk = FALSE,
    min_cells = 10
) {
  library(Seurat)
  library(limma)
  library(Matrix)
  comparison_type <- match.arg(comparison_type)
  
  # Get cluster and sample info
  meta <- seurat_obj@meta.data
  clusters <- if (!is.null(cluster_column)) meta[[cluster_column]] else Idents(seurat_obj)
  samples <- if (!is.null(sample_column)) meta[[sample_column]] else rep("sample1", ncol(seurat_obj))
  
  # Filter low cell clusters
  cluster_counts <- table(clusters)
  valid_clusters <- names(cluster_counts[cluster_counts >= min_cells])
  keep_cells <- which(clusters %in% valid_clusters)
  
  clusters <- factor(clusters[keep_cells], levels = valid_clusters)
  samples <- samples[keep_cells]
  
  # Get raw counts or normalized expression
  if (pseudobulk) {
    # Use raw counts (for summing across cells)
    expr_matrix <- as.matrix(GetAssayData(seurat_obj, assay = assay, slot = "counts")[, keep_cells])
    
    # Build pseudobulk matrix: sum counts per cluster-sample
    comb_labels <- paste(clusters, samples, sep = "_")
    pseudobulk_mat <- rowsum(t(expr_matrix), group = comb_labels)
    pseudobulk_mat <- t(pseudobulk_mat)
    
    # Get design matrix
    cluster_labels <- sapply(strsplit(colnames(pseudobulk_mat), "_"), `[`, 1)
    group <- factor(cluster_labels)
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    
  } else {
    # Single-cell mode (each cell is a sample)
    expr_matrix <- as.matrix(GetAssayData(seurat_obj, assay = assay, slot = "data")[, keep_cells])
    group <- factor(clusters)
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
  }
  
  # Make contrasts
  if (comparison_type == "one_vs_rest") {
    contrast_strings <- sapply(levels(group), function(cl) {
      others <- setdiff(levels(group), cl)
      paste0(cl, " - (", paste(others, collapse = " + "), ")/", length(others))
    })
    names(contrast_strings) <- paste0(levels(group), "_vs_rest")
    
  } else {
    combs <- combn(levels(group), 2, simplify = FALSE)
    contrast_strings <- sapply(combs, function(pair) {
      paste0(pair[1], " - ", pair[2])
    })
    names(contrast_strings) <- paste0(sapply(combs, paste, collapse = "_vs_"))
  }
  
  contrast.matrix <- makeContrasts(contrasts = contrast_strings, levels = design)
  
  # Run limma
  fit <- lmFit(expr_matrix, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Extract DEGs
  degs_list <- lapply(colnames(contrast.matrix), function(cn) {
    topTable(fit2, coef = cn, number = Inf, adjust.method = "fdr")
  })
  names(degs_list) <- colnames(contrast.matrix)
  
  return(degs_list)
}
