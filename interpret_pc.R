interpret_pc <- function(seurat_obj, pc = 3, nfeatures = 20, plot = TRUE) {
  # Check PCA is available
  if (!"pca" %in% names(seurat_obj@reductions)) {
    stop("No PCA reduction found. Run RunPCA() first.")
  }
  
  # Extract loadings (gene weights)
  loadings <- seurat_obj@reductions$pca@feature.loadings
  
  if (pc > ncol(loadings)) {
    stop(paste("PC", pc, "not found. PCA has only", ncol(loadings), "components."))
  }
  
  # Sort by absolute contribution
  pc_loadings <- loadings[, pc, drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::arrange(desc(abs(.data[[1]])))
  
  # Top positive and negative genes
  top_pos <- head(pc_loadings[order(-pc_loadings[,2]), ], nfeatures)
  top_neg <- head(pc_loadings[order(pc_loadings[,2]), ], nfeatures)
  
  # Combine
  results <- list(
    top_positive = top_pos,
    top_negative = top_neg,
    all_loadings = pc_loadings
  )
  
  # Optional plot
  if (plot) {
    plot_data <- rbind(
      dplyr::mutate(top_pos, direction = "positive"),
      dplyr::mutate(top_neg, direction = "negative")
    )
    
    p <- ggplot(plot_data, aes(x = reorder(gene, !!sym(names(plot_data)[2])), 
                               y = !!sym(names(plot_data)[2]), 
                               fill = direction)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste0("Top loadings for PC", pc),
           y = "Loading weight", x = "Gene") +
      theme_minimal()
    
    print(p)
  }
  
  return(results)
}
