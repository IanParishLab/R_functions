process_celltypist_results <- function(seu, 
                                       pred, 
                                       encyclopedia_file = "/researchers/nicole.saw/references/celltypist/encyclopedia_table.xlsx", 
                                       reduction = 'umap',
                                       save.loc, 
                                       sample_name) {
  
  # 2. Read CellTypist encyclopedia
  invisible(
    enc <- xlsx::read.xlsx(file = encyclopedia_file, sheetIndex = 1)
  )
  
  # 3. Create summary table for plotting
  df <- data.frame(table(pred$majority_voting)) %>% 
    setNames(c("Low.hierarchy.cell.types", "count")) %>% 
    left_join(., enc[, 1:2], by = c("Low.hierarchy.cell.types" = "Low.hierarchy.cell.types")) %>% 
    select(High.hierarchy.cell.types, Low.hierarchy.cell.types, count) %>% 
    mutate(count = ifelse(is.na(count), 0, count))
  
  # 4. Barplot
  bar_plot <- ggbarplot(df, sort.val = "desc",
                        x = "Low.hierarchy.cell.types", y = "count", 
                        fill = "High.hierarchy.cell.types",
                        label = TRUE, lab.pos = "out", lab.size = 3, 
                        lab.vjust = 0.5, lab.hjust = -0.5) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 7))
  
  # 5. Add CellTypist labels to Seurat object
  colnames(pred)[which(colnames(pred) == "majority_voting")] <- "Low.hierarchy.cell.types"
  keep_row_names <- rownames(pred)
  pred <- left_join(pred, enc[, 1:2], by = "Low.hierarchy.cell.types")
  rownames(pred) <- keep_row_names
  seu <- AddMetaData(seu, pred)
  
  # 6. DimPlot of High hierarchy labels
  umap_plot <- DimPlot(seu, group.by = "High.hierarchy.cell.types", reduction = reduction)
  
  png(here(save.loc, "plots", paste0(sample_name, ".celltypist_results.png")), width = 10, height = 10, units = 'in', res = 300)
  print(bar_plot)
  print(umap_plot)
  dev.off()
  
  return(seu)
}