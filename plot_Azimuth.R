plot_Azimuth <- function(seu.obj, level = 1, suffix, save_name, save.loc){
  suppressPackageStartupMessages({
    library(ggpubr)
    library(tidyverse)
  })

# make output directory
dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)


  ident <- paste0("predicted.celltype.l", level, ".", suffix)
  score <- paste0("predicted.celltype.l", level, ".score.", suffix)
  map <- paste0("mapping.score.", suffix)

  Idents(seu.obj) <- ident
  cell_types <- unique(seu.obj[[ident]]) %>% unlist() %>% unname()

  message("[MSG] Generating Azimuth plots for: ", ident)

  pred_pl <- setNames(lapply(cell_types, function(cell_type) {
    p1 <- DimPlot(
      seu.obj,
      pt.size = 3,
      sizes.highlight = 3,
      cells.highlight = WhichCells(seu.obj, idents = cell_type),
      order = TRUE,
      repel = TRUE,
      raster = TRUE
    ) + NoAxes() + NoLegend() + ggtitle(cell_type)

    p2 <- gghistogram(seu.obj[[map]][seu.obj[[ident]] == cell_type]) + ggtitle(map)
    p3 <- gghistogram(seu.obj[[score]][seu.obj[[ident]] == cell_type]) + ggtitle(score)

    return(p1 + p2 + p3 & theme(aspect.ratio = 1, title = element_text(size = 10)))
  }), cell_types)

  pred_pl <- pred_pl[sort(cell_types)]

  pdf_path <- file.path(plot_dir, paste0(save_name, "_", ident, ".pdf"))
  message("[MSG] Saving Azimuth plots to: ", pdf_path)
  pdf(pdf_path, width = 10, height = 4)
  print(pred_pl)
  dev.off()
}
