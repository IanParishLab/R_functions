plot_Azimuth <- function(obj, level = 1, suffix, save_name, save.loc){
  
  library(ggpubr)
  library(tidyverse)
  
  ident = paste0("predicted.celltype.l",level,".",suffix)
  score = paste0("predicted.celltype.l",level,".score.",suffix)
  map = paste0("mapping.score.",suffix)
  
  Idents(obj) <- ident
  cell_types <- unique(obj[[ident]]) %>% unlist %>% unname
  pred_pl <- setNames(lapply(cell_types, function(cell_type){
    p1 <- DimPlot(obj,pt.size = 3, sizes.highlight = 3, 
                  cells.highlight = WhichCells(obj, idents = cell_type), 
                  order=T,repel=T, raster = TRUE) +
      NoAxes() + NoLegend() +
      ggtitle(cell_type)
    p2 <- gghistogram(obj[[map]][obj[[ident]] == cell_type]) + 
      ggtitle(map)
    p3 <- gghistogram(obj[[score]][obj[[ident]] == cell_type]) + 
      ggtitle(score)
    
    return(p1 + p2 + p3 & 
             theme(aspect.ratio = 1, title = element_text(size = 10)))
  }), cell_types)
  
  pred_pl <- pred_pl[sort(cell_types)]
 
  ifelse(!dir.exists(file.path(save.loc, "plots")),
         dir.create(file.path(save.loc, "plots"), recursive = TRUE), paste0(save.loc," directory exists"))

  pdf(file.path(save.loc, 'plots', paste0(save_name, ".",ident,".pdf")), width = 10, height = 4)
  print(pred_pl)
  dev.off()
}
