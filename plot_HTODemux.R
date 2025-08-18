plot_HTODemux <- function(seu.obj, 
                          sample_name, 
                          save.loc = "plots", 
                          HTO_AssayName = "HTO", 
                          RNA_AssayName = "RNA", 
                          hash_class_column = "HTO_classification.global", 
                          hash.ident = "HTO_maxID",
                          cols = NULL,
                          ncol = 4, nrow = 3,
                          plot.width = 12, plot.height = 6){
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
    require(ggpubr)
  })
  
  # make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)

  # print hash class  
  print(paste0("[MSG] Checking hash classification..."))
  if(isFALSE(all(seu.obj[[hash_class_column]] == "Singlet"))){
    print("[WARNING] not all droplets are hash singlets")
  }
  
  # make plots
  DefaultAssay(seu.obj) <- HTO_AssayName
  
  df1 <- as.data.frame(t(seu.obj[[HTO_AssayName]]@counts))
  distr_plot <- df1 %>%
    pivot_longer(cols = names(df1)) %>%
    mutate(logged = log(value + 1)) %>%
    ggplot(aes(x = logged)) +
    xlab("log(counts)") +
    xlim(0.1,8) +
    geom_density(adjust = 3) +
    facet_wrap(~name, scales = "fixed", ncol = 4) &
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1)
  
  
  # RidgePlot
  Idents(seu.obj) <- hash.ident
  ridge_plot <- RidgePlot(seu.obj, 
                          assay = HTO_AssayName, 
                          features = rownames(seu.obj[[HTO_AssayName]]),
                          ncol = ncol,
                          cols = cols) &
    theme(title = element_text(size = 7), 
          text = element_text(size = 7), 
          axis.text = element_text(size = 7))
  # FeatureScatter plots
  p <- list()
  for (j in 2:length(rownames(seu.obj))){
    p[[j]] <- FeatureScatter(seu.obj, 
                             feature1 = rownames(seu.obj[[HTO_AssayName]])[1], 
                             feature2 = rownames(seu.obj[[HTO_AssayName]])[j],
                             cols = cols) &
      NoLegend() &
      theme(aspect.ratio=1,
            title = element_text(size = 7), 
            text = element_text(size = 7), 
            axis.text = element_text(size = 7))
  }
  ft_scatter <- ggarrange(plotlist = p[-1], ncol = ncol, nrow = nrow)
  
  vln <- VlnPlot(seu.obj,
                 features = paste0(HTO_AssayName,"_margin"), 
                 split.by = hash_class_column,
                 pt.size = 0.1, 
                 cols = scales::hue_pal()(3)) &
    theme(title = element_text(size = 7), 
          text = element_text(size = 7), 
          axis.text = element_text(size = 7))
  
  # save plots
  pdf(file.path(save.loc, "plots", paste0(sample_name,".pdf")), width = plot.width, height = plot.height)
  distr_plot %>% print
  ridge_plot %>% print
  ft_scatter %>% print
  vln %>% print
  dev.off()
  
}
