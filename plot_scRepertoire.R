plot_scRepertoire <- function(so, 
                              seu.obj.ident = "PatientID", 
                              hash.ident = "sampleSource", 
                              clonotype.ident = "cloneSize", 
                              cluster.ident = "seuratClusters",
                              clonotype.ident.levels, 
                              cloneTypes = c(Rare = 0.001, Small = 0.005, Medium = 0.01, Large = 0.05, Expanded = 0.1, Hyperexpanded = 1),
                              reduction.name = "umap.scvi", 
                              hash.superset.ident = "sampleSourceSuperset",
                              circle_filter_clusters = NULL, 
                              save.loc = "./scRepertoire", 
                              save.name = "CD8"){
  
  suppressMessages(library(ggraph))
  suppressMessages(library(ggalluvial))
  suppressMessages(library(ggpubr))
  suppressMessages(library(scRepertoire))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(patchwork))
  suppressMessages(library(cowplot))
  suppressMessages(library(circlize))
  
  # Split seurat object for per patient clonotype analysis 
  ## Basic quantification of clonotypes
  # Visualise distribution of clonotype bins on the UMAP:
  clonotype.ident.levels <- levels(so$cloneSize)

  Idents(so) <- cluster.ident
  # split object
  split_so <- SplitObject(so, split.by = seu.obj.ident)
  colorblind_vector <- setNames(
    c(colorRampPalette(viridisLite::viridis(option="viridis",n = 10), bias=0.7)(length(clonotype.ident.levels)-1),'grey'),
    clonotype.ident.levels)
  
  colorblind_vector <- colorblind_vector[clonotype.ident.levels]
  colorblind_vector <- colorblind_vector[which(!is.na(colorblind_vector[clonotype.ident.levels]))]
  
  so@meta.data[[clonotype.ident]] <- factor(so@meta.data[[clonotype.ident]], levels = clonotype.ident.levels)
  
  Idents(so) <- clonotype.ident
  q<-list()
  for(i in 1:length(colorblind_vector)){
    q[[i]] <- DimPlot(so, split.by = seu.obj.ident,
                      reduction = reduction.name, group.by = clonotype.ident, order =TRUE,
                      cells.highlight = list(WhichCells(so, idents = clonotype.ident.levels[i])),
                      cols.highlight = list(colorblind_vector[i]),
                      sizes.highlight = 0.3,
                      cols='grey') + 
      ggtitle(clonotype.ident.levels[i]) +
      NoLegend() + 
      NoAxes() +
      theme(aspect.ratio = 1, text = element_text(size = 10), plot.title = element_text(hjust = 0.5, size = 10))
  }
  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonetypeLocationComposition.pdf")), width = 8, height = 8)
  ggarrange(plotlist = q, nrow = 2, ncol = 3) %>% print#length(colorblind_vector)) %>% print
  dittoSeq::dittoBarPlot(so, split.by = seu.obj.ident, group.by = cluster.ident, var = clonotype.ident, xlab = "Cluster", split.ncol = 2, split.nrow = 3,
                         color.panel = colorblind_vector, legend.show = TRUE, scale = "percent", retain.factor.levels = TRUE) %>% print
  dev.off()
  
  # Visualize clonal frequency onto UMAPs to show clonal expansion:
  Idents(so) <- cluster.ident
  pdf(file.path(save.loc, paste0(save.name,".clonalOverlay.pdf")), width = 10, height = 8)
  lapply(cloneTypes[-length(cloneTypes)], function(cutoff){
    clonalOverlay(so,
                  reduction = reduction.name,
                  freq.cutpoint = cutoff,
                  bins = 20,
                  facet = seu.obj.ident) +
      ggtitle(cutoff) +
      scale_color_manual(values = cluster_col) +
      theme(aspect.ratio = 1) +
      guides(color = guide_legend(override.aes = list(size = 2)))
  }) %>% print
  dev.off()

  # Quantify numbers of unique clonotypes, here it is scaled to the total number of clonotypes:
  ## use split seurat object
  # A chain
  p <- lapply(split_so, function(seu){
    quantContig(seu, cloneCall="strict", scale = TRUE,
                split.by = hash.ident,chain = "TRA") +
      ylim(0,100) +
      labs(subtitle = unique(seu[[seu.obj.ident]]), x = "Sample source") +
      theme(axis.text.x = element_text(hjust = 1, angle = 45),
            legend.position = "none")
  })
  a_chain <- ggarrange(plotlist = p, nrow=3,ncol=2)
  title <- ggdraw() + draw_label("A chain", fontface = 'bold')
  a_chain <- plot_grid(title, a_chain, ncol = 1, rel_heights = c(0.1, 1))

  # B chain
  p <- lapply(split_so, function(seu){
    quantContig(seu, cloneCall="strict", scale = TRUE,
                split.by = hash.ident,chain = "TRB") +
      ylim(0,100) +
      labs(subtitle = unique(seu[[seu.obj.ident]]), x = "Sample source") +
      theme(axis.text.x = element_text(hjust = 1, angle = 45),
            legend.position = "none")
  })
  b_chain <- ggarrange(plotlist = p, nrow=3,ncol=2)
  title <- ggdraw() + draw_label("B chain", fontface = 'bold')
  b_chain <- plot_grid(title, b_chain, ncol = 1, rel_heights = c(0.1, 1))

  # combined A & B chain plots of unique clonotypes
  ab_chain <- a_chain + b_chain

  # save plot
  pdf(file.path(save.loc, paste0(save.name,".uniqueClonotypes.pdf")), width = 10, height = 8)
  ab_chain %>% print
  dev.off()

  # Visualize the number of clonotypes at specific frequencies by sample source:
  p <- lapply(split_so, function(seu){
    abundanceContig(seu, split.by = hash.ident,
                    cloneCall = "strict",
                    scale = FALSE,
                    chain = "TRA") +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  contig_freq_a <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = TRUE)
  title <- ggdraw() + draw_label("A chain", fontface = 'bold')
  contig_freq_a <- plot_grid(title, contig_freq_a, ncol = 1, rel_heights = c(0.1, 1))

  p <- lapply(split_so, function(seu){
    abundanceContig(seu, split.by = hash.ident,
                    cloneCall = "strict",
                    scale = FALSE,
                    chain = "TRB") +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  contig_freq_b <- ggarrange(plotlist = p, nrow=3, ncol=2, common.legend = TRUE)
  title <- ggdraw() + draw_label("B chain", fontface = 'bold')
  contig_freq_b <- plot_grid(title, contig_freq_b, ncol = 1, rel_heights = c(0.1, 1))

  # combined A & B chain plots of clonetype abundance
  contig_freq_ab <- contig_freq_a + contig_freq_b

  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonetypeAbundance.pdf")), width = 8, height = 6)
  contig_freq_ab %>% print
  dev.off()

  ## More advanced quantification of clonotypes
  # Count of cells in each cluster, assigned into specific frequency ranges:
  p <- lapply(split_so, function(seu){
    occupiedscRepertoire(seu, x.axis = hash.ident, proportion = TRUE, facet.by = seu.obj.ident) +
      scale_color_manual(values = colorblind_vector) +
      theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  # plots of clonotype freq
  clonotype_freq <- ggarrange(plotlist = p, nrow=3,ncol=2, legend.grob = get_legend(p[[length(p)]]), legend = "right")
  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonetypeFreqCount.pdf")), width = 8, height = 10)
  clonotype_freq %>% print
  dev.off()


  #Clonal space homeostasis asks what % of total immune receptor sequencing run is filled by clones in distinct proportions:
  p <- lapply(split_so, function(seu){
    clonalHomeostasis(seu, split.by = hash.ident,
                      cloneCall = "strict",
                      cloneTypes = cloneTypes) +
      scale_color_manual(values = colorblind_vector) +
      theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })

  # plots of clonal space occupied in immune receptor sequencing
  clonotype_freq <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = TRUE)

  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonetypeImmuneFreq.pdf")), width = 10, height = 8)
  clonotype_freq %>% print
  dev.off()
  
  # `clonalProportion()` ranks the clones by total number and place them into bins, instead of looking at the relative proportion of the clone
  p <- lapply(split_so, function(seu){
    clonalProportion(seu, split.by = hash.ident,
                     cloneCall = "strict",
                     split = c(1,5,10,20,50,100,1000,10000,30000)) +
      theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })

  # Relative proportional space occupied by specific clonotypes
  clonotype_prop <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = TRUE)

  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonetypeRelativeProportion.pdf"), width = 6, height = 4))
  clonotype_prop
  dev.off()
  
  # Visualise two-sample, direct comparison of clonotypes using `scatterClonotype()`
  # * The clonotypes are categorized by counts into singlets or multiplets, and are either exclusive or shared between the selected samples
  p <- lapply(split_so, function(seu){
    scatterClonotype(seu, split.by = hash.superset.ident,
                     cloneCall = "strict",
                     x.axis = "PBMC",
                     y.axis = "LN",
                     dot.size = "total",
                     graph = "proportion") +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })

  # Relative proportional space occupied by specific clonotypes
  two_sample_scatter <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = TRUE)
  # save plot
  pdf(file.path(save.loc, paste0(save.name,".twoSampleScatter.pdf"), width = 6, height = 4))
  two_sample_scatter
  dev.off()
  
  # Clonotypes between samples and changes in dynamics, visualized in an alluvial plot
  p <- lapply(split_so, function(seu){
    compareClonotypes(seu,
                      split.by = hash.superset.ident,
                      numbers = 3,
                      samples = c("PBMC", "LN"),
                      cloneCall = "strict",
                      chain = "TRA",
                      graph = "alluvial") +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  a_alluv <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = FALSE)
  a_alluv

  p <- lapply(split_so, function(seu){
    compareClonotypes(seu,
                      split.by = hash.superset.ident,
                      numbers = 3,
                      samples = c("PBMC", "LN"),
                      cloneCall = "strict",
                      chain = "TRB",
                      graph = "alluvial") +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  b_alluv <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = FALSE)
  b_alluv
  
  # Visualise clonotypes across multiple categories, also with an alluvial plot
  p <- lapply(split_so, function(seu){
    alluvialClonotypes(so,
                       cloneCall = "strict",
                       chain = "TRA",
                       y.axes = c(seu.obj.ident, cluster.ident, hash.ident),
                       color = cluster.ident) +
      scale_fill_manual(values = cluster_col) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })

  a_alluv_multi <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = TRUE)

  p <- lapply(split_so, function(seu){
    alluvialClonotypes(so,
                       cloneCall = "strict",
                       chain = "TRB",
                       y.axes = c(seu.obj.ident, cluster.ident, hash.ident),
                       color = cluster.ident) +
      scale_fill_manual(values = cluster_col) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  b_alluv_multi <- ggarrange(plotlist = p, nrow=3,ncol=2, common.legend = TRUE)

  alluv_multi <- a_alluv_multi + b_alluv_multi
  
  # Length distribution of the CDR3 sequences - should reveal a multimodal curve, unless single chain selected.
  # the multimodal curve is a product of using the NA for the unreturned chain sequence and multiple chains within a single barcode.
  p <- lapply(split_so, function(seu){
    p1 <- lengthContig(seu,
                       split.by = hash.superset.ident ,
                       cloneCall = "nt",
                       chain = "both",
                       scale = TRUE) +
      labs(subtitle = "Both")

    p2 <- lengthContig(seu,
                       split.by = hash.superset.ident ,
                       cloneCall = "nt",
                       chain = "TRA",
                       scale = TRUE)  +
      labs(subtitle = "TRA chain")

    p3 <- lengthContig(seu,
                       split.by = hash.superset.ident ,
                       cloneCall = "nt",
                       chain = "TRB",
                       scale = TRUE)  +
      labs(subtitle = "TRB chain")

    multi_p <- ggarrange(plotlist = list(p1,p2,p3), ncol=3, legend = "none")
    title <- ggdraw() + draw_label(unique(seu[[seu.obj.ident]]))
    multi_p <- plot_grid(title, multi_p, ncol = 1, rel_heights = c(0.1, 1))
    return(multi_p)
  })
  length_dist_plot <- ggarrange(plotlist = p, nrow = 5, common.legend = TRUE)

  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonetypeImmuneFreq.pdf")), width = 8, height = 10)
  length_dist_plot %>% print
  dev.off()
  
  # Visualizing the distribution of any VDJ and C gene of the TCR, showing relative usage of TCR genes
  p <- lapply(split_so, function(seu){
    q <- lapply(c("V","J","C"), function(gn){
      vizGenes(seu,
               split.by = hash.superset.ident,
               gene = gn,
               chain = "TRA",
               plot = "heatmap",
               order = "variance",
               scale = TRUE) +
        labs(subtitle = paste0(gn, " genes"))
    })
    q <- ggarrange(plotlist = q, ncol = 3)
    title <- ggdraw() + draw_label(unique(seu[[seu.obj.ident]]))
    multi_p <- plot_grid(title, q, ncol = 1, rel_heights = c(0.1, 1))
  })

  # plot of relative TCR gene usage
  vizgene_plots <- ggarrange(plotlist = p, ncol = 1)
  vizgene_plots

  # # single chain
  # ind <- grep("LN",names(combined_TAB),value=TRUE)
  # vizGenes(combined_TAB[ind],
  #          gene = "V",
  #          chain = "TRB",
  #          y.axis = "J",
  #          plot = "heatmap",
  #          scale = TRUE,
  #          order = "gene")
  
  ## Network analysis
  # Generates a network based on clonal proportions of sample source, visualise network interaction of clonotypes shared between clusters
  # with Identity filter
  p <- lapply(split_so, function(seu){
    clonalNetwork(seu,
                  reduction = reduction.name,
                  identity = hash.ident,
                  filter.clones = NULL,
                  # filter.identity = c("PBMC", "LN"),
                  cloneCall = "strict") +
      # scale_colour_manual(values = cluster_col) +
      theme(aspect.ratio = 1) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  network_plot <- ggarrange(plotlist = p, nrow = 3, ncol = 2)
  # save plot
  pdf(file.path(save.loc, paste0(save.name,".clonalNetwork.pdf")), width = 10, height = 10)
  network_plot %>% print
  dev.off()
  
  # Using `circlize` to visualize the interconnection of clusters
  # Each chord represents the number of clone unique and shared across the multiple Clusters
  pdf(file.path(save.loc, paste0(save.name,".circlePlot.pdf")), width = 8, height = 6)
  par(mfrow = c(2,3))
  lapply(split_so, function(seu){
    circles <- getCirclize(seu, group.by = cluster.ident)
    # Graphing the chord diagram
    circ <- chordDiagram(circles,
                         self.link = 1, 
                         grid.col = cluster_col)
    text(0,0, unique(seu[[seu.obj.ident]]), cex=2)
    return(circ)
  })
  
  if(!is.null(circle_filter_clusters)){
    par(mfrow = c(2,3))
    lapply(split_so, function(seu){
      circles <- getCirclize(seu, group.by = cluster.ident)
      circles <- circles %>% filter(from %in% circle_filter_clusters | to %in% circle_filter_clusters) %>% filter(value != 0)
      # circles <- getCirclize(seu, group.by = hash.ident)
      # circles <- circles[which(circles$from %in% c("idLN","udLN","undLN") | circles$to %in% c("idLN","udLN","undLN")),]
      # Graphing the chord diagram
      circ <- chordDiagram(circles,
                           self.link = 1, 
                           grid.col = cluster_col)
      text(0,0, unique(seu[[seu.obj.ident]]), cex=2)
      return(circ)
    })
  }
  dev.off()
  
  ## Clustering Clonotypes
  # Directly cluster TCR cells within the single-cell object using the same approach as above. Uses edit distances of either the nt or aa of the CDR3 to cluster similar TCRs together.
  # * Edit distance, or Levenshtein distance measures the similarity of the V- and J-region used to assemble the TCR beta-chain.
  so <- clusterTCR(so,
                   chain = "TRA", # choose one of TRA, TRB, TRG, TRD
                   group.by = cluster.ident,
                   sequence = "aa",
                   threshold = 0.8)

  pdf(file.path(save.loc, paste0(save.name,".clusterTCR.pdf"), width = 12, height = 8))
  DimPlot(so, group.by = "TRA_cluster") +
    scale_color_manual(values = colorblind_vector(length(unique(so@meta.data[,"TRA_cluster"])))) +
    theme(aspect.ratio = 1)
  dev.off()
  
  ## Diversity Analysis
  # Calculate traditional measures of diversity. The former two metrics are generally used to estimate baseline diversity, Chao/ACE indices are used to estimate the richness of the samples.
  p <- lapply(split_so, function(seu){
    clonalDiversity(seu,
                    cloneCall = "strict",
                    chain = "TRA",
                    split.by = hash.ident,
                    n.boots = 100) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  a_div <- ggarrange(plotlist = p, nrow=5,ncol=1, common.legend = TRUE)
  title <- ggdraw() + draw_label("A chain", fontface = 'bold')
  a_div <- plot_grid(title, a_div, ncol = 1, rel_heights = c(0.1, 1))
  
  p <- lapply(split_so, function(seu){
    clonalDiversity(seu,
                    cloneCall = "strict",
                    chain = "TRB",
                    split.by = hash.ident,
                    n.boots = 100) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
  })
  b_div <- ggarrange(plotlist = p, nrow=5,ncol=1, common.legend = TRUE)
  title <- ggdraw() + draw_label("B chain", fontface = 'bold')
  b_div <- plot_grid(title, b_div, ncol = 1, rel_heights = c(0.1, 1))
  
  ab_div <- a_div + b_div
  # save plots
  pdf(file.path(save.loc, paste0(save.name,".clonalDiversity.pdf")), width = 12, height = 10)
  ab_div %>% print
  dev.off()
  
  # ## omit
  # # Split by specific samples
  # ind <- c(grep("PBMC",names(combined_TAB),value=TRUE),
  #          grep("LN",names(combined_TAB),value=TRUE))
  # clonalDiversity(combined_TAB[ind],
  #                 cloneCall = "gene",
  #                 x.axis = "sample",
  #                 n.boots = 100) +
  #   theme(axis.text.x = element_text(hjust = 1, angle = 45))
  
  # [StartracDiversity](https://www.nature.com/articles/s41586-018-0694-x) - methods for looking at clonotypes by cellular origins and cluster identification.
  # * STARTRAC-expa, STARTRAC-migr and STARTRAC-tran, are designed to measure the degree of clonal expansion, tissue migration, and state transition of TRUE cell clusters upon TCR tracking, respectively
  # * Calculates diversity on a clonotype level instead of a sample level (traditional)
  Idents(so) <- hash.ident
  s_div <- StartracDiversity(so,
                             type = hash.ident,
                             sample = seu.obj.ident,
                             by = "overall") +
    theme(axis.text.x = element_text(hjust = 1, angle = 45)) +
    scale_color_manual(values = cluster_col)
  # save plots
  pdf(file.path(save.loc, paste0(save.name,".StartracDiversity.pdf")), width = 6, height = 4)
  s_div %>% print
  dev.off()
  
  ## Overlap analysis
  # We can also visualise measures of similarity between the samples, here we use the overlap coefficient - which shows the overlap of clonotypes scaled to the length of unique clonotypes:
  pal = colorRamp2(seq(0,0.4,0.1), c("#000004","#270C4C","#872168","#EC6726","#FCFFA4"))
  
  #All patients
  p1 <- (clonalOverlap(so, split.by = hash.ident,
                       cloneCall = "strict",
                       method = "overlap") +
           theme(axis.text.x = element_text(hjust = 1, angle = 30)) +
           update_geom_defaults("text", list(size = 2)) +
           labs(subtitle = "All patients")) +
    (clonalOverlap(so, split.by = hash.superset.ident,
                   cloneCall = "strict",
                   method = "overlap") +
       theme(axis.text.x = element_text(hjust = 1, angle = 30)) +
       update_geom_defaults("text", list(size = 2)) +
       labs(subtitle = "All patients"))
  
  # split by patient, sample source
  p <- lapply(split_so, function(seu){
    a <- clonalOverlap(seu, split.by = hash.ident,
                       cloneCall = "strict",
                       method = "overlap") +
      theme(axis.text.x = element_text(hjust = 1, angle = 30)) +
      update_geom_defaults("text", list(size = 2)) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
    return(a)
  })
  overlap_plot1 <- ggarrange(plotlist = p , nrow=2, ncol=3)
  
  # split by patient, clusters
  p <- lapply(split_so, function(seu){
    a <- clonalOverlap(seu, split.by = cluster.ident,
                       cloneCall = "strict",
                       method = "overlap") +
      update_geom_defaults("text", list(size = 2)) +
      theme(axis.text.x = element_text(hjust = 1, angle = 30)) +
      labs(subtitle = unique(seu[[seu.obj.ident]]))
    return(a)
  })
  overlap_plot2 <- ggarrange(plotlist = p, nrow=2, ncol=3)
  
  pdf(file.path(save.loc, paste0(save.name,".clonalOverlapPlot_sampleSource_clusters.pdf")), width = 12, height = 5)
  overlap_plot1 %>% print
  overlap_plot2 %>% print
  dev.off()
}
