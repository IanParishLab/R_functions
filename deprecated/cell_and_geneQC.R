cell_and_geneQC <- function(sce.obj, libsize_nmads = 3, mito_nmads = 2, lowcount_mean_threshold = -8, genome, cellQC = TRUE, geneQC = TRUE){
  require(SingleCellExperiment)
  require(tidyverse)
  
  par(mfrow = c(1,2))
  
  if(!is.null(cellQC)){
    print("[MSG] addPerCellQCMetrics...")
    
    if(genome == "human"){
      sce.obj <- addPerCellQCMetrics(x = sce.obj,subsets = list(Mito = grep("MT-",rownames(sce.obj))))
    } else if (genome == "mouse"){
      sce.obj <- addPerCellQCMetrics(x = sce.obj,subsets = list(Mito = grep("mt-",rownames(sce.obj))))
    }
    
    # cell QC
    print(paste0("[MSG] cellQC with libsize_nmads = ",libsize_nmads,"..."))
    
    libsize_drop <- isOutlier(
      metric = sce.obj$total,
      nmads = libsize_nmads,
      type = "lower",
      log = TRUE)
    colData(sce.obj)$libsize_drop <- libsize_drop
    
    print(paste0("[MSG] cellQC with mito_nmads = ",mito_nmads,"..."))
    mito_drop <- isOutlier(
      metric = colData(sce.obj)$subsets_Mito_percent,
      nmads = mito_nmads,
      type = "higher")
    colData(sce.obj)$mito_drop <- mito_drop
    
    plot_df <- data.frame(logtotal=log(sce.obj$total),
                          libsize_drop=factor(libsize_drop),
                          mito_drop=factor(mito_drop),
                          logdetected=log(sce.obj$detected))
    
    if(length(table(plot_df$libsize_drop)) > 1){
      plot_df$libsize_drop<-relevel(plot_df$libsize_drop, "TRUE")
    } else {
      print(paste0("[MSG] libsize_drop returns all FALSE, no cells to be dropped..."))
    }
    if(length(table(plot_df$mito_drop)) > 1){
      plot_df$mito_drop<-relevel(plot_df$mito_drop, "TRUE")
    } else {
      print(paste0("[MSG] mito_drop returns all FALSE, no cells to be dropped..."))
    }
    
    #A histogram showing the distribution of flagged cells
    p <- plot_df %>%
      ggplot( aes(x=logtotal, fill=libsize_drop)) +
      geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=30) +
      scale_fill_manual(values=c( "#69b3a2","#404080")) +
      # theme_ipsum() +
      labs(fill="") +
      xlab("Cell level log total count") +
      ylab("Frequency") +
      ggtitle('The distribution of cell-level log total counts\n  flagging cells with low library size')
    print(p)
  }
  
  # gene QC
  if(!is.null(geneQC)){
    print("[MSG] addPerFeatureQCMetrics...")
    
    sce.obj <- addPerFeatureQCMetrics(x = sce.obj)
    #Remove genes with zero counts for each gene
    sce.obj <- subset(sce.obj, rowData(sce.obj)$mean > 0)
    
    # detect low abundant genes
    lowcount_drop <- log(rowData(sce.obj)$mean) < lowcount_mean_threshold
    
    # mean count of genes across all cells
    plot_df2 <- data.frame(mean_genecount=log(rowData(sce.obj)$mean), 
                           lowcount_drop=factor(lowcount_drop))
    q <- plot_df2 %>%
      arrange(desc(lowcount_drop)) %>%
      ggplot(aes(x=mean_genecount, fill=lowcount_drop)) +
      geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity',bins=50) +
      scale_fill_manual(values=c("#404080", "#69b3a2")) +
      labs(fill="") + xlab("Log mean count across all cells") + ylab("Frequency") +
      ggtitle('The distribution of log mean count across\n all cells, flagging those with a low mean count')
    
    print(q)
  }
  
  return(list(libsize_drop,mito_drop,lowcount_drop))
}