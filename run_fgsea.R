# geneListDF = "/researchers/nicole.saw/projects/Christina_Scheffler/R/allmarkers.pool_filt2_cd8.logfc0.rds"
# geneListDF = presto::wilcoxauc(pool_filt3_cd4)
# refGMT = "/researchers/nicole.saw/references/gmt/beavis_gene_sig.gmt"
# clusters = seq(from=0, to=7)
# method = "wilcoxauc"

run_fgsea <- function(geneListDF, refGMT, readGMT = TRUE, seed.use,
                      cluster_column = "cluster", logfc_column = "avg_log2FC", gene_column = "gene",
                      nPerms = 20000, return_rnk = FALSE, plot_fgsea_res = TRUE){
  
  require(fgsea)
  require(tidyverse)
  
  set.seed(seed.use)
  
  if(readGMT == TRUE) {
    refGMT <- gmtPathways(refGMT)
  } else {
    refGMT <- refGMT
  }
  
  res <- list()
  ranks<-list()
  clusters <- unique(geneListDF[cluster_column]) %>% deframe
  for (j in clusters){
    unique_clusters <- gtools::mixedsort(unique(geneListDF[cluster_column]))
    
    ind <- which(geneListDF[cluster_column]==j)
    ranks[[j]] <- geneListDF[[logfc_column]][ind]
    names(ranks[[j]]) <- geneListDF[[gene_column]][ind]
    
    ranks[[j]] <- geneListDF %>%
      dplyr::filter(!!as.name(cluster_column) == j) %>%
      arrange(desc(!!as.name(logfc_column))) %>%
      dplyr::select(!!as.name(gene_column),!!as.name(logfc_column)) %>%
      deframe
    
    ranks[[j]] <- ranks[[j]][which(!is.na(ranks[[j]]))]
    
    res[[j]] <- fgseaMultilevel(pathways = refGMT, stats = ranks[[j]], nPermSimple = nPerms, nproc = 1)
    res[[j]][[cluster_column]] <- j
  }
  res_1 <-  plyr::rbind.fill(res)
  
  # for(i in 1:nrow(res_1)){
  #   res_1$leadingEdge[i] <- paste0(unlist(res_1$leadingEdge[i]), collapse = ",")
  # }

  # print plot  
  if(plot_fgsea_res == TRUE){
    p <- ggplot(res_1, aes(x = !!as.name(cluster_column), y = pathway)) +
      geom_point(aes(color = NES, size = -log(pval)), alpha = 0.7) +
      scale_colour_gradientn(colours = c("blue", "blue","white","red","red", "red"), oob = scales::squish) +
      scale_size(range = c(1, 10)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      scale_x_discrete(limits = c("", unique_clusters))
    print(p)
  }
  
  if(isTRUE(return_rnk)){
    return(ranks)
  }else {
    return(res_1)
  }
}
