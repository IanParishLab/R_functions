run_topGO <- function(geneListDF, gene_column = "gene", foldchange_cutoff_column = "avg_log2FC", p_val_adj_column = "p_val_adj", foldchange_cutoff = NULL, p_val_adj_cutoff = 0.05,  mapping, topNodes) {
  
  # geneListDF = markers
  # gene_column = "gene"
  # p_val_adj_column = "p_val_adj"
  # mapping="org.Mm.eg.db"
  # topNodes=20
  
  require(topGO)
  require(tidyverse)
  require(stringr)
  
  assayed.genes <- geneListDF[[gene_column]]
  de.genes <- geneListDF[[gene_column]][ which(geneListDF[[p_val_adj_column]] < p_val_adj_cutoff & abs(geneListDF[[foldchange_cutoff_column]]) < foldchange_cutoff) ]
  genes <-factor(as.integer(geneListDF[[p_val_adj_column]] < 0.05))
  names(genes) <- sub("[|].*", "", geneListDF[[gene_column]])
  genes <- genes[!is.na(genes)]
  names(genes)
  
  onts = c( "MF", "BP", "CC" )
  tab = as.list(onts)
  names(tab) = onts
  
  for(j in 1:3){
    tgd <- new("topGOdata", ontology=onts[j], allGenes = genes, geneSelectionFun = function(x)(x == 1), nodeSize=20, annot=annFUN.org, mapping=mapping, ID = "symbol" )
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    tab[[j]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, Fisher.classic = resultTopGO.classic, orderBy = "Fisher.classic" , topNodes = topNodes)
    genes_GO <- genesInTerm(tgd)
    genes_annot = lapply(genes_GO,function(x) x[x %in% de.genes] )
    check  <- apply(tab[[j]],1, function(x) genes_annot[[x[1]]])
    tab[[j]]$genes <- unlist(lapply(check, function(x) paste0(x, collapse=',')), use.names=FALSE)
  }
  topGOResults <- plyr::rbind.fill(tab)
  topGOResults <- topGOResults[order(as.numeric(topGOResults$Fisher.elim)),]
  
  return(topGOResults)
}
