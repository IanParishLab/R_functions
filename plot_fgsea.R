replot_fgsea <- function(fgsea_res_df, ranks_stat, pathways){
  
  # fgsea_res_df=human
  # ranks_stat=ranks_stat
  # pathways= IA_vs_Ctrl
  # plot_pathway='IA_vs_Ctrl_Up'
  
  gsea.rnk<-as.data.frame(ranks_stat)
  gsea.rnk$hgnc.symbol<-rownames(gsea.rnk)
  rownames(gsea.rnk)<-NULL
  colnames(gsea.rnk) <- c("metric","hgnc.symbol")
  metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
  gsea.rnk<-gsea.rnk[order(gsea.rnk$metric),]
  gsea.metric<-"None"
  
  # # Get template name
  gsea.template <- ""
  gsea.gene.set <- plot_pathway
  # Get enrichment score
  gsea.enrichment.score<- fgsea_res_df$NES[which(fgsea_res_df$pathway == plot_pathway)]
  gsea.enrichment.score<-round(gsea.enrichment.score, 4)
  # # Get nominal p-value
  gsea.p.value<-fgsea_res_df$padj[which(fgsea_res_df$pathway == plot_pathway)]
  gsea.p.value<-as.numeric(gsea.p.value)
  gsea.p.value<-format.pval(gsea.p.value,eps=gsea.p.value,digits = 4)
  # # Get hit indices
  rnk <- rank(-ranks_stat)
  ord <- order(rnk)
  gseaParam<-1
  statsAdj <- ranks_stat[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathways[[plot_pathway]], names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  
  # Get gene set name
  gsea.normalized.enrichment.score<-fgsea_res_df$NES[which(fgsea_res_df$pathway == gsea.gene.set)]
  gsea.normalized.enrichment.score<-round(gsea.normalized.enrichment.score, 4)
  print(gsea.normalized.enrichment.score)
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  diff <- (max(tops) - min(bottoms)) / 8
  gsea.hit.indices <- toPlot$x
  
  # # Get ES profile
  if (gsea.normalized.enrichment.score<0) {
    toPlot$y[sapply(toPlot$y, is.numeric)] <- toPlot$y[sapply(toPlot$y, is.numeric)] * -1
    gsea.es.profile <- toPlot$y
    enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
  }else{
    gsea.es.profile <- toPlot$y
    enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
  }
  
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3)), heights = c(1.7, 0.5, 0.2))
  # layout.show(gsea.layout)
  
  # Create plots
  
  # Enrichment curve
  par(mar = c(0, 5, 3, 2))
  if (gsea.normalized.enrichment.score<0) {
    plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
         c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 5, xaxt = "n",
         xaxs = "i", xlab = "", ylab = "", 
         ylim = enrichment.score.range, 
         cex.axis = 3, cex.main = 1, cex.lab= 3, 
         main = list(paste0("GENESET: ",plot_pathway," NES: ",gsea.normalized.enrichment.score), font = 1),
         panel.first = {
           abline(h = seq(round(enrichment.score.range[1], digits = 1),
                          enrichment.score.range[2], 0.1),
                  col = "gray95", lty = 2)
           abline(h = 0, col = "gray50", lty = 2, lwd = 5)
         })
  }else{
    plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
         c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 5, xaxt = "n",
         xaxs = "i", xlab = "", ylab = "", 
         ylim = enrichment.score.range, 
         cex.axis = 3, cex.main = 1, cex.lab= 3, 
         main = list(paste0("GENESET: ",plot_pathway," NES: ",gsea.normalized.enrichment.score), font = 1),
         panel.first = {
           abline(h = seq(round(enrichment.score.range[1], digits = 1),
                          enrichment.score.range[2], 0.1),
                  col = "gray95", lty = 2)
           abline(h = 0, col = "gray50", lty = 2, lwd = 5)
         })
  }
  
  # Hit indices
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))
  abline(v = gsea.hit.indices, lwd = 0.75)
  par(mar = c(0, 5, 0, 2))
  rank.colors <- gsea.rnk$metric - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea.rnk$metric)))
  box()
  text(length(gsea.rnk$metric) / 2, 0.7,labels = gsea.template, cex = 2.5)
  text(length(gsea.rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA),cex=2.5)
  text(length(gsea.rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA),cex=2.5)
  box()
  
  # Reset to default
  par(def.par)
}
