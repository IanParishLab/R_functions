run_emptydrops <- function(counts, limit = 100, add.limit=100, method, plot.path, remove.genes=NULL) {
  require(DropletUtils)
  require(scater)
  
  # counts <- data[[1]]$`Gene Expression`
  # plot.path <- paste0("emptydrops/",capture[i],".pdf")
  # method <- "emptyDropsCellRanger"
  # limit = 100
  
  br.out <- barcodeRanks(counts, lower = limit)
  if (metadata(br.out)$knee < 8^3){
    print(paste0("knee = ",metadata(br.out)$knee, ", removing genes if any..."))
    counts <- counts[-which(rownames(counts) %in% c(remove.genes)),]
  }
  
  pdf(file = plot.path, width = 15, height = 5)
  par(mfrow=c(1,3))
  
  # barcodeRank plot
  if (metadata(br.out)$inflection < 150) {
    limit <- limit + add.limit
  }
  br.out <- barcodeRanks(counts, lower = limit)
  
  # print final barcode rank result
  print(metadata(br.out))
  
  uniq <- !duplicated(br.out$rank)
  plot(br.out$rank[uniq], br.out$total[uniq], log="xy",
       xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
  abline(h=metadata(br.out)$inflection, col="darkgreen", lty=2)
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  legend("bottomleft", legend=c(paste0("Inflection=",metadata(br.out)$inflection), paste0("Knee=",metadata(br.out)$knee)),
         col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
  
  # p-value distribution histograms
  # also choosing a emptyDrops method
  if (method == "emptyDrops") {
    print("Using emptyDrops...")
    e.out <- emptyDrops(counts, test.ambient = TRUE, niters = 15000)
    hist(e.out$PValue[e.out$Total <= limit & e.out$Total > 0],
         xlab="P-value", col="grey80", main = "emptyDrops")
    
  } else if (method == "emptyDropsCellRanger"){
    print("Using emptyDropsCellRanger...")
    e.out <- emptyDropsCellRanger(counts, niters = 15000) 
    hist(e.out$PValue[e.out$Total > 0], main = "emptyDropsCellRanger")
    
  } else if (metadata(br.out)$knee < 8^3){
    print("Using emptyDrops with retain = Inf...")
    e.out <- emptyDrops(counts, test.ambient = TRUE, niters = 15000, retain = Inf)
    hist(e.out$PValue[which(e.out$Total <= limit & e.out$Total > 0)],
         xlab="P-value", col="grey80", main = "emptyDrops with retain = Inf")
    
  } else {
    stop("what is going on?")
  }
  
  # stop if no non-empty droplets
  if (table(e.out$FDR < 0.001)[["TRUE"]] == 0){
    stop("no non-empty droplets")
  }
  
  print(summary(e.out$FDR < 0.001))
  
  # MA plot
  is.cell <- e.out$FDR <= 0.001
  filtered.counts <- counts[,which(is.cell)]
  cells <- calculateAverage(filtered.counts)
  ambient.cells <- which(is.na(e.out$FDR))
  ambient <- rowSums(counts[, ambient.cells])
  
  colour <- rep('black', length(cells))
  colour[names(cells) %in% remove.genes] <- "dodgerblue"
  # colour[names(cells) %in% tcr_bcr_genes] <- "red"
  
  edgeR::maPlot(ambient, cells, normalize = TRUE, col = colour)
  
  print(table(is.cell, useNA='ifany'))
  dev.off()
  
  return(e.out)
}


# tcr_bcr_genes <-c(grep("^TR[ABDG][VJC]", rownames(counts),value=T), grep("^IG[HKL]V|^IG[HKL]J|^IG[KL]C|^IGH[ADEGM]", rownames(counts),value=T))
