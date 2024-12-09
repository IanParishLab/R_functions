prep_HTODemux <- function(seu.obj, HTO.counts.path, sampleName, HTO_AssayName = "HTO", RNA_AssayName = "RNA", exclude_samples = NULL, demux = FALSE){
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
  })
  
  HTO.counts <- Read10X(HTO.counts.path)$`Antibody Capture`
  # HTO.counts <- HTO.counts$decontaminated_counts
  rownames(HTO.counts) <- gsub(" ","-",rownames(HTO.counts))
  colnames(HTO.counts) <- gsub(paste0(sampleName,"_"),"", colnames(HTO.counts))
  
  if(!is.null(exclude_samples)){
    print(paste0("[MSG] Removing these samples from HTO assay: ", paste(exclude_samples, collapse = ",")))
    ind <- which(rownames(HTO.counts) %in% exclude_samples)
    seu.obj[[HTO_AssayName]] <- CreateAssayObject(counts = HTO.counts[-ind,colnames(seu.obj)]) 
  } else {
    seu.obj[[HTO_AssayName]] <- CreateAssayObject(counts = HTO.counts[,colnames(seu.obj)])
  }
  
  # print hashes remaining
  print(paste0("[MSG] Hashes remaining and their `rowSums()`..."))
  print(rowSums(seu.obj[[HTO_AssayName]]))
  
  if(isTRUE(demux)) {
    
    DefaultAssay(seu.obj) <- HTO_AssayName
    seu.obj <- NormalizeData(seu.obj, assay = HTO_AssayName, normalization.method = 'CLR')
    seu.obj <- ScaleData(seu.obj, assay = HTO_AssayName)
    seu.obj <- HTODemux(seu.obj, assay = HTO_AssayName, positive.quantile = 0.99)
    
  }
  
  return(seu.obj)
}
