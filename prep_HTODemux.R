prep_HTODemux <- function(seu.obj, HTO.counts.path, sample_name, HTO_AssayName = "HTO", RNA_AssayName = "RNA", exclude_samples = NULL, demux = FALSE){
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
  })
  
  HTO.counts <- Read10X(HTO.counts.path)$`Antibody Capture`
  rownames(HTO.counts) <- gsub(" ","-",rownames(HTO.counts))
  colnames(HTO.counts) <- gsub(paste0(sample_name,"_"),"", colnames(HTO.counts))
  
  if(isTRUE(demux)) {
    
    if(!is.null(exclude_samples)){
      message("[MSG] Removing these samples from HTO assay: ", paste(exclude_samples, collapse = ","))
      ind <- which(rownames(HTO.counts) %in% exclude_samples)
      seu.obj[[HTO_AssayName]] <- CreateAssayObject(counts = HTO.counts[-ind, Cells(seu.obj)]) 
      
      # print hashes rowSums
      message("[MSG] Hashes remaining & their `rowSums()`...")
      print(rowSums(seu.obj[[HTO_AssayName]]))
      
    } else {
      seu.obj[[HTO_AssayName]] <- CreateAssayObject(counts = HTO.counts[, Cells(seu.obj)])
    }
    
    DefaultAssay(seu.obj) <- HTO_AssayName
    seu.obj <- NormalizeData(seu.obj, assay = HTO_AssayName, normalization.method = 'CLR')
    seu.obj <- ScaleData(seu.obj, assay = HTO_AssayName)
    seu.obj <- HTODemux(seu.obj, assay = HTO_AssayName, positive.quantile = 0.99)
    
  }
  
  return(seu.obj) 
}
