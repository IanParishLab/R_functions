process_CITEseqCount <- function(seu, assay = "HTO", counts_path, hash_names = NULL, exclude_samples = NULL){
  # counts_path = here("../CITE-seq-Count/TotalSeqA/")
  # assay = "HTO"
  
  # counts_path = here("../CITE-seq-Count/MULTIseq/")
  # assay = "LMO"
  # 
  # seu <- multi

  # set paths
  umi_path <- file.path(counts_path,"umi_count")
  htos_path <- file.path(counts_path,"read_count")
  
  umis <- Read10X(umi_path, gene.column = 1)
  umis <- umis[rownames(umis) != "unmapped",]
  colnames(umis) <- paste0(colnames(umis),"-1")
  
  htos <- Read10X(htos_path, gene.column = 1)
  htos <- htos[rownames(htos) != "unmapped",]
  colnames(htos) <- paste0(colnames(htos),"-1")
  
  if(is.null(hash_names)){
    rownames(umis) <- gsub('-..*',"",rownames(umis))
    rownames(htos) <-  gsub('-..*',"",rownames(htos))
  } else {
    rownames(umis) <- hash_names
    rownames(htos) <-  hash_names
  }
  
  # Select cell barcodes detected by both RNA and HTO and are "Singlets"
  joint.bc <- Reduce(intersect, list(colnames(htos),colnames(umis),colnames(seu)))
  print(length(joint.bc))
  # Subset RNA and HTO counts by joint cell barcodes
  umis <- umis[rownames(htos), joint.bc]
  htos <- as.matrix(htos[, joint.bc])
  seu <- subset(seu, cells = joint.bc)
  
  # Add HTO data as a new assay independent from RNA
  if(!is.null(exclude_samples)){
    message("[MSG] Removing these samples from HTO assay: ", paste(exclude_samples, collapse = ","))
    ind <- which(rownames(htos) %in% exclude_samples)
    seu[[assay]] <- CreateAssayObject(counts = htos[-ind,]) 
    
    # print hashes rowSums
    message("[MSG] Hashes remaining & their `rowSums()`...")
    print(rowSums(seu[[assay]]))
    
  } else {
    seu[[assay]] <- CreateAssayObject(counts = htos)
  }
  
  DefaultAssay(seu) <- assay
  seu <- NormalizeData(seu, assay = assay, normalization.method = "CLR", margin = 1)
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- HTODemux(seu, assay = assay, positive.quantile = 0.99)
  
  return(seu) # returns a demux'd CITEseqCount assay
}
