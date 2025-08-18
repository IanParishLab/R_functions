process_CITEseqCount <- function(seu, assay = "HTO", counts_path){
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
  rownames(umis) <- gsub('-..*',"",rownames(umis))
  colnames(umis) <- paste0(colnames(umis),"-1")
  
  htos <- Read10X(htos_path, gene.column = 1)
  htos <- htos[rownames(htos) != "unmapped",]
  rownames(htos) <-  gsub('-..*',"",rownames(htos))
  colnames(htos) <- paste0(colnames(htos),"-1")
  
  # Select cell barcodes detected by both RNA and HTO and are "Singlets"
  joint.bc <- Reduce(intersect, list(colnames(htos),colnames(umis),colnames(seu)))
  print(length(joint.bc))
  # Subset RNA and HTO counts by joint cell barcodes
  umis <- umis[rownames(htos), joint.bc]
  htos <- as.matrix(htos[, joint.bc])
  
  # Add HTO data as a new assay independent from RNA
  seu[[assay]] <- CreateAssayObject(counts = htos)
  DefaultAssay(seu) <- assay
  seu <- NormalizeData(seu, assay = assay, normalization.method = "CLR", margin = 1)
  seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- HTODemux(seu, assay = assay, positive.quantile = 0.99)
  
  return(seu) # returns a demux'd CITEseqCount assay
}
