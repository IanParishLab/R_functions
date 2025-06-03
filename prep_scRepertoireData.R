# Prepare files
prep_scRepertoireData <- function(seu.obj, contigs, sample_names = NULL, 
                                  seu.obj.ident = "PatientID", hash.ident = "sampleSource",
                                  # cloneTypes = c(Single=1, Small=5, Medium=20, Large=50, Larger=100, Hyperexpanded = 1000),
                                  cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.05, Hyperexpanded = 1),
                                  Tcell.type = "T-AB", save.name = "renamed", save.loc){
  
  library(scRepertoire)
  library(Seurat)
  library(tidyverse)
  
  # # test parameters
  # seu.obj = scvi
  # contig.dir = contig.dir
  # sample_names = capture
  # seu.obj.ident = "PatientID"
  # hash.ident = "sampleSource"
  # cloneTypes = cloneTypes
  # Tcell.type = "T-AB"
  # save.name = save.name
  # save.loc = save.loc
  # # end test parameters
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  # Match barcodes between contig and Seurat data
  seu.obj <- SplitObject(seu.obj, split.by = seu.obj.ident)
  
  for(i in 1:length(seu.obj)){
    seu.obj[[i]]@project.name <- sample_names[i]
    
    # check if contig and Seurat data is the same
    if(isFALSE(names(contigs)[i] == seu.obj[[i]]@project.name)){
      stop("contig data and Seurat object don't match")
    } 
    
    print(paste0("[MSG] ", 
                 paste0(table(Cells(seu.obj[[i]]) %in% contigs[[i]]$barcode)["TRUE"],
                        "/",ncol(seu.obj[[i]])),
                 " common GEX barcodes in TCR data..."))
  }
  contigs_all <- plyr::rbind.fill(contigs)
  saveRDS(contigs_all, file.path(save.loc, paste0(save.name, ".processed.filtered_contig_annotations.rds")))
  
  # Create HTO contig list
  if(length(seu.obj) > 1){
    print("[MSG] Merging Seurat objects for createHTOContigList...")
    seu.obj <- merge(seu.obj[[1]], seu.obj[2:length(seu.obj)])
  }
  
  if(is.list(seu.obj) & length(seu.obj) == 1) {
    seu.obj <- seu.obj[[1]]
  }
  
  DefaultAssay(seu.obj) <- "RNA"
  hto_contig_list <- createHTOContigList(contigs_all, seu.obj, group.by = c(hash.ident, seu.obj.ident))
  
  ind <- which(lapply(hto_contig_list, nrow) == 0)
  if(length(ind > 0)){
    print(paste0("[MSG] Removing ", length(ind), " HTO contigs because nrow == 0..."))
    hto_contig_list <- hto_contig_list[-ind]
  }
  
  # Add TCR expression data to Seurat object with combineTCR & combineExpression
  combined <- combineTCR(hto_contig_list, 
                         samples = names(hto_contig_list)
  )
  seu.obj$contigCellBarcode <- Cells(seu.obj)
  new.hto.contig.barcodes <- paste0(seu.obj@meta.data[[hash.ident]],".",
                                    seu.obj@meta.data[[seu.obj.ident]],"_",
                                    seu.obj@meta.data$contigCellBarcode)
  
  seu.obj <- RenameCells(seu.obj, new.names = new.hto.contig.barcodes)
  
  combined_seu.obj <- combineExpression(combined, seu.obj, 
                                        cloneCall = "strict", 
                                        cloneSize = cloneTypes,
                                        group.by = "sample",
                                        proportion = TRUE,
                                        #filterNA = TRUE,
                                        addLabel = TRUE
  )
  print("`cloneSize`:")
  table(combined_seu.obj$cloneSize) %>% print
  # range(combined_seu.obj$Frequency[which(combined_seu.obj$cloneType == "Small (0.001 < X <= 0.005) ")])
  saveRDS(combined_seu.obj, file.path(save.loc, paste0(save.name,".combineExpression.",Tcell.type,".rds")))
  
  # extract important metadata to add to the seu.obj
  add_this_metadata <- combined_seu.obj@meta.data %>% select(contigCellBarcode:cloneSize)
  rownames(add_this_metadata) <- combined_seu.obj@meta.data$cellBarcode
  
  return(add_this_metadata)
}
