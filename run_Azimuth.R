run_Azimuth <- function(obj, reference, save.loc, avail.refs = FALSE, ...){
  suppressPackageStartupMessages({
    # library(Biostrings)
    library(Azimuth)
    library(SeuratData)
    library(Seurat)
    library(tidyverse)
  })
  
  # reference must be a named list with path or any references from the following ref$Dataset
  if(isTRUE(avail.refs)){
    ref <- AvailableData() %>% data.frame %>% .[grep("Azimuth", .$Summary),]
    return(ref)
  }
  
  # run Azimuth
  for(i in 1:length(reference)){
    
    ref_suffix = paste0(".", names(reference)[i])
    umap.name = paste0("umap", ref_suffix)
    obj <- RunAzimuth(obj, assay = "RNA", reference = reference[i], umap.name = umap.name, ...)
    
    # extract and rename columns for cell type labels
    cols <- c(grep("^predicted.celltype..*[1-9]$", colnames(obj@meta.data), value=TRUE),
              grep("^predicted.celltype..*score$", colnames(obj@meta.data), value=TRUE),
              "mapping.score")
    # rename columns
    for(col in cols){
      obj[[paste0(col, ref_suffix)]] <- obj[[col]]
      obj[[col]] <- NULL
    }
    
  }
  
  # save cell labels
  ref_suffix <- paste0(".", names(reference), collapse = "|") 
  save_name <- paste0(".", names(reference), collapse = "|") %>% 
    gsub("\\.","",.) %>% gsub("\\|","_",.)
  pred_cell_labels <- obj@meta.data[grep(ref_suffix, colnames(obj@meta.data), value = TRUE)]

  ifelse(!dir.exists(file.path(save.loc, "int_obj")),
         dir.create(file.path(save.loc, "int_obj"), recursive = TRUE), paste0(save.loc," directory exists"))


  saveRDS(pred_cell_labels, file.path(save.loc, "int_obj", paste0(obj@project.name,".",save_name,".predicted_cell_labels.azimuth.rds")))
  
  return(obj)
  
}

