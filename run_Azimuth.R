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
  for(i in seq_along(reference)){
    
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
  
  # make output directory
  dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)
  
  # save cell labels
  ref_suffix <- paste0(".", names(reference), collapse = "|") 
  save_name <- paste0(".", names(reference), collapse = "|") %>% 
    gsub("\\.","",.) %>% gsub("\\|","_",.)
  pred_cell_labels <- obj@meta.data[grep(ref_suffix, colnames(obj@meta.data), value = TRUE)]
  saveRDS(pred_cell_labels, file.path(save.loc, "int_obj", paste0(obj@project.name,".",save_name,".predicted_cell_labels.azimuth.rds")))
  
  return(obj)
  
}

