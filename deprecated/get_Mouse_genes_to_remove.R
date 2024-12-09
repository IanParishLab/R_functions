get_Mouse_genes_to_remove <- function(rna_counts){
  library(tidyverse)
  library(AnnotationHub)
  
  mito_rp_genes <- unique(c(grep("^mt-", rownames(rna_counts), v=T), grep("^Mtmr", rownames(rna_counts), v=T),grep("^Rp[sl]", rownames(rna_counts), v=T)))
  tcr_bcr_genes <-c(grep("^Tr[abdg][vjc]", rownames(rna_counts),value=T), grep("^Ig[hkl]v | ^Ig[hkl]j | ^Ig[kl]c | ^Igh[adegm]", rownames(rna_counts),value=T))
  
  # Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:
  
  cell_cycle_genes <- read.csv(text = RCurl::getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") )
  
  # Connect to AnnotationHub
  ah <- AnnotationHub()
  
  # Access the Ensembl database for organism
  ahDb <- query(ah, 
                pattern = c("Mus musculus", "EnsDb"), 
                ignore.case = TRUE)
  
  # Acquire the latest annotation files
  id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n = 1)
  
  # Download the appropriate Ensembldb database
  edb <- ah[[id]]
  
  # Extract gene-level information from database
  annotations <- genes(edb, 
                       return.type = "data.frame")
  
  # Select annotations of interest
  annotations <- annotations %>%
    dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
  
  # Get gene names for Ensembl IDs for each gene
  cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
  
  # Acquire the S phase genes
  s_genes <- cell_cycle_markers %>%
    dplyr::filter(phase == "S") %>%
    pull("gene_name")
  
  # Acquire the G2M phase genes        
  g2m_genes <- cell_cycle_markers %>%
    dplyr::filter(phase == "G2/M") %>%
    pull("gene_name")
  
  mouseGenes <- setNames(mito_rp_genes,tcr_bcr_genes,cell_cycle_markers$gene_name, c("mito_rp_genes","tcr_bcr_genes","cell_cycle_markers"))
  
  return(mouseGenes)
}