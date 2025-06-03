get_unique_genes <- function(gene_lists) {
  # Flatten the list into a single vector
  all_genes <- unlist(gene_lists)
  
  # Count occurrences of each gene
  gene_counts <- table(all_genes)
  
  # Filter genes that appear only once
  unique_genes <- names(gene_counts[gene_counts == 1])
  
  # Find which list contains each unique gene
  result <- lapply(gene_lists, function(gene_list) {
    intersect(gene_list, unique_genes)
  }) %>% unlist
  
  return(result)
}
# # example usage
# SR <- readRDS("FindMarkers.p14.clustered.RNA.LR.rds")
# SR <- split(SR, SR$cluster)
# SRgn <- lapply(SR, function(df){ df%>% filter(avg_log2FC > 3 & p_val_adj < 0.05) %>% pull(var = gene) })

# get_unique_clones <- function(clon) {
#   # Flatten the list into a single vector
#   all_clones <- unlist(clon)
#   
#   # Count occurrences of each gene
#   clone_counts <- table(all_clones)
#   
#   # Filter genes that appear only once
#   unique_clones <- names(clone_counts[clone_counts == 1])
#   
#   # Find which list contains each unique gene
#   result <- lapply(clon, function(clon_list) {
#     intersect(clon_list, unique_clones)
#   }) %>% unlist
#   
#   return(result)
# }

