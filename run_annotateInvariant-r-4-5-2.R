# args <- commandArgs(trailingOnly = TRUE)
# 
# # Quick check to ensure an input file was provided
# if (length(args) == 0) {
#   stop("Error: No input file provided. Usage: ~/R/functions/run_annotateInvariant-r-4-5-2.R <input_file.qs>", call. = FALSE)
# }
# 
# suppressPackageStartupMessages({
#   library(scRepertoire)
#   library(immApex)
#   library(janitor)
#   library(Seurat)
#   library(tidyverse)
#   library(qs2)
# })
# 
# seu <- qs_read(args[1])
# seu <- annotateInvariant(seu, type = "MAIT", species = 'human')
# seu <- annotateInvariant(seu, type = "iNKT", species = 'human')
# 
# print(tabyl(seu@meta.data, MAIT.score))
# print(tabyl(seu@meta.data, iNKT.score))