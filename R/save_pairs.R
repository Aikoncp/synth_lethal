library(tictoc)
library(progressr)
library(tidyverse)
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
parallel <- FALSE
if (parallel) {
  library(foreach)
  library(doParallel)
  cl <- makeCluster(max(detectCores() - 5, 1))
  registerDoParallel(cl)
}

if (!exists("gene_effect")) {
  
  essentiality_path <-
    file.path(datapath, 'CRISPR_gene_effect.csv')
  expression_path <- file.path(datapath, 'CCLE_expression.csv')
  
  gene_effect <- read_csv(essentiality_path)
  expression <- read_csv(expression_path)
  
  cell_lines_1 <- gene_effect$DepMap_ID
  cell_lines_2 <- expression[[1]]
  cell_lines <- intersect(cell_lines_1, cell_lines_2)
  
  genes <- colnames(gene_effect)
  
  gene_effect <-
    gene_effect[gene_effect$DepMap_ID %in% cell_lines, ]
  gene_effect <- gene_effect[order(gene_effect$DepMap_ID), ]
  expression <- expression[expression[[1]] %in% cell_lines, ]
  expression <- expression |> rename(DepMap_ID = ...1)
  expression <- expression[order(expression$DepMap_ID), ]
}

SL_pairs_binarize_expression_path <-
  file.path(datapath, 'SL_pairs_binarize_expression.csv')
tic()
with_progress(SL_pairs_binarize_expression <- obtain_pairs("binarize_expression", n=10))
toc()
write_csv(SL_pairs_binarize_expression, SL_pairs_binarize_expression_path)



SL_pairs_binarize_essentiality_path <-
  file.path(datapath, 'SL_pairs_binarize_essentiality.csv')
tic()
SL_pairs_binarize_essentiality <- obtain_pairs("binarize_essentiality", max_p_value = max_p_value)
toc()
write_csv(SL_pairs_binarize_essentiality, SL_pairs_binarize_essentiality_path)