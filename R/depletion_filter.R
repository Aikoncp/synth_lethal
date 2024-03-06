#Take in the existing pairs, the ccle expression, and perform hyoergeometric test on the pairs, and add the p-value to the tibble. 
library(tidyverse)

datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
drug_targets_path <- file.path(datapath, 'drug_targets.csv')
targeted_genes <- read_csv(drug_targets_path) |>
  transmute(genes = str_split(genes, ",")) |>
  unlist() |>
  unique()


expression_path <- file.path(datapath, 'CCLE_expression.csv')

if (SL_pairs_method == "SynLethDB") {
  SL_pairs_path <- file.path(datapath, 'Human_SL.csv')
} else if (SL_pairs_method == "binarize_expression") {
  SL_pairs_path <-
    file.path(datapath, 'SL_pairs_binarize_expression.csv')
} else if (SL_pairs_method == "binarize_essentiality") {
  SL_pairs_path <-
    file.path(datapath, 'SL_pairs_binarize_essentiality.csv')
  
}
SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)

if (!exists("expression")) {
  gene_effect <- read_csv(essentiality_path)
  expression <- read_csv(expression_path)
  
  cell_lines_1 <- gene_effect$DepMap_ID
  cell_lines_2 <- expression[[1]]
  cell_lines <- intersect(cell_lines_1, cell_lines_2)
  
  expression <- expression[expression[[1]] %in% cell_lines, ]
  expression <- expression |> rename(DepMap_ID = ...1)
  expression <- expression[order(expression$DepMap_ID), ]
  colnames(expression) <- word(colnames(expression))
  expression <- as_tibble(expression)
  
}

for row in 
slice_min()

overlap <- intersect(cell_line_low_expressed_A, cell_line_low_expressed_B)
depletion_p_value <- phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE)  #If the p-value is small there has been a significant depletion


