


compute_essentiality_association <- function(datapath, outpath = ".", method = "ISLE", n = 10, tissue = 'None', bimodal = FALSE) {
  essentiality_path <- file.path(datapath, 'CRISPR_gene_effect.csv')
  expression_path <- file.path(datapath, 'CCLE_expression.csv')
  
  gene_effect <- read_csv(essentiality_path)
  expression <- read_csv(expression_path)
  
  cell_lines_1 <- gene_effect$DepMap_ID
  cell_lines_2 <- expression[[1]]
  cell_lines <- intersect(cell_lines_1, cell_lines_2)
  # bimodals <- read.csv('../data/bimodal.csv')
  
  genes <- colnames(gene_effect)
  # bimodal_filtered <- c('DepMap_ID')
  # bimodals <- bimodals$Symbol
  # for (gene in genes) {
  #   spliced <- sub(' \\(.*', '', gene)
  #   if (spliced %in% bimodals) {
  #     bimodal_filtered <- c(bimodal_filtered, gene)
  #   }
  # }
  # 
  # if (bimodal) {
  #   gene_effect <- gene_effect[bimodal_filtered]
  # }
  
  gene_effect <- gene_effect[gene_effect$DepMap_ID %in% cell_lines, ]
  expression <- expression[expression[[1]] %in% cell_lines, ]
  
  # if (tissue != "None") {
  #   genes_types <- tissue_specific_info(datapath)
  #   if (tissue == 'pan') {
  #     cat('pan\n')
  #     filtered_genes <- genes_types[genes_types$primary_disease == 'Pancreatic Cancer', ]
  #   } else if (tissue == 'lung') {
  #     filtered_genes <- genes_types[genes_types$primary_disease == 'Lung Cancer', ]
  #     cat('lung\n')
  #   } else if (tissue == 'breast') {
  #     cat('breast\n')
  #     filtered_genes <- genes_types[genes_types$primary_disease == 'Breast Cancer', ]
  #   } else {
  #     cat('invalid tissue type! Available: "None"/"pan"/"lung"/"breast"\n')
  #     return(NULL)
  #   }
  
  # filtered <- filtered_genes[, 1]
  # 
  # gene_effect <- gene_effect[gene_effect$DepMap_ID %in% filtered, ]
  # expression <- expression[expression$X1 %in% filtered, ]
  
  if (method == "ISLE") {
    wilcoxon_results <- perform_isle_method(gene_effect, expression, n)
    json_name <- paste("wilcoxon_results_", tissue, "_n=", n, ".json", sep = '')
    write_json(wilcoxon_results, file = file.path(outpath, json_name))
    return(wilcoxon_results)
  }
} else {  # all tissues
  if (method == "ISLE") {
    wilcoxon_results <- perform_isle_method(gene_effect, expression, n)
    json_name <- paste("wilcoxon_results_pancan_n=", n, ".json", sep = '')
    write_json(wilcoxon_results, file = file.path(outpath, json_name))
    return(wilcoxon_results)
  }
}

cat('done!\n')

return(NULL)
}

perform_isle_method <- function(gene_effect, expression, n, useParallel = FALSE) {
  total_cell_lines <- nrow(gene_effect)
  cat(total_cell_lines, '\n')
  benchmark <- total_cell_lines %/% 3
  
  if (n == "all") {
    total_gene_As <- ncol(gene_effect)
    total_gene_Bs <- ncol(expression)
  } else {
    total_gene_As <- n + 1
    total_gene_Bs <- n + 1
  }
  
  if (!useParallel) {
    wilcoxon_results <- list()
    for (i in 2:total_gene_As) {
      if (i %% 10 == 0) {
        cat(i, '\n')
      }
      x <- gene_effect[, c(1, i)]
      gene_A <- colnames(x)[2]
      wilcoxon_results[[gene_A]] <- list()
      for (j in 2:total_gene_Bs) {
        y <- expression[, c(1, j)]
        gene_B <- colnames(y)[2]
        y <- y[order(y[[2]], decreasing = FALSE), ]
        bottom_third <- y[1:benchmark, 1]
        top_third <- y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
        bottom_third <- bottom_third[[1]]
        top_third <- top_third[[1]]
        essentiality_top <- x[x$DepMap_ID %in% top_third, 2]
        essentiality_bottom <- x[x$DepMap_ID %in% bottom_third, 2]
        essentiality_top <- as.numeric(unlist(essentiality_top))
        essentiality_bottom <- as.numeric(unlist(essentiality_bottom))
        wilcoxon <- wilcox.test(essentiality_bottom, essentiality_top)$p.value
        wilcoxon_results[[gene_A]][[gene_B]] <- wilcoxon
      }
    }
    wilcoxon_results_df <- melt(wilcoxon_results)
    return(list(wilcoxon_results, wilcoxon_results_df))
  } else {  # Use parallel
    # Implement parallel processing here
    return(NULL)
  }
}

isle_submethod <- function(gene_effect, expression, i, j, resDict) {
  y <- expression[, c(1, j)]
  gene_B <- colnames(y)[2]
  y <- y[order(y[, 2]), ]
  bottom_third <- y[1:benchmark, 1]
  top_third <- y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
  bottom_third <- as.character(bottom_third)
  top_third <- as.character(top_third)
  essentiality_top <- x[x$DepMap_ID %in% top_third, 2]
  essentiality_bottom <- x[x$DepMap_ID %in% bottom_third, 2]
  essentiality_top <- as.numeric(essentiality_top)
  essentiality_bottom <- as.numeric(essentiality_bottom)
  wilcoxon <- wilcox.test(essentiality_bottom, essentiality_top)
  resDict[[gene_A]][[gene_B]] <- wilcoxon
}

tissue_specific_info <- function(datapath) {
  sample_info_path <- file.path(datapath, 'sample_info.csv')
  info <- load_dataset(sample_info_path)
  genes_types <- info[, c('DepMap_ID', 'primary_disease')]
  return(genes_types)
}

# Main
datapath <- '../depmap_data/'
essentiality_path <- file.path(datapath, 'CRISPR_gene_effect.csv')
expression_path <- file.path(datapath, 'CCLE_expression.csv')
cat(essentiality_path, '\n')
compute_essentiality_association(datapath, "ISLE", 10, 'pan', TRUE)
