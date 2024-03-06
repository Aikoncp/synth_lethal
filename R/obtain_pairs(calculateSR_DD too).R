#' Obtain synthetic lethal pairs using the DepMap CRISPR datasets
#'
#' @param method A string. How should the pairs be computed: : `"binarize_expression"`, `'binarize_essentiality'`, or
#'
#' @param n A number. How many genes to use for the computation
#' @param max_p_value A number between 0 and 1. Below which p_value should the pairs be saved.
#'
#' @return A tibble with each pairs and its corresponding statistic
#' @export
#'
#' @examples
obtain_pairs <-
  function(method = "binarize_expression",
           n = "all",
           parallel = FALSE) {
    library(tidyverse)
    datapath <-
      "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
    
    essentiality_path <-
      file.path(datapath, 'CRISPR_gene_effect.csv')
    expression_path <- file.path(datapath, 'CCLE_expression.csv')
    
    
    
    if (!exists("gene_effect")) {
      gene_effect <- read_csv(essentiality_path)
      expression <- read_csv(expression_path)
      
      cell_lines_1 <- gene_effect$DepMap_ID
      cell_lines_2 <- expression[[1]]
      cell_lines <- intersect(cell_lines_1, cell_lines_2)
      
      
      
      gene_effect <-
        gene_effect[gene_effect$DepMap_ID %in% cell_lines, ]
      colnames(gene_effect) <- word(colnames(gene_effect))
      gene_effect <- as_tibble(gene_effect)
      gene_effect <- gene_effect[order(gene_effect$DepMap_ID), ]
      expression <- expression[expression[[1]] %in% cell_lines, ]
      expression <- expression |> rename(DepMap_ID = ...1)
      expression <- expression[order(expression$DepMap_ID), ]
      colnames(expression) <- word(colnames(expression))
      expression <- as_tibble(expression)
      
    }
    
    genes <- colnames(gene_effect)
    
    drug_targets_path <- file.path(datapath, 'drug_targets.csv')
    targeted_genes <- read_csv(drug_targets_path) |>
      transmute(genes = str_split(genes, ",")) |>
      unlist() |>
      unique()
    
    
    
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
    
    
    SL_pairs <-
      tibble(gene1 = character(),
             gene2 = character(),
             p_value = numeric())
    
    SR_DD_pairs <-
      tibble(gene1 = character(),
             gene2 = character(),
             p_value = numeric())
    
    
    if (method == "binarize_expression") {
      for (i in 2:total_gene_As) {
        if (i %% 100 == 0) {
          cat(i, '\n')
        }
        
        x <- gene_effect[, c(1, i)]
        gene_A <- colnames(x)[2]
        
        if (word(gene_A) %in% targeted_genes) {  #Only check gene A because it is the gene that has been "removed", which is what a drug would do too (VEGFA included because it is not in gene_effect for some reason)
         for (j in 2:total_gene_Bs) {
          y <- expression[, c(1, j)]
          gene_B <- colnames(y)[2]
          
            y <- y[order(y[[2]], decreasing = FALSE),]
            bottom_third <- y[1:benchmark, 1]
            top_third <-
              y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
            bottom_third <- bottom_third[[1]]
            top_third <- top_third[[1]]
            essentiality_top <-
              x[x$DepMap_ID %in% top_third, 2] #essentiality of gene x for when gene y is highly expressed
            essentiality_bottom <-
              x[x$DepMap_ID %in% bottom_third, 2]
            essentiality_top <- as.numeric(unlist(essentiality_top))
            essentiality_bottom <-
              as.numeric(unlist(essentiality_bottom))
            wilcoxon_sl <-
              wilcox.test(essentiality_bottom,
                          essentiality_top,
                          alternative = "greater")$p.value #Check whether essentiality_bottom is greater than essentiality_top (when gene y is missing the cell is more likely to die when gene x is removed)
            
            wilcoxon_sr_dd <-
              wilcox.test(essentiality_bottom,
                          essentiality_top,
                          alternative = "less")$p.value #Check whether essentiality_bottom is greater than essentiality_top (when gene y is missing the cell is more likely to die when gene x is removed)
            
            
            if (!is.nan(wilcoxon_sl)) {
              SL_pairs <- SL_pairs |>
                add_row(
                  gene1 = word(gene_A),
                  gene2 = word(gene_B),
                  p_value = wilcoxon_sl
                )
            }
            if (!is.nan(wilcoxon_sr_dd)) {
              SR_DD_pairs <- SR_DD_pairs |>
                add_row(
                  gene1 = word(gene_A),
                  gene2 = word(gene_B),
                  p_value = wilcoxon_sr_dd
                )
            }
            
          }
        }
      }
    } else if (method == 'binarize_essentiality') {
      for (i in 2:total_gene_Bs) {
        if (i %% 100 == 0) {
          cat(i, '\n')
        }
        x <- expression[, c(1, i)]
        gene_B <- colnames(x)[2]
        for (j in 2:total_gene_As) {
          y <- gene_effect[, c(1, j)]
          gene_A <- colnames(y)[2]
          if (word(gene_A) %in% targeted_genes ||
              word(gene_B) %in% targeted_genes) {
            y <- y[order(y[[2]], decreasing = FALSE),]
            bottom_third <- y[1:benchmark, 1]
            top_third <-
              y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
            bottom_third <- bottom_third[[1]]
            top_third <- top_third[[1]]
            expression_top <- x[x[[1]] %in% top_third, 2]
            expression_bottom <- x[x[[1]] %in% bottom_third, 2]
            expression_top <- as.numeric(unlist(expression_top))
            expression_bottom <-
              as.numeric(unlist(expression_bottom))
            wilcoxon <-
              wilcox.test(expression_bottom, expression_top, alternative = "greater")$p.value #check whether the expression of genes with low essentiality is higher than the expression of genes with high essentiality (if gene y is not essential to the cell, then gene x is probably very expressed)
            if (!is.nan(wilcoxon) && wilcoxon < max_p_value) {
              SL_pairs <- SL_pairs |>
                add_row(
                  gene1 = word(gene_A),
                  gene2 = (gene_B),
                  p_value = wilcoxon
                )
            }
          }
        }
      }
      #Continuous
    } else {
      for (i in 2:total_gene_As) {
        if (i %% 100 == 0) {
          cat(i, '\n')
        }
        x <- gene_effect[, c(1, i)]
        gene_A <- colnames(x)[2]
        for (j in 2:total_gene_Bs) {
          y <- expression[, c(1, j)]
          gene_B <- colnames(y)[2]
          if (word(gene_A) %in% targeted_genes ||
              #The genes in the DEPMAP have a different format, word extracts the first word
              (gene_B) %in% targeted_genes) {
            cor <- cor.test(as.numeric(unlist(x[2])), as.numeric(unlist(y[2])),
                            method = "pearson")$estimate
            
            
            if (cor > 0.3) {
              SL_pairs <- SL_pairs |>
                add_row(
                  gene1 = word(gene_A),
                  gene2 = word(gene_B),
                  p_value = cor
                )
              plot_data <- tibble(ess = x[[2]], expr = y[[2]])
              ggplot(data = plot_data,
                     mapping = aes(x = expr, y = ess)) +
                geom_point()
              
              
            }
          }
        }
      }
    }
    
    
    
    
    
    
    if (parallel) {
      if (method == "binarize_expression") {
        SL_pairs <-
          foreach (
            i = 2:total_gene_As,
            .packages = c("tidyverse", "stats"),
            .combine = rbind
          ) %dopar% {
            print("si")
            if (i %% 1 == 0) {
              cat(i, '\n')
            }
            SL_pairs_i <-
              tibble(gene1 = character(),
                     gene2 = character(),
                     p_value = numeric())
            x <- gene_effect[, c(1, i)]
            gene_A <- colnames(x)[2]
            for (j in 2:total_gene_Bs) {
              y <- expression[, c(1, j)]
              gene_B <- colnames(y)[2]
              if (word(gene_A) %in% targeted_genes ||
                  word(gene_B) %in% targeted_genes) {
                y <- y[order(y[[2]], decreasing = FALSE),]
                bottom_third <- y[1:benchmark, 1]
                top_third <-
                  y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
                bottom_third <- bottom_third[[1]]
                top_third <- top_third[[1]]
                essentiality_top <- x[x$DepMap_ID %in% top_third, 2]
                essentiality_bottom <-
                  x[x$DepMap_ID %in% bottom_third, 2]
                essentiality_top <-
                  as.numeric(unlist(essentiality_top))
                essentiality_bottom <-
                  as.numeric(unlist(essentiality_bottom))
                wilcoxon <-
                  wilcox.test(essentiality_bottom,
                              essentiality_top,
                              alternative = "greater")$p.value
                
                if (!is.nan(wilcoxon)) {
                  SL_pairs_i <- SL_pairs_i |>
                    add_row(
                      gene1 = word(gene_A),
                      gene2 = word(gene_B),
                      p_value = wilcoxon
                    )
                }
                
              }
            }
            SL_pairs_i
          }
        
      }
    }
    
    SL_pairs_ordered <-
      SL_pairs[order(SL_pairs$p_value, decreasing = FALSE),]
    
    SR_DD_pairs_ordered <-
      SR_DD_pairs[order(SR_DD_pairs$p_value, decreasing = TRUE),]
    
    pairs <- lst()
    pairs[["SL"]] <- SL_pairs_ordered
    pairs[["SR_DD"]] <- SR_DD_pairs_ordered  #The p-value is simply 1-, in the end you are looking for the same I think, so you dont need to do both, just get the DD
    
    return(pairs)
  }