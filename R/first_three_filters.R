
library(tidyverse)
datapath <- "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

essentiality_path <- file.path(datapath, 'CRISPR_gene_effect.csv')
expression_path <- file.path(datapath, 'CCLE_expression.csv')

gene_effect <- read_csv(essentiality_path)
expression <- read_csv(expression_path)

cell_lines_1 <- gene_effect$DepMap_ID
cell_lines_2 <- expression[[1]]
cell_lines <- intersect(cell_lines_1, cell_lines_2)

genes <- colnames(gene_effect)

gene_effect <- gene_effect[gene_effect$DepMap_ID %in% cell_lines, ]
gene_effect <- gene_effect[order(gene_effect$DepMap_ID), ]
expression <- expression[expression[[1]] %in% cell_lines, ]
expression <- expression |> rename(DepMap_ID = ...1) #
expression <- expression[order(expression$DepMap_ID), ]




max_p_value <- 0.05
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



if (method == "binarize_expression") {
  for (i in 2:total_gene_As) {
    if (i %% 10 == 0) {
      cat(i, '\n')
    }
    x <- gene_effect[, c(1, i)]
    gene_A <- colnames(x)[2]
    for (j in 2:total_gene_Bs) {
      y <- expression[, c(1, j)]
      gene_B <- colnames(y)[2]
      y <- y[order(y[[2]], decreasing = FALSE),]
      bottom_third <- y[1:benchmark, 1]
      top_third <-
        y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
      bottom_third <- bottom_third[[1]]
      top_third <- top_third[[1]]
      essentiality_top <- x[x$DepMap_ID %in% top_third, 2]
      essentiality_bottom <- x[x$DepMap_ID %in% bottom_third, 2]
      essentiality_top <- as.numeric(unlist(essentiality_top))
      essentiality_bottom <- as.numeric(unlist(essentiality_bottom))
      wilcoxon <-
        wilcox.test(essentiality_bottom, essentiality_top)$p.value
      if (wilcoxon < max_p_value) {
        SL_pairs <- SL_pairs |>
          add_row(gene1 = gene_A,
                  gene2 = gene_B,
                  p_value = wilcoxon)
      }
    }
  }
} else if (method == 'binarize_essentiality') {
  for (i in 2:total_gene_As) {
    if (i %% 10 == 0) {
      cat(i, '\n')
    }
    x <- expression[, c(1, i)]
    gene_A <- colnames(x)[2]
    for (j in 2:total_gene_Bs) {
      y <- gene_effect[, c(1, j)]
      gene_B <- colnames(y)[2]
      y <- y[order(y[[2]], decreasing = FALSE),]
      bottom_third <- y[1:benchmark, 1]
      top_third <-
        y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
      bottom_third <- bottom_third[[1]]
      top_third <- top_third[[1]]
      expression_top <- x[x[[1]] %in% top_third, 2]
      expression_bottom <- x[x[[1]] %in% bottom_third, 2]
      expression_top <- as.numeric(unlist(expression_top))
      expression_bottom <- as.numeric(unlist(expression_bottom))
      wilcoxon <-
        wilcox.test(expression_bottom, expression_top)$p.value
      if (wilcoxon < max_p_value) {
        SL_pairs <- SL_pairs |>
          add_row(gene1 = gene_A,
                  gene2 = gene_B,
                  p_value = wilcoxon)
      }
    }
  }
  #Continuous
} else {
  for (i in 2:total_gene_As) {
    if (i %% 10 == 0) {
      cat(i, '\n')
    }
    x <- gene_effect[, c(1,i)]
    gene_A <- colnames(x)[2]
    for (j in 2:total_gene_Bs) {
      y <- expression[, c(1, j)]
      gene_B <- colnames(y)[2]
      cor <- cor.test(as.numeric(unlist(x[2])), as.numeric(unlist(y[2])), 
               method = "pearson")$estimate
  
      
      if (cor > 0.3) {
        SL_pairs <- SL_pairs |>
          add_row(gene1 = gene_A,
                  gene2 = gene_B,
                  p_value = cor)
        plot_data <- tibble(ess = x[[2]], expr = y[[2]])
        ggplot(
          data = plot_data,
          mapping = aes(x = expr, y = ess)
        ) +
          geom_point()


      }
    }
  }
}

SL_pairs_ordered <-
  SL_pairs[order(SL_pairs$p_value, decreasing = TRUE),]