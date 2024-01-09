#' Create a synthetic lethal network
#' 
#' createEssentialNetwork takes as arguments a path to Cancer Dependency Map essentiality data
#' and parameters defining the network properties. It then computes a gene essentiality 
#' network. All the variations produce a test statistic and a p-value.
#' 
#' @param datapath  Path to directory containing Dependency Map data
#' @param outpath   Output path to save .rds object, optionally including filename.
#' @param method    Character string describing which method to use. One of c("ISLE")
#' @param n         The number of genes to compute associations for, useful for debugging and testing.
#' Default is 0, which computes the entire network. 
#' @param tissue    Character string defining which tissue to use. Default is "all"
#' 
#' @returns         A dataframe containing a matrix of effect sizes (effectSizes), a matrix of p-values
#' (pvals), and a list of network parameters. 
#' @export
#' @importFrom stats wilcox.test
#' 

createEssentialNetwork <- function(datapath=".", outpath=".", method="ISLE", n=0, 
                                   tissue="all"){
  


  

    essentiality_path <- file.path(datapath, 'CRISPR_gene_effect.csv')
    expression_path <- file.path(datapath, 'CCLE_expression.csv')
    
    gene_effect <- read_csv(essentiality_path)
    expression <- read_csv(expression_path)
    
    cell_lines_1 <- gene_effect$DepMap_ID
    cell_lines_2 <- expression[[1]]
    cell_lines <- intersect(cell_lines_1, cell_lines_2)
    
    genes <- colnames(gene_effect)

    gene_effect <- gene_effect[gene_effect$DepMap_ID %in% cell_lines, ]
    expression <- expression[expression[[1]] %in% cell_lines, ]
    max_p_value <- 0.05



  
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
    
  
      wilcoxon_results <- tibble(gene1 = character(), gene2 = character(), p_value = numeric())
      for (i in 2:total_gene_As) {
        if (i %% 10 == 0) {
          cat(i, '\n')
        }
        x <- gene_effect[, c(1, i)]
        gene_A <- colnames(x)[2]
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
          if (wilcoxon < max_p_value) {
            wilcoxon_results <- wilcoxon_results |> 
              add_row(gene1 = gene_A, gene2 = gene_B, p_value = wilcoxon)
          }
        }
      }
      wilcoxon_ordered <- wilcoxon_results[order(wilcoxon_results$p_value), ]

     
  }


  

  
  
  
}