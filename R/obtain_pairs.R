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
    
    expression_genes <- colnames(expression)[-1] #Remove the DepMapID
    effect_genes <- colnames(gene_effect)[-1]
    
    n_expression_genes <- ncol(expression)
    
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
      tibble(
        gene1 = character(),
        gene2 = character(),
        p_value = numeric(),
        depletion_p_value = numeric()
      )
    
    
    p <- progressor(length(targeted_genes) + 1)
    
    
    top_thirds <- lst()
    bottom_thirds <- lst()
    expression_tops <- lst()
    expression_bottoms <- lst()
    
    
    thirds_path <-
      file.path(datapath, 'thirds.RData')
    
    load(thirds_path) #Get the bottom_thirds and top_thirds
    
    p()
    if (!exists("bottom_thirds") || !exists("top_thirds")) {
      for (gene in expression_genes) {   #Guardaho que mai canviara
        
        y <- expression |> select(DepMap_ID, {{gene}})
        
        bottom_third <-
          y |>  slice_min(order_by = .data[[gene]], prop = 0.33)  #It will get more than the benchmark if the minimum value is the same for more than benchmark cell lines (0) Better than the alternative because otherwise it would be more arbitrary
        top_third <-
          y |>  slice_max(order_by = .data[[gene]], prop = 0.33) |> filter(.data[[gene]] > 0)  #given that sometimes anything above 0 is selected, bimodality could be implemented
        
        
        
        top_thirds[[gene]] <- top_third
        bottom_thirds[[gene]] <- bottom_third
        
      }
    }
    
    
    if (method == "binarize_expression") {
      for (gene_A in targeted_genes) {
        
        p()
        
        
        
        if (gene_A %in% effect_genes) {
          x <- gene_effect |> select(DepMap_ID, {{gene_A}})
          
          
          if (gene_A %in% colnames(expression)) {
            
            bottom_third_A <- bottom_thirds[[gene_A]]
          }
          #Only check gene A because it is the gene that has been "removed", which is what a drug would do too (VEGFA included because it is not in gene_effect for some reason)
          for (gene_B in expression_genes) {
  
           
            
            bottom_third <- bottom_thirds[[gene_B]]
            top_third <- top_thirds[[gene_B]]
            
            
            essentiality_bottom <-
              x |>  dplyr::filter(DepMap_ID %in% bottom_third$DepMap_ID) |> select({
                {
                  gene_A
                }
              }) |>  unlist() |> as.numeric()
            essentiality_top <-
              x |>  dplyr::filter(DepMap_ID %in% top_third$DepMap_ID) |> select({
                {
                  gene_A
                }
              }) |> unlist() |> as.numeric()
            
            #Do the wilcoxon test
            if(length(top_third$DepMap_ID) > 0) {
              wilcoxon_sl <-
              wilcox.test(essentiality_bottom,
                          essentiality_top,
                          alternative = "greater")$p.value #Check whether essentiality_bottom is greater than essentiality_top (when gene y is missing the cell is more likely to die when gene x is removed)
            
            }
            else {
              wilcoxon_sl <- NA
              }
            #Do the hypergeometric depletion test
            
            if (gene_A %in% colnames(expression)) {
              overlap <-
                length(intersect(
                  bottom_third$DepMap_ID,
                  bottom_third_A$DepMap_ID
                ))
              
              #phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE) (the depletion test)
              
              depletion_p_value <-
                phyper(overlap,
                       length(bottom_third$DepMap_ID),
                       total_cell_lines - length(bottom_third$DepMap_ID),
                       length(bottom_third_A$DepMap_ID),
                       lower.tail = TRUE)
              
            }
            else {
              depletion_p_value <- 2 #fesho millorrrrrrrrrrrrrrrrrrrrrr
            }
            
            if (!is.nan(wilcoxon_sl)) {
              SL_pairs <- SL_pairs |>
                add_row(
                  gene1 = word(gene_A),
                  gene2 = word(gene_B),
                  p_value = wilcoxon_sl,
                  depletion_p_value = depletion_p_value
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
            i = 1:length(targeted_genes),
            .packages = c("tidyverse", "stats"),
            .combine = rbind
          ) %dopar% {
            
            p()
            
            gene_A <- targeted_genes[i]
            
            if (gene_A %in% effect_genes) {
              x <- gene_effect |> select(DepMap_ID, {{gene_A}})
              
              
              if (gene_A %in% colnames(expression)) {
                z <- expression |> select(DepMap_ID, {
                  {
                    gene_A
                  }
                })
                bottom_third_A <-
                  z |>  slice_min(order_by = .data[[gene_A]], n = benchmark)
              }
              for (gene_B in expression_genes) {
                
                y <- expression |> select(DepMap_ID, {{gene_B}})
                
                bottom_third <-
                  y |>  slice_min(order_by = .data[[gene_B]], n = benchmark)
                top_third <-
                  y |>  slice_max(order_by = .data[[gene_B]], n = benchmark)
                
                
                essentiality_bottom <-
                  x |>  dplyr::filter(DepMap_ID %in% bottom_third$DepMap_ID) |> select({
                    {
                      gene_A
                    }
                  }) |>  unlist() |> as.numeric()
                essentiality_top <-
                  x |>  dplyr::filter(DepMap_ID %in% top_third$DepMap_ID) |> select({
                    {
                      gene_A
                    }
                  }) |> unlist() |> as.numeric()
                
                #Do the wilcoxon test
                wilcoxon_sl <-
                  wilcox.test(essentiality_bottom,
                              essentiality_top,
                              alternative = "greater")$p.value #Check whether essentiality_bottom is greater than essentiality_top (when gene y is missing the cell is more likely to die when gene x is removed)
                
                
                #Do the hypergeometric depletion test
                
                if (gene_A %in% colnames(expression)) {
                  overlap <-
                    length(intersect(
                      bottom_third$DepMap_ID,
                      bottom_third_A$DepMap_ID
                    ))
                  depletion_p_value <-
                    phyper(overlap,
                           benchmark,
                           total_cell_lines - benchmark,
                           benchmark,
                           lower.tail = TRUE)
                  
                }
                else {
                  depletion_p_value <- 2 #fesho millorrrrrrrrrrrrrrrrrrrrrr
                }
                
                if (!is.nan(wilcoxon_sl)) {
                  SL_pairs_i <- SL_pairs |>
                    add_row(
                      gene1 = word(gene_A),
                      gene2 = word(gene_B),
                      p_value = wilcoxon_sl,
                      depletion_p_value = depletion_p_value
                    )
                }
                
                
              }
            }
            SL_pairs_i
          }
        
      }
    }
    
    SL_pairs <-
      SL_pairs[order(SL_pairs$p_value, decreasing = FALSE),]
    
    q_value <- p.adjust(SL_pairs$p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(q_value = q_value)
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_p_value = 1 - SL_pairs$p_value)
    
    SR_DD_q_value <- p.adjust(SL_pairs$SR_DD_p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_q_value = SR_DD_q_value)
    
    
    depletion_q_value <-
      p.adjust(SL_pairs$depletion_p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(depletion_q_value = depletion_q_value)
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_depletion_p_value = 1 - depletion_p_value)
    
    SR_DD_depletion_q_value <- p.adjust(SL_pairs$SR_DD_depletion_p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_depletion_q_value = SR_DD_depletion_q_value)
    
    
    return(SL_pairs)
  }