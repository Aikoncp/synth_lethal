

library(progressr)
library(foreach)
library(doFuture)
registerDoFuture() ## tell foreach to use futures
plan(multisession, workers = max(detectCores() - 5, 1))




# xs <- 1:5
# 
# with_progress({
#   p <- progressor(along = xs) ## create a 5-step progressor
#   y <- foreach(x = xs) %dopar% {
#     p()                       ## signal a progression update
#     Sys.sleep(6.0 - x)
#     sqrt(x)
#   }
# })
# 
# 
# n <- 10
# SL_pairs <-
#   foreach (
#     i = 2:n,
#     .packages = c("tidyverse", "stats"),
#     .combine = progress_bar(n)
#   ) %dopar% {
#     if (i %% 1 == 0) {
#       print(i)
#     }
#   }

with_progress({
  p <- progressor(total_gene_As - 1)
  SL_pairs <-
    foreach (
      i = 2:total_gene_As,
      .packages = c("tidyverse", "stats"),
      .combine = rbind
    ) %dopar% {
      p(sprintf("i=%g", i))
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
          y <- y[order(y[[2]], decreasing = FALSE), ]
          bottom_third <- y[1:benchmark, 1]
          top_third <-
            y[(total_cell_lines - benchmark + 1):total_cell_lines, 1]
          bottom_third <- bottom_third[[1]]
          top_third <- top_third[[1]]
          essentiality_top <- x[x$DepMap_ID %in% top_third, 2]
          essentiality_bottom <-
            x[x$DepMap_ID %in% bottom_third, 2]
          essentiality_top <- as.numeric(unlist(essentiality_top))
          essentiality_bottom <-
            as.numeric(unlist(essentiality_bottom))
          wilcoxon <-
            wilcox.test(essentiality_bottom, essentiality_top, alternative = "greater")$p.value
          
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
})
