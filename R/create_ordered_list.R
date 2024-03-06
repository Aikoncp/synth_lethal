expression_genes <- word(colnames(expression))[-1] #Remove the DepMapID
effect_genes <- word(colnames(gene_effect))[-1]

top_thirds <- lst()
bottom_thirds <- lst()
expression_tops <- lst()
expression_bottoms <- lst()


for (gene in expression_genes) {
  
  y <- expression |> select(DepMap_ID, {{gene}})
  
  bottom_third <-
    y |>  slice_min(order_by = .data[[gene]], n = benchmark)
  top_third <-
    y |>  slice_max(order_by = .data[[gene]], n = benchmark)
  
  
  top_thirds[[gene]] <- top_third
  bottom_thirds[[gene]] <- bottom_third
    
}