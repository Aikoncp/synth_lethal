#' Title
#'
#' @param n A number. How many synthetic lethal pairs for each gene
#' @param SL_pairs_method 
#'
#' @return
#' @export
#'
#' @examples
filter_pairs <-
  function(SL_pairs, n = 100,
           SL_pairs_method = "binarize_expression") {
    
    datapath <-
      "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
    drug_targets_path <- file.path(datapath, 'drug_targets.csv')
    targeted_genes <- read_csv(drug_targets_path) |>
      transmute(genes = str_split(genes, ",")) |>
      unlist() |>
      unique()
    
    # if (SL_pairs_method == "SynLethDB") {
    #   SL_pairs_path <- file.path(datapath, 'Human_SL.csv')
    # } else if (SL_pairs_method == "binarize_expression") {
    #   SL_pairs_path <-
    #     file.path(datapath, 'SL_pairs_binarize_expression.csv')
    # } else if (SL_pairs_method == "binarize_essentiality") {
    #   SL_pairs_path <-
    #     file.path(datapath, 'SL_pairs_binarize_essentiality.csv')
    #   
    # }
    # SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)
    
##This is for if the first column is not targeted gene
    
    #   SL_pairs_1 <- SL_pairs_unfiltered |>
    #   filter(n1.name %in% targeted_genes) 
    # 
    # 
    # SL_pairs_2 <- SL_pairs_unfiltered |>
    #   filter(n2.name %in% targeted_genes) |>
    #   rename(foo = n1.name) |>
    #   rename(n1.name = n2.name, n2.name = foo)
    # 
    # SL_pairs_unfiltered <- bind_rows(SL_pairs_1, SL_pairs_2) |>  #Put all the targeted_genes in the first row and remove duplicates
    #   filter(!duplicated(paste0(pmax(n1.name, n2.name), pmin(n1.name, n2.name))))
      

    # q <- p.adjust(SL_pairs_unfiltered$p_value, method="fdr")
    # 
    # SL_pairs_unfiltered <- SL_pairs_unfiltered |> add_column(q_value = q) 
    # 
    # 
 
    
    SL_pairs_filtered <-  SL_pairs |> dplyr::filter(q_value < 0.05, depletion_q_value < 0.05, phylo_coefficient < 10.5, phylo_coefficient > 0) |>  #USe the same boundary for the phylogenetic coefficient as ISLE
      group_by(gene1) |> 
      slice_min(survival_coef, n = n) |>  #The ones with the most negative
    ungroup()
    write_csv(SL_pairs_filtered, "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data/SL_pairs_filtered.csv")
    
    SR_pairs_filtered <-  SL_pairs |> dplyr::filter(SR_DD_q_value < 0.05, SR_DD_depletion_q_value < 0.05) |> 
      group_by(gene1) |> 
      slice_max(survival_coef, n = n) |>  #The ones with the most negative 
    ungroup()
    write_csv(SR_pairs_filtered, "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data/SR_pairs_filtered.csv")
 
    
    #Possibly add other ways of filtering
    
    
    
    
    
    #getting only the 25 with less p-value for each gene
    return(list("SR_pairs_filtered" = SR_pairs_filtered, "SL_pairs_filtered" = SL_pairs_filtered))
    
    
  }