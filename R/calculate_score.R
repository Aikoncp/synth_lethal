#' Calculate the score of an ENLIGHT trial
#'
#' @param drug_name The trial that will be scored
#' @param targeted_genes The genes targeted by the drug used in the trial
#' @param parallel Whether to parallelize the algorithm
#' @param synthetic_rescue Whether to also use synthetic rescue pairs
#' @param SL_pairs_method The method used to obtain the synthetic lethal pairs
#'
#' @return
#' @export
#'
#' @examples
calculate_score <-
  function(drug_name,
           targeted_genes,
           SL_pairs_method = "SynLethDB",
           parallel = FALSE,
           synthetic_rescue = TRUE,
           SL_pairs = tibble(),
           SR_pairs = tibble()) {

    # drug_name <- "Anti-PD1_4"
    # targeted_genes <- c("PDCD1", "CD274")
    
    
    
    datapath <-
      "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
# 
#     if (SL_pairs_method == "SynLethDB") {
#       SL_pairs_path <- file.path(datapath, 'Human_SL.csv')
#     } else if (SL_pairs_method == "binarize_expression") {
#       SL_pairs_path <-
#         file.path(datapath, 'SL_pairs_binarize_expression.csv')
#     } else if (SL_pairs_method == "binarize_essentiality") {
#       SL_pairs_path <-
#         file.path(datapath, 'SL_pairs_binarize_essentiality.csv')
#     } else {
#       SL_pairs_path <-
#         file.path(datapath, "SL_pairs_filtered.csv")
#     }
#     
#     SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)
#     
#     if (synthetic_rescue) {
#       #SR_pairs_path <- file.path(datapath, 'Human_SR.csv')
#       SR_pairs_path <- file.path(datapath, 'SR_pairs_filtered.csv')
#       SR_pairs <- read_csv(SR_pairs_path, show_col_types = FALSE)
#     }
    
    drug_responses_path <-
      file.path(datapath,
                'enlight-data-main',
                "drug_response_classifications.csv")
    drug_responses <-
      read_csv(drug_responses_path, show_col_types = FALSE)
    
    
    trial_expression_path <-
      file.path(datapath, 'enlight-data-main', str_c(drug_name, ".csv"))
    trial_expression <-
      read_csv(trial_expression_path, show_col_types = FALSE)
    

    
    drug_responses_drug <- drug_responses |>
      filter(Dataset == drug_name) |>
      rename(sample_ID = `Sample ID`)
    
    n_patients <- nrow(drug_responses_drug)
    
    final_scores <- numeric(n_patients)
    
    for (targeted_gene in targeted_genes) {
      #print(targeted_gene)
      
      # SL_pairs_drug1 <- SL_pairs |>
      #   filter(gene1 == targeted_gene) |>
      #   rename(SL_gene = gene2) |>
      #   select(SL_gene, p_value)
      # 
      # 
      # SL_pairs_drug2 <- SL_pairs |>
      #   filter(gene2 == targeted_gene) |>
      #   rename(SL_gene = gene1) |>
      #   select(SL_gene, p_value)
      
      SL_pairs_drug <- SL_pairs |> filter(gene1 == targeted_gene) |> select(gene2, p_value) |> rename(SL_gene = gene2) #Take the pairs containing the targeted gene and keep its pair
       # bind_rows(SL_pairs_drug1, SL_pairs_drug2) #All the SL pairs of the gene that the drug targets

      
      
      
      if(synthetic_rescue){
      SR_pairs_drug <- SR_pairs |> filter(gene1 == targeted_gene) |> select(gene2, p_value) |> rename(SR_gene = gene2)

      }
      
      ###Calculate each patient's score
      
      SL_avg_gene_expression <- trial_expression |>
        filter(index %in% SL_pairs_drug$SL_gene  | index == targeted_gene) |>
        rowwise() |>
        transmute(
          SL_gene = index,
          upper_third = quantile(c_across(where(is.numeric)), 0.66),
          lower_third = quantile(c_across(where(is.numeric)), 0.33)
        )
      
      if(synthetic_rescue){
      SR_avg_gene_expression <- trial_expression |>
        filter(index %in% SR_pairs_drug$SR_gene) |>
        rowwise() |>
        transmute(
          SR_gene = index,
          upper_third = quantile(c_across(where(is.numeric)), 0.66),
          lower_third = quantile(c_across(where(is.numeric)), 0.33)
        )
      }
      
      if (parallel) {
        scores <-
          foreach (i = 1:n_patients, .packages = "tidyverse") %dopar% {
            score <- 0
            patient <- drug_responses_drug$sample_ID[i]
            
            if(targeted_gene %in% trial_expression$index){
            targeted_gene_expression <- trial_expression |>
              filter(index == targeted_gene) |>
              select(all_of(patient)) |>
              pull()
            }
            
            if (!(targeted_gene %in% trial_expression$index) ||
                targeted_gene_expression > SL_avg_gene_expression |> filter(SL_gene == targeted_gene) |> select(lower_third) |> pull()) {
              #If the targeted gene is not expressed then the drug will not have any effect regardless. if we dontknow we let through??
              for (gene in as.list(SL_avg_gene_expression$SL_gene)) {
                gene_expression_value <- trial_expression |>
                  filter(index == gene) |>
                  select(all_of(patient)) |>
                  pull()
                
                
                
                if (!is.na(gene_expression_value) && gene_expression_value < (SL_avg_gene_expression |> filter(SL_gene == gene) |> select(lower_third) |> pull())) {
                  score <- score + 1
                } else if (gene_expression_value > (SL_avg_gene_expression |> filter(SL_gene == gene) |> select(upper_third)) |> pull()) {
                  #score <- score - 1   ENLIGHT does not do this
                }
                
                
                
                
              }
              
              for (gene in as.list(SR_avg_gene_expression$SR_gene)) {
                if (synthetic_rescue &&
                    gene_expression_value > (SR_avg_gene_expression |> filter(SR_gene == gene) |> select(upper_third)) |> pull()) {
                  score <- score + 1
                }
              }
            }
            
            score #Ara mateix si un gen te moltes synthetic lethal pairs cada un aper individual compta menys
          }
        
      }
      
      else {
        scores <- foreach (i = 1:n_patients) %do% {
          score <- 0
          patient <- drug_responses_drug$sample_ID[i]
          for (gene in as.list(SL_avg_gene_expression$SL_gene)) {
            gene_expression_value <- trial_expression |>
              filter(index == gene) |>
              select(patient) |>
              pull()
            
            if (gene_expression_value < (
              SL_avg_gene_expression |> filter(SL_gene == gene) |> select(lower_third) |> pull()
            )) {
              score <- score + 1
            } else if (gene_expression_value > SL_avg_gene_expression |> filter(SL_gene == gene) |> select(upper_third) |> pull()) {
              score <- score - 1
            }
            
            
          }
          #Add score to patient
          
          print(score)
          score #Ara mateix si un gen te moltes synthetic lethal pairs cada un aper individual compta menys
        }
        
        
        
      }
      final_scores <-  unlist(scores) + final_scores
      
    }
    
    
   
    
    
    #### Normalization of the scores
    #final_scores <- final_scores / length(targeted_genes)
    final_scores <- (final_scores + abs(min(final_scores)))
    final_scores <- final_scores / max(final_scores, 1)
    
    drug_responses_drug <- drug_responses_drug |>
      mutate(score = final_scores)
    
    #return(drug_responses_drug)
    
    # ggplot(data = drug_responses_drug, aes(x = score, fill = Response)) +
    #   geom_histogram(binwidth = 2,
    #                  alpha = 0.5,
    #                  position = "identity") +
    #   theme_minimal()
    
   
    
    
    
    return(drug_responses_drug)
    
    # wilcoxon <-
    #   wilcox.test(as.numeric(drug_responders$score),
    #               as.numeric(drug_non_responders$score))
    
  }

#retornar el que minteressi perfer numeros