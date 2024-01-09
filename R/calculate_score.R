calculate_score <- function(drug_name, targeted_gene) {
  # drug_name <- "Anti-PD1_4"
  # targeted_gene <- "PDCD1"
    
  datapath <-
    "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
  
  SL_pairs_path <- file.path(datapath, 'Human_SL.csv')
  SL_pairs <- read_csv(SL_pairs_path)
  drug_responses_path <-
    file.path(datapath,
              'enlight-data-main',
              "drug_response_classifications.csv")
  drug_responses <- read_csv(drug_responses_path)
  
  
  trial_expression_path <-
    file.path(datapath, 'enlight-data-main', str_c(drug_name, ".csv"))
  trial_expression <-  read_csv(trial_expression_path)

  
  SL_pairs_drug1 <- SL_pairs |>
    filter(n1.name == targeted_gene) |>
    rename(SL_gene = n2.name) |>
    select(SL_gene, r.statistic_score)
  
  
  SL_pairs_drug2 <- SL_pairs |>
    filter(n2.name == targeted_gene) |>
    rename(SL_gene = n1.name) |>
    select(SL_gene, r.statistic_score)
  
  SL_pairs_drug <-
    bind_rows(SL_pairs_drug1, SL_pairs_drug2) #All the SL pairs of the gene that the drug targets
  SL_pairs_drug_list <- SL_pairs_drug[["SL_gene"]]
  
  drug_responses_drug <- drug_responses |>
    filter(Dataset == str_to_title(drug_name)) |>
    rename(sample_ID = `Sample ID`)
  
  ###Calculate each patient's score
  
  avg_gene_expression <- trial_expression |>
    filter(index %in% SL_pairs_drug_list) |>
    rowwise() |>
    transmute(
      SL_gene = index,
      upper_third = quantile(c_across(where(is.numeric)), 0.66),
      lower_third = quantile(c_across(where(is.numeric)), 0.33)
    )
  
  score_list <- list()
  for (patient in as.list(drug_responses_drug$sample_ID)) {
    score <- 0
    for (gene in as.list(avg_gene_expression$SL_gene)) {
      gene_expression_value <- trial_expression |>
        filter(index == gene) |>
        select(patient) |>
        pull()
      
      if (gene_expression_value < (avg_gene_expression |> filter(SL_gene == gene) |> select(lower_third) |> pull())) {
        score <- score + 1
      } else if (gene_expression_value > avg_gene_expression |> filter(SL_gene == gene) |> select(upper_third) |> pull()) {
        score <- score - 1
      }
      
      
    }
    #Add score to patient
    score <- score / length(avg_gene_expression$SL_gene)
    print(score)
    score_list <- append(score_list, score)
  }
  

  drug_responses_drug <- drug_responses_drug |>
    mutate(score = as.numeric(score_list))
  
  
  
  # ggplot(data = drug_responses_drug, aes(x = score, fill = Response)) +
  #   geom_histogram(binwidth = 2,
  #                  alpha = 0.5,
  #                  position = "identity") +
  #   theme_minimal()
  
  drug_responses_drug |>
    group_by(Response) |>
    ggplot(aes(
      x = Response,
      y = score,
      group = Response,
      fill = Response
    )) +
    geom_boxplot()
  
  
  
 
  # wilcoxon <-
  #   wilcox.test(as.numeric(drug_responders$score),
  #               as.numeric(drug_non_responders$score))
  
}

#retornar el que minteressi perfer numeros