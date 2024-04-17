library(tidyverse)


TCGA_path <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data/TCGA"


files <- list.files(path = TCGA_path, pattern = ".rds")

expression <- tibble()
survival <- tibble()


for (file in files) {
  print(file)
if (file != "TCGA_UCEC.rds"){ #Too messed up
    
    dat <- readRDS(file.path(TCGA_path, file))
    
    expression_i <-
      assay(dat[[grep("RNASeq2GeneNorm", names(experiments(dat)))[1]]]) #Get the RNASeq2GeneNorm expression matrix (only get the first onte hat is so)
    
    expression_i <-
      as_tibble(expression_i,  rownames = "gene")
   
    
    old_genes <- expression_i$gene
    
    survival_i <- as_tibble(colData(dat))
    
    colnames(survival_i) <-
      str_replace_all(colnames(survival_i), ".x", "") #Some datasets do weird things I don't know
    survival_i <- as_tibble(survival_i, .name_repair = "unique")
    
    if (file == "TCGA_SKCM.rds"){
      survival_i <- survival_i |> mutate(tumor_tissue_site = "melanoma")
    }
    # else if(file == "TCGA_UCEC.rds") {
    #   survival_i <- survival_i |> mutate(tumor_tissue_site = "endometrium", years_to_birth = patient.age_at_initial_pathologic_diagnosis)
    # }
    
    if(!'pathology_T_stage' %in% names(survival_i)){
      survival_i <- survival_i |>  add_column(pathology_T_stage = NA)
    }
    if(!'pathology_N_stage' %in% names(survival_i)){
      survival_i <- survival_i |>  add_column(pathology_N_stage = NA)
    }
    if(!'pathology_M_stage' %in% names(survival_i)){
      survival_i <- survival_i |>  add_column(pathology_M_stage = NA)
    }
    
  
    
    survival_i <- survival_i |> mutate(file = file)
    

    survival_i <-
      survival_i  |> select(
        c(
          patientID,
          years_to_birth,
          vital_status,
          days_to_death,
          days_to_last_followup,
          tumor_tissue_site,
          gender,
          pathology_T_stage,
          pathology_N_stage,
          pathology_M_stage,
          file
          
          
        )
      )
    
    if (nrow(expression) == 0) {
      expression <- expression_i
      survival <- survival_i
    } else {
      expression <-
        expression |> add_column(expression_i |> select(-gene))
      survival <- survival |> add_row(survival_i)
      
    }
    if ((file != files[1])  &&
        (identical(old_genes, as.vector(expression_i$gene)))) {
      #check the genes are the same for each file, otherwise big problems
      print("toma")
    }
}
  
}



#Sometimes there's repeated expressions and no repeated in survival


duplicates <- duplicated(str_sub(colnames(expression), 1, 12))
#I will remove any repeats completely, both from expression and survival


to_be_removed <-
  as.character(unique(colnames(expression)[duplicates]))
expression <-
expression|> select(-any_of(to_be_removed))

colnames(expression) <- str_sub(colnames(expression), 1, 12)
expression <- as_tibble(expression)

survival <- survival  |> filter(patientID %in% colnames(expression))


#All of the time columns must be joined into one, the information is in the status. Also days being 0 does not work with many models, and changing it by half should not have a significant effect on its validity
survival <- survival |> mutate(days = coalesce(days_to_death,days_to_last_followup)) |> dplyr::mutate(days = ifelse(days == 0, 0.5, days)) |> mutate(gender = as.factor(gender))#At 0 all patients must be alive for the models to work

#Remove patients with negative days or na days
#negative_days <- survival |> dplyr::filter(is.na(days) | days < 0) |> select(patientID)
#expression <- expression |> select(-any_of(negative_days$patientID))

survival <- survival |> dplyr::filter(!is.na(days), days > 0)

#Remove infrequent tissue sites
#low_tissue_site <- survival |> dplyr::filter(is.na(days) | days < 0) |> select(patientID)
survival <- survival |> mutate(tumor_tissue_site = if_else(is.na(tumor_tissue_site), "missing", tumor_tissue_site)) |>  group_by(tumor_tissue_site) |> filter(n() > 10) 
survival <- ungroup(survival)


#Standardize the stages:
survival <- survival |> mutate(pathology_T_stage = str_sub(pathology_T_stage,1,2)) |> mutate(pathology_T_stage = if_else(is.na(pathology_T_stage) | pathology_T_stage == "ti" | pathology_T_stage == "tx", "missing", pathology_T_stage)) 

survival <- survival |> mutate(pathology_N_stage = str_sub(pathology_N_stage,1,2)) |> mutate(pathology_N_stage = if_else(is.na(pathology_N_stage) | pathology_N_stage == "nx", "missing", pathology_N_stage)) 

survival <- survival |> mutate(pathology_M_stage = str_sub(pathology_M_stage,1,2)) |> mutate(pathology_M_stage = if_else(is.na(pathology_M_stage) | pathology_M_stage == "mx", "missing", pathology_M_stage)) 


expression <- expression |> select(gene, all_of(survival$patientID))







#Make them unordered factors

gender <- factor(survival$gender, ordered = FALSE)
tumor_tissue_site <- factor(survival$tumor_tissue_site, ordered = FALSE)
pathology_T_stage <- factor(survival$pathology_T_stage, ordered = FALSE)
pathology_N_stage <- factor(survival$pathology_N_stage, ordered = FALSE)
pathology_N_stage <- factor(survival$pathology_M_stage, ordered = FALSE)

survival <- survival |> mutate(gender = gender, tumor_tissue_site = tumor_tissue_site)


datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

TCGA_expression_path <-
  file.path(datapath, 'TCGA_expression.csv')
write_csv(expression, TCGA_expression_path)

TCGA_survival_path <-
  file.path(datapath, 'TCGA_survival.csv')
write_csv(survival, TCGA_survival_path)







#Create an expression file with the proportion of the expression of each gene for each patient


#expression_percentile <- colnames(expression) |>  purrr::map_dfc(setNames, object = list(numeric()))

# for (i in nrows(expression)){
#   distr <- ecdf(dplyr::slice(expression,i)[-1] |> as.numeric()) #take out the gene
#   new_row <- distr(dplyr::slice(expression,i)[-1] |> as.numeric())
#   new_row <- c(dplyr::slice(expression,i)[1],new_row)
#
# }



t_expression <- expression |>
  pivot_longer(cols = -gene, names_to = 'patientID') |>
  pivot_wider(names_from = gene, values_from = value)

t_expression <- as_tibble(t_expression)



TCGA_t_expression_path <-
  file.path(datapath, 'TCGA_t_expression.csv')
write_csv(t_expression, TCGA_t_expression_path)

i <- 0
#Create a transposed expression tibble with the percentile of the gene expression as values. This can later be used both for continuous and discrete expression values
for(gene in colnames(t_expression)[-1]) {
i <- i + 1
  distr <- ecdf(select(t_expression, {{gene}}) |> unlist() |> as.numeric())
  t_expression <- t_expression |>  mutate(!!gene := distr(select(t_expression, {{gene}}) |> unlist() |> as.numeric()))
if (i %% 100 == 0){
  print (i)
}

}



TCGA_t_expression_percentile_path <-
  file.path(datapath, 'TCGA_t_expression_percentile.csv')
write_csv(t_expression, TCGA_t_expression_percentile_path)