enlight_validation <- function(SL_pairs, synthetic_rescue = TRUE, parallel = TRUE, SL_pairs_method = "filtered", tissue = "all"){
library(rstatix)
library(tictoc)
library(tidyverse)
tic()
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"


pairs <- filter_pairs(SL_pairs)

SL_pairs <- pairs$SL_pairs_filtered
SR_pairs <- pairs$SR_pairs_filtered



drug_targets_path <- file.path(datapath, 'drug_targets.csv')
drug_targets <- read_csv(drug_targets_path) |>
  mutate(genes = str_split(genes, ","))





if (parallel) {
  library(doParallel)
  cl <- makeCluster(max(detectCores() - 2, 1))
  registerDoParallel(cl)
}
plots <- list()
trials_scores <- list()
tests <-  list()
for (i in 2:nrow(drug_targets)) {
  #First must be fixed because it is two drugs
  drug <- drug_targets$drug[i]
  genes <- drug_targets$genes[i][[1]]
  targeted_tissue <- drug_targets$tissue[i]
  print(drug)
  print(genes)
  #make list of all matching dataset files and make for loop for all paths
  
  if (i != 22 &
      i != 23 & i!=6  & i!=7  & i!=8 & i!=9 & (tissue == "all" | tissue == targeted_tissue)) { #VEGFA DOESNT EXIST in the gene_effect, could make it so it considers VEGFA for expression
    #Tipifarnib has problems with the patient indices 
    trial_scores <-
      calculate_score(drug, genes, SL_pairs_method, parallel, synthetic_rescue = TRUE, SL_pairs = SL_pairs, SR_pairs = SR_pairs)
    print(i)
    
    
    test <- t.test(score ~ Response, data = trial_scores, alternative = "less", var.equal = TRUE)
    
    
    #Plot
    
    plot <- trial_scores |>
      group_by(Response) |>
      ggplot(aes(
        x = Response,
        y = score,
        group = Response,
        fill = Response
      )) +
      geom_boxplot() +
      stat_summary(
        fun.data = get_box_stats,
        geom = "text",
        hjust = 0.5,
        vjust = 0.9
      ) +
      geom_dotplot(binaxis = 'y',
                   stackdir = 'center',
                   dotsize = 1) +
      # geom_text(data = means, aes(label = weight, y = weight + 0.08)) +
      ggtitle(
        str_c(
          "Plot of the Clinical Trial ",
          drug, "    p-value = ", round(test$p.value,3)
        )
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      # labs(caption = 
      #      str_c("Synthetic lethal pairs obtained via the method ",
      #      "Binarize expression")) +
      xlab("Response to treatment") + ylab("Synthetic Lethality Score")
    
    
    print(plot)

    #Tests 

    # trial_scores  |>   #Check normality (Comprova que funcioni)
    #   group_by(Response)  |>
    #   shapiro_test(score)
    # 
    # 
    # trial_scores |>       #Check outliers
    #   group_by(Response) |>
    #   identify_outliers(score)
    # 
    # trial_scores |>  levene_test(score ~ Response)  #Check equal variances

   




    tests[[drug]] <- test
    plots[[drug]] <- plot
    trials_scores[[drug]] <- trial_scores
    
  }
}

if (parallel){
  stopCluster(cl)
}

toc()

return(list("tests" = tests, "plots" = plots, "trials_scores" = trials_scores))
}