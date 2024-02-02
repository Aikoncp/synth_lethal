

library(tictoc)
library(tidyverse)
tic()
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

drug_targets_path <- file.path(datapath, 'drug_targets.csv')
drug_targets <- read_csv(drug_targets_path) |>
  mutate(genes = str_split(genes, ","))


parallel <- TRUE
SL_pairs_method <- "SynLethDB"

if (parallel) {
  library(doParallel)
  cl <- makeCluster(max(detectCores() - 2, 1))
  registerDoParallel(cl)
}
plots <- list()
trials_scores <- list()
for (i in 2:nrow(drug_targets)) {
  #First must be fixed because it is two drugs
  drug <- drug_targets$drug[i]
  genes <- drug_targets$genes[i][[1]]
  print(drug)
  print(genes)
  #make list of all matching dataset files and make for loop for all paths
  
  if (i != 22 &
      i != 23) {
    #Tipifarnib has problems with the patient indices
    trial_scores <-
      calculate_score(drug, genes, SL_pairs_method, parallel)
    print(i)
    
    
    
    
    
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
          "Plot of the dataset ",
          drug,
          "(Synthetic lethal pairs obtained via the method ",
          SL_pairs_method,
          ")"
        )
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("Response to treatment") + ylab("Synthetic Lethality Score")
    
    
    

    #Tests 
    
    trial_scores  |>   #Check normality (Comprova que funcioni)
      group_by(score)  |> 
      shapiro.test(Response)
    
    
    trial_scores |>       #Check outliers
      group_by(Response) |> 
      identify_outliers(score)
    
    trial_scores |>  levene_test(score ~ Response)  #Check equal variances
    
    test <- t.test(score ~ Response, data = trial_scores, var.equal = TRUE)
    
    
    
    
    tests[[drug]] <- test
    plots[[drug]] <- plot
    trials_scores[[drug]] <- trial_scores
  }
}

if (parallel)
  stopCluster(cl)

toc()

#return(plots)