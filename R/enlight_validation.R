

datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

drug_targets_path <- file.path(datapath, 'drug_targets.csv')
drug_targets <- read_csv(drug_targets_path)

for (i in 1:length(drug_targets$drug)) {
  drug <- drug_targets[i]$drug
  gene <- drug_targets[i]$gene
  #make list of all matching dataset files and make for loop for all paths
  calculate_score(drug, gene)
  
  
  
  
}