library(tidyverse)


datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"



CTRP_exps_path <- file.path(datapath, "CTRP-exps.csv")

CTRP_exps <- read_csv(CTRP_exps_path)

CTRP_exps <- CTRP_exps |> dplyr::select(treatmentid, sampleid, aac_recomputed, ic50_recomputed)
treatment_ids <- unique(CTRP_exps$treatmentid)
sample_ids <- unique(CTRP_exps$sampleid)

CTRP_exps <- CTRP_exps |> dplyr::filter(!is.na(ic50_recomputed), !is.na(aac_recomputed)) 