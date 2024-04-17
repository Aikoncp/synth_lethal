create_network <-  ##Create a network of SL and SR candidates
  function(cancer_type = "pan-cancer",
           crispr_method = "binarize_expression", n = "all", synthetic_rescue = TRUE, save = TRUE) {
    
    library(progressr)
    
    print("Undergoing the CRISPR test:")
    SL_pairs <- with_progress(obtain_pairs(method = crispr_method, n = n, cancer_type = cancer_type))
    
    print("Undergoing the survival test:")
    SL_pairs <- with_progress(survival_filter(cancer_type = cancer_type))
    
    print("Undergoing the phylogenetic test:")
    SL_pairs <- with_progress(phylogetic_test())
    
    
    if(save){
    SL_pairs_path <-
      file.path(datapath, str_c("SL_pairs_", cancer_type, "_", crispr_method, "_", n, "_", ifelse(synthetic_rescue, "SR", ""), ".csv"))
    
    write_csv(SL_pairs, SL_pairs_path)
    }
    
    return(SL_pairs)
    
  }