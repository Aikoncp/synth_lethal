evaluate_network <- function(SL_pairs, tissue = "all"){
  
  pairs <- filter_pairs(SL_pairs)
  
  results <- enlight_validation(SL_pairs = pairs$SL_pairs_filtered, SR_pairs = pairs$SR_pairs_filtered, tissue = tissue)
  
  
} 