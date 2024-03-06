#' Predict using a pmodel
#' 
#' pmodelApply takes an existing pmodel and applies it to an input set of samples. 
#' 
#' @param mymodel   An existing pmodel
#' @param exprData  Expression data matrix
#' 
#' @returns         A vector of predicted efficacies.
#' @export
#' @import caret pROC

pmodelApply <- function(mymodel, exprData){
  
  genes <- colnames(mymodel$trainingData)[1:length(mymodel$trainingData) - 1]
  
  x <- t(as.matrix(exprData$ex)[genes, ]) 

  x_hat_pre <- predict(mymodel, x)
  
  res <- exprData$res
  
  obs <- factor(res$response)
    
  eval <- pROC::roc(obs, x_hat_pre)
  print(eval)
  plot(eval)
    
  out <- list()
  out$eval <- eval
  out$predictions <- x_hat_pre
  return(out)
}