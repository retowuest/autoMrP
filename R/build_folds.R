

build_folds <- function(survey, 
                        L2.unit,  
                        ebma.size = 1/3, 
                        k.folds = 5, 
                        cv.sampling = "L2 units"){
  # EBMA hold-out fold
  ebma.size <- round(nrow(survey) * ebma.size, digits = 0)
  
  if(ebma.size>0){
    ebma_folding_out <- ebma_folding(data = survey,
                                     L2.unit = L2.unit,
                                     ebma.size = ebma.size)
    ebma_fold <- ebma_folding_out$ebma_fold
    cv_data <- ebma_folding_out$cv_data
  } else{
    ebma_fold <- NULL
    cv_data <- survey
  }
  
  # K folds for cross-validation
  cv_folds <- cv_folding(data = cv_data,
                         L2.unit = L2.unit,
                         k.folds = k.folds,
                         cv.sampling = cv.sampling)
} 

