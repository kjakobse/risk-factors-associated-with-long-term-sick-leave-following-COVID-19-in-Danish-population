DiscreteCovariateNames <- function(covariates, 
                                   discrete_covariates = NULL) {
  if(!is.character(covariates))  stop("covariates must be a character vector")
  
  if(is.null(discrete_covariates)) {
    return(discrete_covariates)
  } else {
    if(!is.character(discrete_covariates)) {
      stop("discrete_covariates must be NULL or a character vector")
    }
    
    return(str_subset(covariates, str_c("^(", 
                                        str_c(discrete_covariates,
                                              collapse = "|"),
                                        ")")))
  }
}