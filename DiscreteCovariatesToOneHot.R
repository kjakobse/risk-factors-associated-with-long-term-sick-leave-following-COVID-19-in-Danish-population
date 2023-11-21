DiscreteCovariatesToOneHot <- function(kovariater) {
  if(!is.data.frame(kovariater)) {
    stop("Covariates must be a data.frame or data.frame like object.")
  }
  
  for(i in seq_along(kovariater)) {
    if(!is.factor(kovariater[[i]])) {
      stop("Each covariate must be a factor.")
    }
  }
  X <- suppressMessages(map(names(kovariater),
                            ~ model.matrix(~ kovariater[[.x]] + 0)) %>%
                          map_dfc(~ as_tibble(.x)))
  
  X_names <- c()
  for(i in seq_along(names(kovariater))) {
    for(j in seq_along(levels(kovariater[[names(kovariater)[i]]]))) {
      X_names <- c(X_names, str_c(names(kovariater)[i],
                                  "_",
                                  levels(kovariater[[names(kovariater)[i]]])[j]))
    }
  }
  names(X) <- X_names
  return(X)
}