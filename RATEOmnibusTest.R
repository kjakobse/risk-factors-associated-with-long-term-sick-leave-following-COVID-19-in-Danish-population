RATEOmnibusTest <- function(cf, 
                            level = 0.95,
                            target = "AUTOC", 
                            q = seq(0.1, 1, 0.1),
                            R = 500,
                            num_threads = 2,
                            sample.weights = NULL,
                            clusters = NULL,
                            ...) {
  require(policytree)
  sample <- sample(seq_along(cf$Y.orig), length(cf$Y.orig) / 2, replace = FALSE)
  cf_1 <- causal_forest(
    X = cf$X.orig[sample,],
    Y = cf$Y.orig[sample],
    W = cf$W.orig[sample],
    num.trees  = cf$`_num_trees`,
    num.threads = num_threads,
    ...
  )
  cf_2 <- causal_forest(
    X = cf$X.orig[-sample,],
    Y = cf$Y.orig[-sample],
    W = cf$W.orig[-sample],
    num.trees  = cf$`_num_trees`,
    num.threads = num_threads,
    ...
  )
  tau_1 <- predict(cf_2, newdata = cf_1$X.orig)$predictions
  tau_2 <- predict(cf_1, newdata = cf_2$X.orig)$predictions
  tau <- vector(mode = "double", length = length(cf$Y.orig))
  tau[sample] <- tau_1
  tau[-sample] <- tau_2
  tmp <- double_robust_scores(cf_1)
  dr_1 <- vector(mode = "double", length = length(cf_1$W.orig))
  for (i in seq_along(cf_1$W.orig)) {
    dr_1[i] <- tmp[i, cf_1$W.orig[i] + 1]
  }
  tmp <- double_robust_scores(cf_2)
  dr_2 <- vector(mode = "double", length = length(cf_2$W.orig))
  for (i in seq_along(cf_2$W.orig)) {
    dr_2[i] <- tmp[i, cf_2$W.orig[i] + 1]
  }
  dr <- vector(mode = "double", length = length(cf$Y.orig))
  dr[sample] <- dr_1
  dr[-sample] <- dr_2
  
  rate <- rank_average_treatment_effect.fit(
    DR.scores = dr,
    priorities = tau,
    target = target,
    q = q,
    R = R,
    sample.weights = sample.weights,
    clusters = clusters
  )
  confint <- rate$estimate + 
    tibble(
      estimate = 0,
      lower = -qnorm(1 - (1 - level) / 2) * rate$std.err,
      upper =  qnorm(1 - (1 - level) / 2) * rate$std.err
    )
  pval <- 2 * pnorm(-abs(rate$estimate) / rate$std.err)
  return(
    list(
      rate = rate,
      confint = confint,
      pval = pval
    )
  )
}