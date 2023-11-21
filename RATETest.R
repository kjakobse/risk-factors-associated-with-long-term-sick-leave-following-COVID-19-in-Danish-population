RATETest <- function(cf, 
                     cov,
                     cov_type = c("continuous", "discrete"),
                     target = "AUTOC",
                     q = seq(0.1, 1, by = 0.1),
                     R = 500, 
                     level = 0.95) {
  if (!hasArg(q) && cov_type[1] == "discrete") {
    q <- c(
      0.001,
      cumsum(rev(table(cf[["X.orig"]][[cov]]))) / 
        length(cf[["X.orig"]][[cov]])
    )
  }
  rate <- rank_average_treatment_effect(
    cf,
    cf[["X.orig"]][[cov]],
    q = q,
    target = target,
    R = R
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