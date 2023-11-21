CausalForestDynamicSubgroups <- function(cf,
                                         num_rankings = 3,
                                         num_folds = 5,
                                         ...) {
  # save dots arguments
  args <- list(...)
  # Methods for forest with and without clustering
  if (length(cf$clusters) > 0) {
    # partition data
    folds <- rep(0, length(cf$clusters))
    for (i in unique(cf$clusters)) {
      folds[cf$clusters == i] <- sample(sort(seq_along(folds[cf$clusters == i]) %% 
                                               num_folds) + 1)
    }
    
    n <- length(cf$Y.orig)
    indices <- split(seq(n), folds)
    result <- map_dfr(indices, function(idx, args) {
      # Fit outcome model to predict on held-out data. Note fitting a causal 
      # forest only gives outcome predictions on data included in training data.
      forest_m <- do.call(
        regression_forest,
        c(
          list(
            X = cf$X.orig[-idx,], 
            Y = cf$Y.orig[-idx], 
            clusters = cf$clusters[-idx]
          ),
          args
        )
      )
      m_hat <- predict(forest_m, cf$X.orig[idx,])$predictions
      
      # Fit exposure model to predict on held-out data. Note fitting a causal 
      # forest only gives exposure predictions on data included in training data.
      forest_e <- do.call(
        regression_forest,
        c(
          list(
            X = cf$X.orig[-idx,], 
            Y = cf$W.orig[-idx], 
            clusters = cf$clusters[-idx]
          ),
          args
        )
      )
      e_hat <- predict(forest_e, cf$X.orig[idx,])$predictions
      
      # train forest with a held-out fold and original clustering
      cf_rank <- do.call(
        causal_forest,
        c(
          list(
            X = cf$X.orig[-idx,], 
            Y = cf$Y.orig[-idx], 
            W = cf$W.orig[-idx], 
            Y.hat = forest_m$predictions,
            W.hat = forest_e$predictions,
            clusters = cf$clusters[-idx]
          ),
          args
        )
      )
      
      # Estimate cate's in held-out fold
      tau_hat <- predict(object = cf_rank, newdata = cf$X.orig[idx,])$predictions
      
      # aipw scores in held-out fold
      mu_hat_0 <- m_hat - e_hat * tau_hat       
      mu_hat_1 <- m_hat + (1 - e_hat) * tau_hat
      
      aipw_scores <- 
        tau_hat +
        cf$W.orig[idx] / e_hat * (cf$Y.orig[idx] -  mu_hat_1) -
        (1 - cf$W.orig[idx]) / (1 - e_hat) * (cf$Y.orig[idx] - mu_hat_0)
      
      # rank observations by cate in held-out fold
      tau_hat_quantiles <- quantile(x = tau_hat, 
                                    probs = seq(0, 1, by = 1/num_rankings)
      )
      # if quantiles are not unique, manually sort and cut into appropriate 
      # groups
      if (length(tau_hat_quantiles) == length(unique(tau_hat_quantiles))) {
        ranking <- cut(
          x = tau_hat, 
          breaks = tau_hat_quantiles, 
          include.lowest = TRUE,
          labels = seq_len(num_rankings)
        )
      } else {
        len <- length(tau_hat)
        ranking <- tau_hat |>
          (\(x) {
            tibble(
              id = seq_along(x),
              tau_hat = x
            )
          }
          )() |>
          arrange(tau_hat) |>
          mutate(
            id_2 = seq_along(tau_hat),
            rank = cut(x = id_2,
                       breaks = c(seq(0, len %% num_rankings, by = 1) *
                                    (len %/% num_rankings + 1),
                                  seq(len %% num_rankings + 1, num_rankings,
                                      length.out = num_rankings - len %% num_rankings) *
                                    len %/% num_rankings + len %% num_rankings),
                       include.lowest = TRUE,
                       labels = seq_len(num_rankings))
          ) |>
          arrange(id) |>
          pull(rank)
      }
      
      # collect ranking and aipw scores
      res <- tibble(id = idx,
                    tau_hat = tau_hat,
                    ranking = ranking,
                    aipw_scores = aipw_scores)
      
      res
    }, 
    args)
    result <- result |> arrange(id)
  } else {
    # partition data into folds
    folds <- sort(seq_along(cf$Y.orig) %% num_folds) + 1
    
    # train forest using folds as clusters
    cf_rank <- do.call(
      causal_forest,
      c(
        list(
          X = cf$X.orig, 
          Y = cf$Y.orig, 
          W = cf$W.orig, 
          clusters = folds
        ),
        args
      )
    )
    
    # estimate cate's. Note that the cluster containing the sample to predict is
    # left out for prediction.
    tau_hat <- predict(object = cf_rank)$predictions
    
    # rank observations within folds
    ranking <- rep(NA, length(cf$Y.orig))
    for (fold in seq_len(num_folds)) {
      tau_hat_quantiles <- quantile(
        x = tau_hat[folds == fold], 
        probs = seq(0, 1, by = 1/num_rankings)
      )
      # if quantiles are not unique, manually sort and cut into appropriate 
      # groups
      if (length(tau_hat_quantiles) == length(unique(tau_hat_quantiles))) {
        ranking[folds == fold] <- cut(
          x = tau_hat[folds == fold], 
          breaks = tau_hat_quantiles, 
          include.lowest = TRUE,
          labels = seq_len(num_rankings)
        )
      } else {
        len <- length(tau_hat[folds == fold])
        ranking[folds == fold] <- tau_hat[folds == fold] |>
          (\(x) {
            tibble(
              id = seq_along(x),
              tau_hat = x
            )
          }
          )() |>
          arrange(tau_hat) |>
          mutate(
            id_2 = seq_along(tau_hat),
            rank = cut(x = id_2,
                       breaks = c(seq(0, len %% num_rankings, by = 1) *
                                    (len %/% num_rankings + 1),
                                  seq(len %% num_rankings + 1, num_rankings,
                                      length.out = num_rankings - len %% num_rankings) *
                                    len %/% num_rankings + len %% num_rankings),
                       include.lowest = TRUE,
                       labels = seq_len(num_rankings))
          ) |>
          arrange(id) |>
          pull(rank)
      }
    }
    
    # aipw scores
    mu_hat_0 <- cf_rank$Y.hat - cf_rank$W.hat * tau_hat       
    mu_hat_1 <- cf_rank$Y.hat + (1 - cf_rank$W.hat) * tau_hat
    
    aipw_scores <- 
      tau_hat + 
      cf$W.orig / cf_rank$W.hat * (cf$Y.orig -  mu_hat_1) - 
      (1 - cf$W.orig) / (1 - cf_rank$W.hat) * (cf$Y.orig - mu_hat_0)
    
    # collect ranking and aipw scores
    result <- tibble(id = seq_along(cf$Y.orig),
                     tau_hat = tau_hat,
                     ranking = ranking,
                     aipw_scores = aipw_scores)
  }
  
  # fit linear model of aipw scores to find average in each rank
  ols <- lm(result$aipw_scores ~ 0 + factor(result$ranking))
  cf_rank_ate <- tibble(method = "aipw", 
                        ranking = str_c("Q", seq_len(num_rankings)), 
                        estimate = coeftest(ols, vcov=vcovHC(ols, "HC3"))[,1],
                        std_err = coeftest(ols, vcov=vcovHC(ols, "HC3"))[,2])
  
  # plot with estimates and 95 % confidence intervals within each ranking:
  cf_rank_ate_plot <- ggplot(cf_rank_ate, aes(x = ranking, y = estimate)) + 
    geom_point() +
    geom_errorbar(aes(ymin = estimate - 2 * std_err, 
                      ymax = estimate + 2 * std_err), 
                  width = 0.2) +
    ylab("") + xlab("") +
    ggtitle("AIPW score within each ranking (as defined by predicted CATE)")
  
  # table with tests for differences between ranking groups
  res <- tibble()
  for (i in seq_len(num_rankings - 1)) {
    lev <- seq_len(num_rankings)
    ols <- lm(result$aipw_scores ~ 1 + factor(result$ranking,
                                              levels = c(lev[i], lev[-i])))
    res <- as_tibble(
      coef(summary(ols))[seq(i + 1, num_rankings), 
                         c(1, 2, 4), 
                         drop = FALSE]
    ) |>
      mutate(id = paste("Rank",  seq(i + 1, num_rankings), "- Rank ", i)) |> 
      (\(x) rbind(res, x))()
  }
  
  # Adjust for multiple testing using the Benjamini-Hockberg procedure
  res <- res |>
    rename(`Orig. p-value` = `Pr(>|t|)`) |>
    mutate(
      `95 % CI` = paste0("(",
                         sprintf(
                           fmt = "%.3f",
                           Estimate - qnorm(0.975) * `Std. Error`
                         ),
                         ", ",
                         sprintf(
                           fmt = "%.3f",
                           Estimate + qnorm(0.975) * `Std. Error`
                         ),
                         ")"),
      `Orig. p-value` = sprintf(
        fmt = "%.3f",
        `Orig. p-value`
      ),
      `Adj. p-value` = sprintf(
        fmt = "%.3f",
        p.adjust(`Orig. p-value`, method = "BH")
      )
    ) |>
    select(id, Estimate, `Std. Error`, `95 % CI`, everything())
  
  # Plot heatmap with average of covariates in each group
  df <- map_dfr(
    colnames(cf$X.orig),
    function(covariate) {
      # beregn gennemsnit og standardafvigelse for hver kovariat
      fmla <- formula(str_c("`", covariate, "`", "~ 0 + ranking"))
      data_cf <- as_tibble(cf$X.orig) |>
        mutate(
          Y = cf$Y.orig,
          W = cf$W.orig,
          ranking = factor(result$ranking)
        )
      ols <- lm(fmla, data = data_cf)
      ols_res <- coeftest(ols, vcov = vcovHC(ols, "HC3"))
      
      # results
      avg <- ols_res[, 1]
      stderr <- ols_res[, 2]
      
      # collect results in table
      tibble(
        covariate = covariate,
        avg = avg,
        stderr = stderr,
        ranking = paste0("Q", seq_len(num_rankings)),
        scaling = pnorm((avg - mean(avg)) / sd(avg)),
        variation = sd(avg) / sd(data_cf[[!!covariate]]),
        labels = paste0(sprintf("%.3f", avg), 
                        "\n", 
                        "(", 
                        sprintf("%.3f", stderr), 
                        ")")
      )
    }
  ) |>
    mutate(covariate = fct_reorder(covariate, variation))
  
  # plot heatmap
  heatmap <- ggplot(df) +
    aes(ranking, covariate) +
    geom_tile(aes(fill = scaling)) + 
    geom_text(aes(label = labels)) +
    scale_fill_gradient(low = "#0000FF", high = "#FF8000") +
    ggtitle("Average covariate values within group (based on CATE estimate ranking)") +
    theme_minimal() + 
    ylab("") + xlab("CATE estimate ranking") +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11))
  
  return(
    list(
      cf_subgroups = result,
      cf_rank_ate = cf_rank_ate,
      cf_rank_ate_plot = cf_rank_ate_plot,
      cf_rank_diff_test = res,
      heatmap_data = df,
      heatmap = heatmap
    )
  )
}