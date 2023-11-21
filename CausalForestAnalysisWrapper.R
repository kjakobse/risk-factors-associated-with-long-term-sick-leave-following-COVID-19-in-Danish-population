CausalForestAnalysisWrapper <- function (data, 
                                         cov = NULL, 
                                         binary = seq_along(cov), 
                                         ...) {
  if (!exists("DiscreteCovariatesToOneHot")) {
    stop ("The function DiscreteCovariatesToOneHot must be available.")
  }
  cov_bin <- cov[binary]
  if (is.null(binary) || (is.integer(binary) & length(binary) == 0)) {
    cov_dis <- cov
  } else {
    cov_dis <- cov[-binary]
  }
  return(
    causal_forest(
      X = data |>
        select(
          "age",
          "sex",
          "charlson_score_cat",
          "chronic_asthma",
          "chronic_diabetes",
          "chronic_high_blood_pressure",
          "chronic_copd_lung_disease",
          "chronic_headache",
          all_of(cov_bin)
        ) |>
        mutate(across(2:last_col(), ~ as.numeric(.x) - 1)) |>
        bind_cols(
          data |>
            select(
              "status_education",
              all_of(cov_dis)
            ) |>
            DiscreteCovariatesToOneHot()
        ),
      Y = as.numeric(data$sick_leave_outcome) - 1,
      W = as.numeric(data$test_result) - 1,
      ...
    )
  )
}