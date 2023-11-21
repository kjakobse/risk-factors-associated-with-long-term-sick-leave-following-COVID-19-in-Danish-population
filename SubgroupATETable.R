SubgroupATETable <- function(x, y, cf, level, subset = NULL) {
  require(stringr)
  require(dplyr)
  ci_names <- c(
    paste0(100 * level, "% CI - lower"),
    paste0(100 * level, "% CI - upper")
  )
  if(is.null(x) && is.null(y) && is.null(subset)) {
    ate <- average_treatment_effect(cf)
    return(
      tibble(
        subgroup = "Full population",
        n = length(cf$predictions),
        estimate = ate[1],
        !!sym(ci_names[1]) := ate[1] - qnorm(1 - (1 - level) / 2) * ate[2],
        !!sym(ci_names[2]) := ate[1] + qnorm(1 - (1 - level) / 2) * ate[2]
      )
    )
  } else if (is.null(subset)) {
    subset <- cf[["X.orig"]][[y]] %in% x
    ate <- average_treatment_effect(cf, subset = subset)
  } else {
    if (length(y) != 2) {
      stop(glue::glue("y must be a character vector of length 2, ",
                      "not {length(y)}, containing a column name ",
                      "and a name to use in the column"))
    }
    y <- c("custom", y)
    ate <- average_treatment_effect(cf, subset = subset)
  }
  tbl <- tibble(
    estimate = ate[1],
    n = sum(subset),
    !!sym(ci_names[1]) := ate[1] - qnorm(1 - (1 - level) / 2) * ate[2],
    !!sym(ci_names[2]) := ate[1] + qnorm(1 - (1 - level) / 2) * ate[2]
  )
  switch(
    EXPR = y[1],
    age = {
      tbl <- bind_cols(
        tibble(subgroup = paste0(y, " - ", x[1], "-", x[length(x)])),
        tbl
      )
    },
    sex = {
      tbl <- bind_cols(
        tibble(subgroup = ifelse(x == 0, "sex - female", "sex - male")),
        tbl
      )
    },
    charlson_score_cat = {
      tbl <- bind_cols(
        tibble(subgroup = case_when(
          x == 3 ~ paste0(y, " - 3+"),
          TRUE ~ paste0(y, " - ", as.character(x))
        )),
        tbl
      )
    },
    chronic_asthma = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(x == 0, "chronic asthma - no", "chronic asthma - yes")
        ),
        tbl
      )
    },
    chronic_diabetes = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(x == 0, "chronic diabetes - no", "chronic diabetes - yes")
        ),
        tbl
      )
    },
    chronic_high_blood_pressure = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(
            x == 0, 
            "chronic high blood pressure - no", 
            "chronic high blood pressure - yes"
          )
        ),
        tbl
      )
    },
    chronic_copd_lung_disease = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(
            x == 0, 
            "COPD or other lung disease - no", 
            "COPD or other lung disease - yes"
          )
        ),
        tbl
      )
    },
    chronic_headache = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(x == 0, "chronic headache - no", "chronic headache - yes")
        ),
        tbl
      )
    },
    fibromyalgia = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(x == 0, "fibromyalgia - no", "fibromyalgia - yes before test")
        ),
        tbl
      )
    },
    obesity_no = {
      tbl <- bind_cols(
        tibble(subgroup = "obesity - no"),
        tbl
      )
    },
    obesity_unknown = {
      tbl <- bind_cols(
        tibble(subgroup = "obesity - unknown"),
        tbl
      )
    },
    obesity_yes = {
      tbl <- bind_cols(
        tibble(subgroup = "obesity - yes"),
        tbl
      )
    },
    high_bmi_no = {
      tbl <- bind_cols(
        tibble(subgroup = "high BMI - no"),
        tbl
      )
    },
    high_bmi_unknown = {
      tbl <- bind_cols(
        tibble(subgroup = "high BMI - unknown"),
        tbl
      )
    },
    high_bmi_yes = {
      tbl <- bind_cols(
        tibble(subgroup = "high BMI - yes"),
        tbl
      )
    },
    chronic_fatigue_syndrome = {
      tbl <- bind_cols(
        tibble(
          subgroup = ifelse(
            x == 0, 
            "chronic fatigue syndrome - no", 
            "chronic fatigue syndrome - yes before test")
        ),
        tbl
      )
    },
    anxiety = {
      tbl <- bind_cols(
        tibble(subgroup = ifelse(x == 0, "anxiety - no", "anxiety - yes before test")),
        tbl
      )
    },
    depression = {
      tbl <- bind_cols(
        tibble(subgroup = ifelse(x == 0, "depression - no", "depression - yes before test")),
        tbl
      )
    },
    ptsd = {
      tbl <- bind_cols(
        tibble(subgroup = ifelse(x == 0, "ptsd - no", "ptsd - yes before test")),
        tbl
      )
    },
    custom = {
      tbl <- bind_cols(
        tibble("{y[2]}" := y[3]),
        tbl
      )
    }
  )
  return(tbl)
}