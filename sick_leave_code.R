# ===   setup   ================================================================
### set working path
setwd("path/to/data")
### load packages
library(glue)
library(grf)
library(lmtest)
library(sandwich)
library(tidyverse)
library(zoo)
### source helper functions
source("CausalForestAnalysisWrapper.R")
source("CausalForestCATETable.R")
source("CausalForestDynamicSubgroups.R")
source("CovariateBalance.R")
source("DiscreteCovariateNames.R")
source("DiscreteCovariatesToOneHot.R")
source("RATEOmnibusTest.R")
source("RATETest.R")
source("RegressionForestAnalysisWrapper.R")
source("SubgroupATETable.R")
### choice of seed
seed <- 949564166
### set number of threads used for training
num_threads <- 2
### set flag for running sensitivity analyses
pcr_sens_train <- FALSE
tun_sens_train <- FALSE
### set ggplot theme
theme_set(theme_pubr())
### helper function for forestploter
get_scale <- function(plot,
                      width_wanted,
                      height_wanted,
                      unit = "in"){
  h <- convertHeight(sum(plot$heights), unit, TRUE)
  w <- convertWidth(sum(plot$widths), unit, TRUE)
  max(c(w/width_wanted,  h/height_wanted))
}
### Load retrospective dataset
data_retrospective_baseline_final_9m <- 
  readRDS("retrospective dataset") %>%
  as_tibble()
### remove data about health after test
data_retrospective_baseline_final_9m <- 
  data_retrospective_baseline_final_9m |>
  mutate(
    fibromyalgia = factor(case_when(
      fibromyalgia == "yes_after_test" ~ "no",
      fibromyalgia == "yes_before_and_after_test" ~ "yes_before_test",
      TRUE ~ as.character(fibromyalgia)
    )),
    chronic_fatigue_syndrome = factor(case_when(
      chronic_fatigue_syndrome == "yes_after_test" ~ "no",
      chronic_fatigue_syndrome == "yes_before_and_after_test" ~ "yes_before_test",
      TRUE ~ as.character(chronic_fatigue_syndrome)
    )),
    anxiety = factor(case_when(
      anxiety == "yes_after_test" ~ "no",
      anxiety == "yes_before_and_after_test" ~ "yes_before_test",
      TRUE ~ as.character(anxiety)
    )),
    depression = factor(case_when(
      depression == "yes_after_test" ~ "no",
      depression == "yes_before_and_after_test" ~ "yes_before_test",
      TRUE ~ as.character(depression)
    )),
    ptsd = factor(case_when(
      ptsd == "yes_after_test" ~ "no",
      ptsd == "yes_before_and_after_test" ~ "yes_before_test",
      TRUE ~ as.character(ptsd)
    )),
    # mode impute missing education status (1 person)
    status_education = replace_na(status_education, "higher_medium")
  ) |>
  rename("high_bmi" = "obesity")

# ===   distribution of full time sick leave   =================================
data_retrospective_baseline_final_9m |> 
  count(sick_leave_full_after_cat_duration_retrospective) |>
  mutate(percent = sprintf("%.2f", 100 * n / sum(n)))

data_retrospective_baseline_final_9m |> 
  filter(sex == "female") |>
  count(sick_leave_full_after_cat_duration_retrospective) |>
  mutate(percent = sprintf("%.2f", 100 * n / sum(n)))

data_retrospective_baseline_final_9m |> 
  filter(sex == "male") |>
  count(sick_leave_full_after_cat_duration_retrospective) |>
  mutate(percent = sprintf("%.2f", 100 * n / sum(n)))

data_retrospective_baseline_final_9m |> 
  filter(age < 50) |>
  count(sick_leave_full_after_cat_duration_retrospective) |>
  mutate(percent = sprintf("%.2f", 100 * n / sum(n)))

data_retrospective_baseline_final_9m |> 
  filter(age >= 50) |>
  count(sick_leave_full_after_cat_duration_retrospective) |>
  mutate(percent = sprintf("%.2f", 100 * n / sum(n)))

data_retrospective_baseline_final_9m |>
  count(test_result, sick_leave_full_after_cat_duration_retrospective) |>
  group_by(test_result) |>
  mutate(percent = sprintf("%.1f", 100 * n / sum(n))) |>
  ungroup()

# ===   test difference between exposure groups   ==============================
# Pearson's chi square test to test hypothesis of independence between exposure
# and discrete participant characteristic. Student's t-test to test hypothesis
# of equal group means for continuous participant characteristics within each
# exposure group. 
chi_test_sex <- chisq.test(matrix(c(34085, 17251, 23002, 14480), ncol = 2))
chi_test_charlson <- 
  chisq.test(
    matrix(c(45984, 2854, 1926, 572, 33988, 1967, 1144, 383), ncol = 2)
  )
chi_test_vaccine <- 
  chisq.test(
    matrix(c(51065, 270, *, 37366, *, *), ncol = 2) 
  )
chi_test_bmi <- 
  chisq.test(
    matrix(c(38696, 4143, 8497, 28599, 2680, 6203), ncol = 2)
  )
chi_test_fibro <- chisq.test(matrix(c(50862, 474, 37153, 329), ncol = 2))
chi_test_fatigue <- chisq.test(matrix(c(50473, 863, 37000, 482), ncol = 2))
chi_test_anxiety <- chisq.test(matrix(c(46905, 4431, 34488, 2994), ncol = 2))
chi_test_depres <- chisq.test(matrix(c(44870, 6466, 33233, 4249), ncol = 2))
chi_test_ptsd <- chisq.test(matrix(c(51299, 1037, 36762, 720), ncol = 2))
chi_test_asthma <- chisq.test(matrix(c(47840, 3496, 34643, 2839), ncol = 2))
chi_test_diabet <- chisq.test(matrix(c(49583, 1753, 36255, 1227), ncol = 2))
chi_test_hbp <- chisq.test(matrix(c(45412, 5924, 33589, 3893), ncol = 2))
chi_test_copd <- chisq.test(matrix(c(50565, 771, 37056, 426), ncol = 2))
chi_test_headac <- chisq.test(matrix(c(49212, 2124, 36013, 1469), ncol = 2))

t_test_age <- t.test(
  age ~ test_result, 
  data = data_retrospective_baseline_final_9m
)

test_pvals <- setNames(
  c(
    t_test_age$p.value, chi_test_sex$p.value, chi_test_charlson$p.value, 
    chi_test_vaccine$p.value, chi_test_bmi$p.value, chi_test_fibro$p.value, 
    chi_test_fatigue$p.value, chi_test_anxiety$p.value, chi_test_depres$p.value, 
    chi_test_ptsd$p.value, chi_test_asthma$p.value, chi_test_diabet$p.value, 
    chi_test_hbp$p.value, chi_test_copd$p.value, chi_test_headac$p.value
  ),
  c(
    "Age", "Sex", "Charlson Comorbidity Index", "Vaccination status", 
    "High BMI", "Fibromyalgia", "Chronic fatigue syndrome", "Anxiety", 
    "Depression", "PTSD", "Chronic asthma", "Diabetes", "High blood pressure", 
    "COPD or other lung disease", "Chronic or frequent headaches or migraines"
  )
)

# ===   train causal forest   ==================================================
set.seed(seed)
forest_Y <- RegressionForestAnalysisWrapper(
  data = data_retrospective_baseline_final_9m,
  cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", 
          "anxiety", "depression", "ptsd"),
  binary = c(1, 3, 4, 5, 6),
  out = "Y", 
  num.trees = 2000, num.threads = num_threads,
  seed = seed,
  alpha = 0.01, min.node.size = 10
)
forest_W <- RegressionForestAnalysisWrapper(
  data = data_retrospective_baseline_final_9m,
  cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", 
          "anxiety", "depression", "ptsd"),
  binary = c(1, 3, 4, 5, 6),
  out = "W", 
  num.trees = 2000, num.threads = num_threads,
  seed = seed,
  alpha = 0.01, min.node.size = 10
)
Y_hat <- predict(forest_Y)$predictions
W_hat <- predict(forest_W)$predictions
cf <- CausalForestAnalysisWrapper(
  data = data_retrospective_baseline_final_9m,
  cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", 
          "anxiety", "depression", "ptsd"),
  binary = c(1, 3, 4, 5, 6),
  Y.hat = Y_hat, W.hat = W_hat,
  num.trees = 2000, num.threads = num_threads,
  seed = seed, 
  alpha = 0.01, min.node.size = 10
)

#  ===   Check forest   ========================================================
### variable importance
variable_importance <- tibble(
  variable_name = names(as_tibble(cf$X.orig)),
  variable_importance = as.numeric(variable_importance(cf))
) |>
  arrange(desc(variable_importance)) |>
  mutate(
    variable_importance_num = variable_importance,
    variable_importance = sprintf("%.3f", variable_importance)
  )

### Covariate balance
X_names_reference <- c(
  "age", "sex", "charlson_score_cat", "chronic_asthma", "chronic_diabetes",
  "chronic_high_blood_pressure", "chronic_copd_lung_disease", 
  "chronic_headache", "fibromyalgia", "chronic_fatigue_syndrome", "anxiety",
  "depression", "ptsd", "status_education", "high_bmi"
)
names(X_names_reference) <- c(
  "age", "sex", "Charlson Comorbidity Index", "chronic asthma", "diabetes", 
  "high blood pressure", "COPD or other chronic lung disease", 
  "chronic or frequent headaches/migraines", "fibromyalgia",                       
  "chronic fatigue syndrome", "anxiety", "depression", 
  "post-traumatic stress disorder",  "education", "high BMI"
)
balance <- CovariateBalance(
  cf, 
  plots = "all",
  balance_table = TRUE,
  names = X_names_reference,
  factor = list(
    sex = c("female" = 0, "male" = 1),
    `Charlson Comorbidity Index` = c("0" = 0, "1" = 1, "2" = 2, "3 or more" = 3),
    `chronic asthma` = c("no" = 0, "yes" = 1),
    diabetes = c("no" = 0, "yes" = 1),
    `high blood pressure` = c("no" = 0, "yes" = 1),
    `COPD or other chronic lung disease` = c("no" = 0, "yes" = 1),
    `chronic or frequent headaches/migraines` = c("no" = 0, "yes" = 1),
    `fibromyalgia` = c("no" = 0, "yes" = 1),
    `chronic fatigue syndrome` = c("no" = 0, "yes" = 1),
    `anxiety` = c("no" = 0, "yes" = 1),
    `depression` = c("no" = 0, "yes" = 1),
    `post-traumatic stress disorder` = c("no" = 0, "yes" = 1),
    education = c(
      "primary" = "primary",
      "secondary" = "secondary",
      "vocational training" = "vocational_training",
      "higher 1-2 years" = "higher_short",
      "higher  2-4 years" = "higher_medium",
      "higher > 4 years" = "higher_long",
      "unknown" = "do_not_know"
    )
  ),
  treatment_name = "RT-PCR test result:",
  cd_ncol = 3,
  ec_ncol = 3,
  love_breaks = c(0, 0.1),
  love_xlim = c(0, 0.12)
)
X_names_reference_manu <- c(
  "age", "sex", "charlson_score_cat", "chronic_asthma", "chronic_diabetes",
  "chronic_high_blood_pressure", "chronic_copd_lung_disease", 
  "chronic_headache", "fibromyalgia", "chronic_fatigue_syndrome", "anxiety",
  "depression", "ptsd", "status_education", "high_bmi"
)
names(X_names_reference_manu) <- c(
  "age", "sex", "Charlson Comorbidity Index", "chronic asthma", "diabetes", 
  "high blood pressure", "COPD", 
  "chronic headaches/migraines", "fibromyalgia",                       
  "chronic fatigue syndrome", "anxiety", "depression", 
  "ptsd",  "education", "high BMI"
)
balance_manu_plot <- CovariateBalance(
  cf, 
  sick_leave = TRUE,
  plots = "ecdf",
  balance_table = FALSE,
  names = X_names_reference_manu,
  factor = list(
    sex = c("female" = 0, "male" = 1),
    `Charlson Comorbidity Index` = c("0" = 0, "1" = 1, "2" = 2, "3+" = 3),
    `chronic asthma` = c("no" = 0, "yes" = 1),
    diabetes = c("no" = 0, "yes" = 1),
    `high blood pressure` = c("no" = 0, "yes" = 1),
    `COPD` = c("no" = 0, "yes" = 1),
    `chronic headaches/migraines` = c("no" = 0, "yes" = 1),
    `fibromyalgia` = c("no" = 0, "yes" = 1),
    `chronic fatigue syndrome` = c("no" = 0, "yes" = 1),
    `anxiety` = c("no" = 0, "yes" = 1),
    `depression` = c("no" = 0, "yes" = 1),
    `ptsd` = c("no" = 0, "yes" = 1),
    education = c(
      "a" = "primary",
      "b" = "secondary",
      "c" = "vocational_training",
      "d" = "higher_short",
      "e" = "higher_medium",
      "f" = "higher_long",
      "g" = "do_not_know"
    )
  ),
  treatment_name = "RT-PCR test result:"
)

### c-for-benefit statistic
# matching on CATE with random tie breaking and replacement
matched_on_cate <-
  Matching::Match(
    Tr = cf$W.orig, 
    X = cf$predictions, 
    M = 1, 
    ties = FALSE, 
    replace = TRUE
  )
ind.A <- matched_on_cate$index.control
ind.B <- matched_on_cate$index.treatedd
# calculate average predicted benefit in matched pairs
pred.harm.A <- cf$predictions[ind.A]
pred.harm.B <- cf$predictions[ind.B]
pred.harm.avg <- (pred.harm.A + pred.harm.B) / 2
# calculate observed benefit in matched pairs
obs.out.A <- cf$Y.orig[ind.A]
obs.out.B <- cf$Y.orig[ind.B]
obs.harm <- obs.out.B - obs.out.A
# Benefit c-statistic
cindex <- Hmisc::rcorr.cens(pred.harm.avg, obs.harm)
c.harm <- cindex["C Index"][[1]]
c.harm.se <- cindex["S.D."][[1]] / 2
c_for_benefit <- tibble(
  c_for_harm = c.harm,
  c_for_harm_lower = c.harm - 1.96 * c.harm.se,
  c_for_harm_upper = c.harm + 1.96 * c.harm.se
)

# ===   Test heterogeneity   ===================================================
### Omnibus test based on linear model using demeaned covariates and outcomes
test_calibration <- test_calibration(cf) |>
  structure(class = NULL) |>
  as_tibble() |>
  (\(x) bind_cols(
    tibble(` ` = c("mean forest prediction", "differential forest prediction")), 
    x
  ))()

### Test based on linear approximation of CATE estimates
blp <- best_linear_projection(cf, cf$X.orig)

### RATE based omnibus test
rate_omnibus_test <- RATEOmnibusTest(cf, alpha = 0.01, min.node.size = 10)

### RATE based tests along covariates
# age prioritization score
rate_age <- RATETest(cf, "age")

# sex prioritization score
rate_sex <- RATETest(cf, "sex", "discrete")

# education prioritization score
education_scores <- tibble(
  education = case_when(
    cf$X.orig$status_education_primary == 1 ~ "primary",
    cf$X.orig$status_education_secondary == 1 ~ "secondary",
    cf$X.orig$status_education_vocational_training == 1 ~ "vocational_training",
    cf$X.orig$status_education_higher_short == 1 ~ "higher_short",
    cf$X.orig$status_education_higher_medium == 1 ~ "higher_medium",
    cf$X.orig$status_education_higher_long == 1 ~ "higher_long",
    cf$X.orig$status_education_do_not_know == 1 ~ "do_not_know",
  ),
  dr_scores = get_scores(cf),
) |>
  mutate(
    priority = case_when(
      education == "primary" ~ 6,
      education == "secondary" ~ 5,
      education == "vocational_training" ~ 4,
      education == "higher_short" ~ 3,
      education == "higher_medium" ~ 2,
      education == "higher_long" ~ 1,
      TRUE ~ NA_real_
    )
  ) |>
  filter(education != "do_not_know")
education_q <- c(
  0.001,
  cumsum(rev(table(education_scores[["priority"]]))) / 
    length(education_scores[["priority"]])
)
rate_education <- list(
  rate = rank_average_treatment_effect.fit(
    DR.scores = education_scores[["dr_scores"]],
    priorities = education_scores[["priority"]],
    target = "AUTOC",
    q = education_q,
    R = 500
  )) %>%
  c(list(
    confint = .$rate$estimate + 
      tibble(
        estimate = 0,
        lower = -qnorm(0.975) * .$rate$std.err,
        upper =  qnorm(0.975) * .$rate$std.err
      ),
    pval = 2 * pnorm(-abs(.$rate$estimate) / .$rate$std.err)
  ))

# Charlson score prioritization score
rate_charlson_score <- 
  RATETest(cf, "charlson_score_cat", "discrete")

# chronic asthma prioritization score
rate_chronic_asthma <- 
  RATETest(cf, "chronic_asthma", "discrete")

# chronic diabetes prioritization score
rate_chronic_diabetes <- 
  RATETest(cf, "chronic_diabetes", "discrete")

# chronic high blood pressure prioritization score
rate_chronic_high_blood_pressure <- 
  RATETest(cf, "chronic_high_blood_pressure", "discrete")

# chronic copd lung disease prioritization score
rate_chronic_copd_lung_disease <- 
  RATETest(cf, "chronic_copd_lung_disease", "discrete")

# chronic headache prioritization score
rate_chronic_headache <- 
  RATETest(cf, "chronic_headache", "discrete")

# fibromyalgia
rate_fibromyalgia <- 
  RATETest(cf, "fibromyalgia", "discrete")

# high_bmi (yes/no)
high_bmi_scores <- tibble(
  high_bmi = case_when(
    cf$X.orig$high_bmi_yes == 1 ~ "yes",
    cf$X.orig$high_bmi_no == 1 ~ "no",
    cf$X.orig$high_bmi_unknown == 1 ~ "unknown",
  ),
  dr_scores = get_scores(cf),
) |>
  mutate(
    priority = case_when(
      high_bmi == "yes" ~ 1,
      high_bmi == "no" ~ 0,
      TRUE ~ NA_real_
    )
  ) |>
  filter(high_bmi != "unknown")
high_bmi_q <- c(
  0.001,
  cumsum(rev(table(high_bmi_scores[["priority"]]))) / 
    length(high_bmi_scores[["priority"]])
)
rate_high_bmi <- list(
  rate = rank_average_treatment_effect.fit(
    DR.scores = high_bmi_scores[["dr_scores"]],
    priorities = high_bmi_scores[["priority"]],
    target = "AUTOC",
    q = high_bmi_q,
    R = 500
  )) %>%
  c(list(
    confint = .$rate$estimate + 
      tibble(
        estimate = 0,
        lower = -qnorm(0.975) * .$rate$std.err,
        upper =  qnorm(0.975) * .$rate$std.err
      ),
    pval = 2 * pnorm(-abs(.$rate$estimate) / .$rate$std.err)
  ))

# chronic fatigue syndrome
rate_chronic_fatigue_syndrome <- 
  RATETest(cf, "chronic_fatigue_syndrome", "discrete")

# anxiety
rate_anxiety <- RATETest(cf, "anxiety", "discrete")

# depression
rate_depression <- RATETest(cf, "depression", "discrete")

# ptsd 
rate_ptsd <- RATETest(cf, "ptsd", "discrete")

## combine RATE estimates
rate_table <- bind_rows(
  as_tibble(
    c(covariate = "omnibus", 
      rate_omnibus_test$confint, 
      pval = rate_omnibus_test$pval)
  ),
  as_tibble(
    c(covariate = "age", 
      rate_age$confint, 
      pval = rate_age$pval)
  ),
  as_tibble(
    c(covariate = "sex", 
      rate_sex$confint, 
      pval = rate_sex$pval)
  ),
  as_tibble(
    c(covariate = "education", 
      rate_education$confint, 
      pval = rate_education$pval)
  ),
  as_tibble(
    c(covariate = "charlson_score", 
      rate_charlson_score$confint, 
      pval = rate_charlson_score$pval)
  ),
  as_tibble(
    c(covariate = "chronic_asthma", 
      rate_chronic_asthma$confint, 
      pval = rate_chronic_asthma$pval)
  ),
  as_tibble(
    c(covariate = "chronic_copd_lung_disease", 
      rate_chronic_copd_lung_disease$confint, 
      pval = rate_chronic_copd_lung_disease$pval)
  ),
  as_tibble(
    c(covariate = "chronic_diabetes", 
      rate_chronic_diabetes$confint, 
      pval = rate_chronic_diabetes$pval)
  ),
  as_tibble(
    c(covariate = "chronic_headache", 
      rate_chronic_headache$confint, 
      pval = rate_chronic_headache$pval)
  ),
  as_tibble(
    c(covariate = "chronic_high_blood_pressure", 
      rate_chronic_high_blood_pressure$confint, 
      pval = rate_chronic_high_blood_pressure$pval)
  ),
  as_tibble(
    c(covariate = "fibromyalgia", 
      rate_fibromyalgia$confint, 
      pval = rate_fibromyalgia$pval)
  ),
  as_tibble(
    c(covariate = "high_bmi", 
      rate_high_bmi$confint, 
      pval = rate_high_bmi$pval)
  ),
  as_tibble(
    c(covariate = "chronic_fatigue_syndrome", 
      rate_chronic_fatigue_syndrome$confint, 
      pval = rate_chronic_fatigue_syndrome$pval)
  ),
  as_tibble(
    c(covariate = "anxiety", 
      rate_anxiety$confint, 
      pval = rate_anxiety$pval)
  ),
  as_tibble(
    c(covariate = "depression", 
      rate_depression$confint, 
      pval = rate_depression$pval)
  ),
  as_tibble(
    c(covariate = "ptsd", 
      rate_ptsd$confint, 
      pval = rate_ptsd$pval)
  )
) |>
  transmute(
    covariate = covariate,
    `RATE (AUTOC)` = sprintf("%.2E", estimate),
    `95% CI - lower` = sprintf("%.2E", lower),
    `95% CI - upper` = sprintf("%.2E", upper),
    `p-value` = sprintf("%.2E", pval)
  )

# ===   Average Treatment Effect   =============================================
### ATE in full population
ate_all <- SubgroupATETable(NULL, NULL, cf, 0.95)
### ATE in subgroups along covariates
# age
ate_age <- c(map(seq(10, 60, 10), ~ .x + 0:9), list(14:49, 50:64)) |>
  map_dfr(\ (x) SubgroupATETable(x, "age", cf, 0.95))
# sex
ate_sex <- 0:1 |>
  map_dfr(\ (x) SubgroupATETable(x, "sex", cf, 0.95))
# education
ate_education <- stringr::str_subset(names(cf$X.orig), "^status_education") |>
  map_dfr(
    \ (y) {
      SubgroupATETable(
        x = NULL,
        y = c(
          "subgroup", 
          paste0(
            "education - ", 
            stringr::str_replace_all(
              stringr::str_remove(y, "status_education_"), "_", " "
            )
          )
        ), 
        cf = cf, 
        level = 0.95, 
        subset = cf[["X.orig"]][[y]] == 1
      )
    }
  )
# Charlson score
ate_charlson_score <- 0:3 |>
  map_dfr(\(x) SubgroupATETable(x, "charlson_score_cat", cf, 0.95))
# chronic_asthma
ate_chronic_asthma <- 0:1 |>
  map_dfr(\(x) SubgroupATETable(x, "chronic_asthma", cf, 0.95))
# chronic_diabetes
ate_chronic_diabetes <- 0:1 |>
  map_dfr(\(x) SubgroupATETable(x, "chronic_diabetes", cf, 0.95))
# chronic_high_blood_pressure
ate_chronic_high_blood_pressure <- 0:1 |>
  map_dfr(\(x) SubgroupATETable(x, "chronic_high_blood_pressure", cf, 0.95))
# chronic_copd_lung_disease
ate_chronic_copd_lung_disease <- 0:1 |>
  map_dfr(\(x) SubgroupATETable(x, "chronic_copd_lung_disease", cf, 0.95))
# chronic_headache
ate_chronic_headache <- 0:1 |>
  map_dfr(\(x) SubgroupATETable(x, "chronic_headache", cf, 0.95))
# fibromyalgia
ate_fibromyalgia <- 0:1 |>
  map_dfr(\(x) SubgroupATETable(x, "fibromyalgia", cf, 0.95))
# high_bmi
ate_high_bmi <- map2_dfr(
  rep(1, 3),
  str_subset(names(cf$X.orig), "high_bmi"),
  \(x, y) SubgroupATETable(x, y, cf, 0.95)
)
# chronic_fatigue_syndrome
ate_chronic_fatigue_syndrome <- 0:1 |>
  map_dfr(\(x) {
    SubgroupATETable(
      x, 
      "chronic_fatigue_syndrome", 
      cf, 
      0.95)
  })
# anxiety
ate_anxiety <- 0:1 |>
  map_dfr(\ (x) SubgroupATETable(x, "anxiety", cf, 0.95))
# depression
ate_depression <- 0:1 |>
  map_dfr(\ (x) SubgroupATETable(x, "depression", cf, 0.95))
# ptsd
ate_ptsd <- 0:1 |>
  map_dfr(\ (x) SubgroupATETable(x, "ptsd", cf, 0.95))
## combine ATE estimates
ate_table <- bind_rows(
  ate_all,
  ate_age,
  ate_sex,
  ate_charlson_score,
  ate_education,
  ate_chronic_asthma,
  ate_chronic_copd_lung_disease,
  ate_chronic_diabetes,
  ate_chronic_headache,
  ate_chronic_high_blood_pressure,
  ate_fibromyalgia,
  ate_high_bmi,
  ate_chronic_fatigue_syndrome,
  ate_anxiety,
  ate_depression,
  ate_ptsd
) |>
  mutate(
    `estimate (%)` = sprintf("%.1f", 100 * estimate),
    `95% CI - lower (%)` = sprintf("%.1f", 100 * `95% CI - lower`),
    `95% CI - upper (%)` = sprintf("%.1f", 100 * `95% CI - upper`)
  )

# ===   CATE estimates   =======================================================
### out-of-bag CATE estimates
oob <- cf$X.orig |>
  bind_cols(
    predict(cf, estimate.variance = TRUE, num.threads = num_threads)
  ) |>
  mutate(
    age_grp_3 = factor((\(age) {
      x <- ceiling((age - 13) / 3)
      y <- 13 + x * 3
      paste0(y - 2, "-", y)
    })(age))
  )

# ===   CATE - age, high bmi, and depression   =================================
### CATE over age, high_bmi, and depression
cate_depression_high_bmi_age <- CausalForestCATETable(
  cf,
  list(
    depression = rep(list("no" = 0, "yes" = 1), each = 3, times = 5),
    high_bmi = rep(as.list(str_subset(names(cf$X.orig), "^high_bmi")), times = 10),
    age = c(rep(list("15-25" = 15:25), 6), 
            rep(list("26-35" = 26:35), 6), 
            rep(list("36-45" = 36:45), 6), 
            rep(list("46-55" = 46:55), 6),
            rep(list("56-65" = 56:65), 6))
  )
)

### data for plots with age, high_bmi, and depression
cate_histogram_data_depression_high_bmi_age <- oob |>
  filter(high_bmi_unknown == 0) |>
  mutate(
    depression_high_bmi_age = factor(case_when(
      depression == 0 & high_bmi_yes == 1 & age < 50 ~ "<50 yrs, no depression, high BMI",
      depression == 0 & high_bmi_yes == 1 & age >= 50 ~ "\u226550 yrs, no depression, high BMI",
      depression == 0 & high_bmi_no == 1 & age < 50 ~ "<50 yrs, no depression, not high BMI",
      depression == 0 & high_bmi_no == 1 & age >= 50 ~ "\u226550 yrs, no depression, not high BMI",
      depression == 1 & high_bmi_yes == 1 & age < 50 ~ "<50 yrs, depression, high BMI",
      depression == 1 & high_bmi_yes == 1 & age >= 50 ~ "\u226550 yrs, depression, high BMI",
      depression == 1 & high_bmi_no == 1 & age < 50 ~ "<50 yrs, depression, not high BMI",
      depression == 1 & high_bmi_no == 1 & age >= 50 ~ "\u226550 yrs, depression, not high BMI",
      TRUE ~ "unknown"
    ), levels = c(
      "<50 yrs, no depression, not high BMI", "\u226550 yrs, no depression, not high BMI",
      "<50 yrs, no depression, high BMI", "\u226550 yrs, no depression, high BMI",
      "<50 yrs, depression, not high BMI", "\u226550 yrs, depression, not high BMI",
      "<50 yrs, depression, high BMI", "\u226550 yrs, depression, high BMI"
    )),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_depression_high_bmi_age <- 
  cate_histogram_data_depression_high_bmi_age |>
  group_by(depression_high_bmi_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.42,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

depression_high_bmi_age_cate_table_lines_data <-
  cate_depression_high_bmi_age |>
  filter(!str_detect(depression_high_bmi_age, "unknown")) |>
  mutate(
    age_grp = case_when(
      str_detect(depression_high_bmi_age, "15-25") ~ 20,
      str_detect(depression_high_bmi_age, "26-35") ~ 30.5,
      str_detect(depression_high_bmi_age, "36-45") ~ 40.5,
      str_detect(depression_high_bmi_age, "46-55") ~ 50.5,
      str_detect(depression_high_bmi_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    depression_high_bmi = factor(case_when(
      str_detect(depression_high_bmi_age, "no_no") ~ "no_no",
      str_detect(depression_high_bmi_age, "yes_no") ~ "yes_no",
      str_detect(depression_high_bmi_age, "no_yes") ~ "no_yes",
      str_detect(depression_high_bmi_age, "yes_yes") ~ "yes_yes",
      TRUE ~ NA_character_
    ),
    levels = c("no_no", "yes_no", "no_yes", "yes_yes"),
    labels = c("no depression, not high BMI", "depression, not high BMI", 
               "no depression, high BMI", "depression, high BMI")
    )
  )

# additional data for result section
histogram_summary_depression_high_bmi_age <- 
  cate_histogram_data_depression_high_bmi_age |>
  group_by(depression_high_bmi_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - age, high BMI, and sex   ========================================
### CATE over age, high BMI, and sex
cate_sex_high_bmi_age <- CausalForestCATETable(
  cf,
  list(
    sex = rep(list("female" = 0, "male" = 1), each = 3, times = 5),
    high_bmi = rep(as.list(str_subset(names(cf$X.orig), "^high_bmi")), times = 10),
    age = c(rep(list("15-25" = 15:25), 6), 
            rep(list("26-35" = 26:35), 6), 
            rep(list("36-45" = 36:45), 6), 
            rep(list("46-55" = 46:55), 6),
            rep(list("56-65" = 56:65), 6))
  )
)

### data for plots with age, high_bmi, and sex
cate_histogram_data_sex_high_bmi_age <- oob |>
  filter(high_bmi_unknown == 0) |>
  mutate(
    sex_high_bmi_age = factor(case_when(
      sex == 0 & high_bmi_yes == 1 & age < 50 ~ "female, <50 yrs, high BMI",
      sex == 0 & high_bmi_yes == 1 & age >= 50 ~ "female, \u226550 yrs, high BMI",
      sex == 0 & high_bmi_no == 1 & age < 50 ~ "female, <50 yrs, not high BMI",
      sex == 0 & high_bmi_no == 1 & age >= 50 ~ "female, \u226550 yrs, not high BMI",
      sex == 1 & high_bmi_yes == 1 & age < 50 ~ "male, <50 yrs, high BMI",
      sex == 1 & high_bmi_yes == 1 & age >= 50 ~ "male, \u226550 yrs, high BMI",
      sex == 1 & high_bmi_no == 1 & age < 50 ~ "male, <50 yrs, not high BMI",
      sex == 1 & high_bmi_no == 1 & age >= 50 ~ "male, \u226550 yrs, not high BMI",
      TRUE ~ "unknown"
    ), levels = c(
      "male, <50 yrs, not high BMI", "female, <50 yrs, not high BMI", 
      "male, <50 yrs, high BMI", "female, <50 yrs, high BMI",  
      "male, \u226550 yrs, not high BMI", "female, \u226550 yrs, not high BMI", 
      "male, \u226550 yrs, high BMI", "female, \u226550 yrs, high BMI"
    )),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_sex_high_bmi_age <- 
  cate_histogram_data_sex_high_bmi_age |>
  group_by(sex_high_bmi_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.5,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

sex_high_bmi_age_cate_table_lines_data <-
  cate_sex_high_bmi_age |>
  filter(!str_detect(sex_high_bmi_age, "unknown")) |>
  mutate(
    age_grp = case_when(
      str_detect(sex_high_bmi_age, "15-25") ~ 20,
      str_detect(sex_high_bmi_age, "26-35") ~ 30.5,
      str_detect(sex_high_bmi_age, "36-45") ~ 40.5,
      str_detect(sex_high_bmi_age, "46-55") ~ 50.5,
      str_detect(sex_high_bmi_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    sex_high_bmi = factor(case_when(
      str_detect(sex_high_bmi_age, "female_no") ~ "female_no",
      str_detect(sex_high_bmi_age, "female_yes") ~ "female_yes",
      str_detect(sex_high_bmi_age, "male_no") ~ "male_no",
      str_detect(sex_high_bmi_age, "male_yes") ~ "male_yes",
      TRUE ~ NA_character_
    ),
    levels = c("male_no", "male_yes", "female_no", "female_yes"),
    labels = c("male, not high BMI", "male, high BMI", 
               "female, not high BMI", "female, high BMI")
    )
  )

# additional data for result section
histogram_summary_sex_high_bmi_age <- 
  cate_histogram_data_sex_high_bmi_age |>
  group_by(sex_high_bmi_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - age, depression, and sex   ======================================
### CATE over age, depression, and sex
cate_sex_depression_age <- CausalForestCATETable(
  cf,
  list(
    sex = rep(list("female" = 0, "male" = 1), each = 2, times = 5),
    depression = rep(list("no" = 0, "yes" = 1), times = 10),
    age = c(rep(list("15-25" = 15:25), 4), 
            rep(list("26-35" = 26:35), 4), 
            rep(list("36-45" = 36:45), 4), 
            rep(list("46-55" = 46:55), 4),
            rep(list("56-65" = 56:65), 4))
  )
)

### data for plots with age, depression, and sex
cate_histogram_data_sex_depression_age <- oob |>
  mutate(
    sex_depression_age = factor(case_when(
      sex == 0 & depression == 1 & age < 50 ~ "female, <50 yrs, depression",
      sex == 0 & depression == 1 & age >= 50 ~ "female, \u226550 yrs, depression",
      sex == 0 & depression == 0 & age < 50 ~ "female, <50 yrs, no depression",
      sex == 0 & depression == 0 & age >= 50 ~ "female, \u226550 yrs, no depression",
      sex == 1 & depression == 1 & age < 50 ~ "male, <50 yrs, depression",
      sex == 1 & depression == 1 & age >= 50 ~ "male, \u226550 yrs, depression",
      sex == 1 & depression == 0 & age < 50 ~ "male, <50 yrs, no depression",
      sex == 1 & depression == 0 & age >= 50 ~ "male, \u226550 yrs, no depression",
      TRUE ~ "unknown"
    ), levels = c(
      "male, <50 yrs, no depression", "female, <50 yrs, no depression", 
      "male, <50 yrs, depression", "female, <50 yrs, depression",
      "male, \u226550 yrs, no depression", "female, \u226550 yrs, no depression", 
      "male, \u226550 yrs, depression", "female, \u226550 yrs, depression"
    )),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_sex_depression_age <- 
  cate_histogram_data_sex_depression_age |>
  group_by(sex_depression_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.455,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

sex_depression_age_cate_table_lines_data <-
  cate_sex_depression_age |>
  mutate(
    age_grp = case_when(
      str_detect(sex_depression_age, "15-25") ~ 20,
      str_detect(sex_depression_age, "26-35") ~ 30.5,
      str_detect(sex_depression_age, "36-45") ~ 40.5,
      str_detect(sex_depression_age, "46-55") ~ 50.5,
      str_detect(sex_depression_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    sex_depression = factor(case_when(
      str_detect(sex_depression_age, "female_no") ~ "female_no",
      str_detect(sex_depression_age, "female_yes") ~ "female_yes",
      str_detect(sex_depression_age, "male_no") ~ "male_no",
      str_detect(sex_depression_age, "male_yes") ~ "male_yes",
      TRUE ~ NA_character_
    ),
    levels = c("male_no", "male_yes", "female_no", "female_yes"),
    labels = c("male, no depression", "male, depression", 
               "female, no depression", "female, depression")
    )
  )

# additional data for result section
histogram_summary_sex_depression_age <- 
  cate_histogram_data_sex_depression_age |>
  group_by(sex_depression_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - high BMI, depression, and sex   =================================
### table with effects across high BMI, depression, and sex
cate_high_bmi_depression_sex <- CausalForestCATETable(
  cf,
  list(
    high_bmi = rep(as.list(str_subset(names(cf$X.orig), "^high_bmi")), each = 4),
    depression = rep(list("no" = 0, "yes" = 1), each = 2, times = 3),
    sex = rep(list("female" = 0, "male" = 1), times = 6)
  )
)

### t.test of difference in means of different groups
aipw_scores <- policytree::double_robust_scores(cf)
aipw_scores <- aipw_scores[,2] - aipw_scores[,1]
# t-test from linear regression for each combination of BMI, sex, and depression
data <- bind_cols(oob, tibble(aipw = aipw_scores)) |>
  mutate(
    high_bmi_depression_sex = paste0(
      case_when(
        high_bmi_no == 1 ~ "not_high_bmi_",
        high_bmi_unknown == 1 ~ "unknown_bmi_",
        high_bmi_yes == 1 ~ "high_bmi_"
      ),
      ifelse(depression == 0, "no_depression_", "depression_"),
      ifelse(sex == 0, "female_", "male_")
    )
  ) |>
  select(high_bmi_depression_sex, aipw)
res <- tibble()
tmp <- data |> 
  mutate(high_bmi_depression_sex = factor(high_bmi_depression_sex))
for(i in seq_along(unique(data$high_bmi_depression_sex))) {
  string <- unique(data$high_bmi_depression_sex)[i]
  tmp <- tmp |> 
    mutate(high_bmi_depression_sex = fct_relevel(high_bmi_depression_sex, string))
  ols <- lm(aipw ~ 1 + high_bmi_depression_sex, data = tmp)
  res <- rbind(
    res, 
    as_tibble(
      coef(summary(ols))[-(1:i), c(1,2,4), drop=FALSE]
    ) %>%
      mutate(
        id = paste(
          str_remove(string, "_$"), 
          "vs.", 
          str_remove(
            str_remove(
              dimnames(summary(ols)$coefficients)[[1]][-(1:i)], 
              "^high_bmi_depression_sex"
            ), 
            "_$"
          )
        )
      )  
  )
}
# Adjust p-values with Benjamini-Hockberg:
cate_high_bmi_depression_sex_ttest <- res %>%
  rename(`Orig. p-value` = `Pr(>|t|)`) %>%
  mutate(
    `95 % CI` = str_c("(",
                      round(Estimate - qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ", ",
                      round(Estimate + qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ")"),
    `Orig. p-value` = round(
      `Orig. p-value`,
      digits = 4),
    `Adj. p-value` = round(
      p.adjust(`Orig. p-value`, method = "BH"),
      digits = 4)) %>%
  select(id, Estimate, `Std. Error`, `95 % CI`, everything())

# two sample t-test for high BMI, no depression vs. not high BMI, depression for males
cate_high_bmi_depression_sex_two_samp_ttest <- t.test(
  data |>
    filter(high_bmi_depression_sex == "not_high_bmi_depression_male_") |>
    pull(aipw),
  data |> 
    filter(high_bmi_depression_sex == "high_bmi_no_depression_male_") |>
    pull(aipw)
)

### data for plots with high BMI, depression, and sex
cate_histogram_data_high_bmi_depression_sex <- oob |>
  filter(high_bmi_unknown == 0) |>
  mutate(
    high_bmi_depression_sex = factor(case_when(
      sex == 0 & depression == 0 & high_bmi_no == 1 ~ "female, no depression, not high BMI",
      sex == 0 & depression == 1 & high_bmi_no == 1 ~ "female, depression, not high BMI",
      sex == 0 & depression == 0 & high_bmi_yes == 1 ~ "female, no depression, high BMI",
      sex == 0 & depression == 1 & high_bmi_yes == 1 ~ "female, depression, high BMI",
      depression == 0 & high_bmi_no == 1 ~ "male, no depression, not high BMI",
      depression == 1 & high_bmi_no == 1 ~ "male, depression, not high BMI",
      depression == 0 ~ "male, no depression, high BMI",
      TRUE ~ "male, depression, high BMI"
    ),
    levels = c(
      "female, no depression, not high BMI", "male, no depression, not high BMI", 
      "female, depression, not high BMI", "male, depression, not high BMI", 
      "female, no depression, high BMI", "male, no depression, high BMI", 
      "female, depression, high BMI", "male, depression, high BMI"
    )),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )

cate_histogram_label_high_bmi_depression_sex <- 
  cate_histogram_data_high_bmi_depression_sex |>
  group_by(high_bmi_depression_sex) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.41,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

high_bmi_depression_sex_cate_table_data <- 
  cate_high_bmi_depression_sex |>
  filter(!str_detect(high_bmi_depression_sex, "unknown")) |>
  mutate(high_bmi_depression_sex = factor(
    high_bmi_depression_sex,
    levels = c(
      "no_no_female", "no_yes_female",
      "yes_no_female", "yes_yes_female",
      "no_no_male", "no_yes_male",
      "yes_no_male", "yes_yes_male"
    ),
    labels = c(
      "female, no depression, not high BMI", "female, depression, not high BMI",
      "female, no depression, high BMI", "female, depression, high BMI",
      "male, no depression, not high BMI", "male, depression, not high BMI",
      "male, no depression, high BMI", "male, depression, high BMI"
    )
  ))

# additional data for result section
histogram_summary_high_bmi_depression_sex <- 
  cate_histogram_data_high_bmi_depression_sex |>
  group_by(high_bmi_depression_sex) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - age, sex, and chronic asthma   ==================================
### table with effects across age, sex, and chronic asthma
cate_sex_chronic_asthma_age <- CausalForestCATETable(
  cf,
  list(
    sex = rep(list("female" = 0, "male" = 1), each = 2, times = 5),
    chronic_asthma = rep(list("asthma" = 1, "no_asthma" = 0), times = 10),
    age = c(rep(list("15-25" = 15:25), 4), 
            rep(list("26-35" = 26:35), 4), 
            rep(list("36-45" = 36:45), 4), 
            rep(list("46-55" = 46:55), 4),
            rep(list("56-65" = 56:65), 4))
  )
)

### t-tests for difference between asthma and no asthma
data <- bind_cols(oob, tibble(aipw = aipw_scores)) |>
  mutate(
    sex_asthma_age = paste0(
      ifelse(sex == 0, "female_", "male_"),
      ifelse(chronic_asthma == 0, "no_asthma_", "asthma_"),
      case_when(
        age %in% 15:25 ~ "15-25",
        age %in% 26:35 ~ "26-35",
        age %in% 36:45 ~ "36-45",
        age %in% 46:55 ~ "46-55",
        age %in% 56:65 ~ "56-65"
      )
    )
  ) |>
  select(sex_asthma_age, aipw)
res <- tibble()
for(i in str_subset(unique(data$sex_asthma_age), "no")) {
  tmp <- data |> 
    filter(
      str_detect(sex_asthma_age, paste0("^",i)) | 
        str_detect(sex_asthma_age, paste0("^", str_remove(i, "_no")))
    )
  ols <- lm(aipw ~ 1 + factor(sex_asthma_age), data = tmp)
  res <- rbind(
    res, 
    as_tibble(
      coef(summary(ols))[2, c(1,2,4), drop=FALSE]
    ) %>%
      mutate(id = paste0(unique(tmp$sex_asthma_age), collapse = " vs. "))  
  )
}
# Adjust p-values with Benjamini-Hockberg:
cate_sex_chronic_asthma_age_ttest <- res %>%
  rename(`Orig. p-value` = `Pr(>|t|)`) %>%
  mutate(
    `95 % CI` = str_c("(",
                      round(Estimate - qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ", ",
                      round(Estimate + qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ")"),
    `Orig. p-value` = round(
      `Orig. p-value`,
      digits = 4),
    `Adj. p-value` = round(
      p.adjust(`Orig. p-value`, method = "BH"),
      digits = 4)) %>%
  select(id, Estimate, `Std. Error`, `95 % CI`, everything())

### data for plots with sex, chronic asthma and age
cate_histogram_data_sex_asthma_age <- oob |>
  mutate(
    sex_asthma_age = factor(case_when(
      sex == 0 & chronic_asthma == 1 & age < 50 ~ "female, <50 yrs, chronic asthma",
      sex == 0 & chronic_asthma == 1 & age >= 50 ~ "female, \u2265 50yrs, chronic asthma",
      sex == 0 & chronic_asthma == 0 & age < 50 ~ "female, <50 yrs, no chronic asthma",
      sex == 0 & chronic_asthma == 0 & age >= 50 ~ "female, \u226550 yrs, no chronic asthma",
      sex == 1 & chronic_asthma == 1 & age < 50 ~ "male, < 50yrs, chronic asthma",
      sex == 1 & chronic_asthma == 1 & age >= 50 ~ "male, \u226550 yrs, chronic asthma",
      sex == 1 & chronic_asthma == 0 & age < 50 ~ "male, <50 yrs, no chronic asthma",
      sex == 1 & chronic_asthma == 0 & age >= 50 ~ "male, \u226550 yrs, no chronic asthma",
      TRUE ~ "unknown"
    ),
    c(
      "male, <50 yrs, no chronic asthma", "female, <50 yrs, no chronic asthma", 
      "male, <50 yrs, chronic asthma", "female, <50 yrs, chronic asthma", 
      "male, \u226550 yrs, no chronic asthma", "female, \u226550 yrs, no chronic asthma", 
      "male, \u226550 yrs, chronic asthma", "female, \u226550 yrs, chronic asthma"
    )
    ),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_sex_asthma_age <- 
  cate_histogram_data_sex_asthma_age |>
  group_by(sex_asthma_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.59,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

sex_chronic_asthma_age_cate_table_lines_data <- 
  cate_sex_chronic_asthma_age |>
  mutate(
    age_grp = case_when(
      str_detect(sex_chronic_asthma_age, "15-25") ~ 20,
      str_detect(sex_chronic_asthma_age, "26-35") ~ 30.5,
      str_detect(sex_chronic_asthma_age, "36-45") ~ 40.5,
      str_detect(sex_chronic_asthma_age, "46-55") ~ 50.5,
      str_detect(sex_chronic_asthma_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    sex_chronic_asthma = factor(case_when(
      str_detect(sex_chronic_asthma_age, "^male_no_asthma") ~ "male_no_asthma",
      str_detect(sex_chronic_asthma_age, "^male_asthma") ~ "male_asthma",
      str_detect(sex_chronic_asthma_age, "^female_no_asthma") ~ "female_no_asthma",
      str_detect(sex_chronic_asthma_age, "^female_asthma") ~ "female_asthma",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_asthma", "male_asthma", "female_no_asthma", "female_asthma"),
    labels = c("male without asthma", "male with asthma", 
               "female without asthma", "female with asthma")
    ),
    patchwork = factor(case_when(
      str_detect(sex_chronic_asthma_age, "^male_no_asthma") ~ "male_no_asthma",
      str_detect(sex_chronic_asthma_age, "^male_asthma") ~ "male_asthma",
      str_detect(sex_chronic_asthma_age, "^female_no_asthma") ~ "female_no_asthma",
      str_detect(sex_chronic_asthma_age, "^female_asthma") ~ "female_asthma",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_asthma", "male_asthma", "female_no_asthma", "female_asthma"),
    labels = c("male without health condition", "male with health condition", 
               "female without health condtion", "female with health condition")
    )
  )

### additional data for result section
histogram_summary_sex_asthma_age <- 
  cate_histogram_data_sex_asthma_age |>
  group_by(sex_asthma_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - age, sex, and headaches   =======================================
### table with effects across age, sex, and headaches
cate_sex_headaches_age <- CausalForestCATETable(
  cf,
  list(
    sex = rep(list("female" = 0, "male" = 1), each = 2, times = 5),
    chronic_headache = rep(list("headaches" = 1, "no_headaches" = 0), times = 10),
    age = c(rep(list("15-25" = 15:25), 4), 
            rep(list("26-35" = 26:35), 4), 
            rep(list("36-45" = 36:45), 4), 
            rep(list("46-55" = 46:55), 4),
            rep(list("56-65" = 56:65), 4))
  )
)

### t-tests for difference between headaches and no headaches
data <- bind_cols(oob, tibble(aipw = aipw_scores)) |>
  mutate(
    sex_headaches_age = paste0(
      ifelse(sex == 0, "female_", "male_"),
      ifelse(chronic_headache == 0, "no_headaches_", "headaches_"),
      case_when(
        age %in% 15:25 ~ "15-25",
        age %in% 26:35 ~ "26-35",
        age %in% 36:45 ~ "36-45",
        age %in% 46:55 ~ "46-55",
        age %in% 56:65 ~ "56-65"
      )
    )
  ) |>
  select(sex_headaches_age, aipw)
res <- tibble()
for(i in str_subset(unique(data$sex_headaches_age), "no")) {
  tmp <- data |> 
    filter(
      str_detect(sex_headaches_age, paste0("^",i)) | 
        str_detect(sex_headaches_age, paste0("^", str_remove(i, "_no")))
    )
  ols <- lm(aipw ~ 1 + factor(sex_headaches_age), data = tmp)
  res <- rbind(
    res, 
    as_tibble(
      coef(summary(ols))[2, c(1,2,4), drop=FALSE]
    ) %>%
      mutate(id = paste0(unique(tmp$sex_headaches_age), collapse = " vs. "))  
  )
}
# Adjust p-values with Benjamini-Hockberg:
cate_sex_headaches_age_ttest <- res %>%
  rename(`Orig. p-value` = `Pr(>|t|)`) %>%
  mutate(
    `95 % CI` = str_c("(",
                      round(Estimate - qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ", ",
                      round(Estimate + qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ")"),
    `Orig. p-value` = round(
      `Orig. p-value`,
      digits = 4),
    `Adj. p-value` = round(
      p.adjust(`Orig. p-value`, method = "BH"),
      digits = 4)) %>%
  select(id, Estimate, `Std. Error`, `95 % CI`, everything())

### data for plots with sex, headaches and age
cate_histogram_data_sex_headaches_age <- oob |>
  mutate(
    sex_headaches_age = factor(case_when(
      sex == 0 & chronic_headache == 1 & age < 50 ~ "female, <50 yrs, chronic or frequent headaches",
      sex == 0 & chronic_headache == 1 & age >= 50 ~ "female, \u226550 yrs, chronic or frequent headaches",
      sex == 0 & chronic_headache == 0 & age < 50 ~ "female, <50 yrs, no chronic or frequent headaches",
      sex == 0 & chronic_headache == 0 & age >= 50 ~ "female, \u226550 yrs, no chronic or frequent headaches",
      sex == 1 & chronic_headache == 1 & age < 50 ~ "male, <50 yrs, chronic or frequent headaches",
      sex == 1 & chronic_headache == 1 & age >= 50 ~ "male, \u226550 yrs, chronic or frequent headaches",
      sex == 1 & chronic_headache == 0 & age < 50 ~ "male, <50 yrs, no chronic or frequent headaches",
      sex == 1 & chronic_headache == 0 & age >= 50 ~ "male, \u226550 yrs, no chronic or frequent headaches",
      TRUE ~ "unknown"
    ),
    c(
      "male, <50 yrs, no chronic or frequent headaches", "female, <50 yrs, no chronic or frequent headaches", 
      "male, <50 yrs, chronic or frequent headaches", "female, <50 yrs, chronic or frequent headaches", 
      "male, \u226550 yrs, no chronic or frequent headaches", "female, \u226550 yrs, no chronic or frequent headaches", 
      "male, \u226550 yrs, chronic or frequent headaches", "female, \u226550 yrs, chronic or frequent headaches"
    )
    ),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_sex_headaches_age <- 
  cate_histogram_data_sex_headaches_age |>
  group_by(sex_headaches_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.59,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

sex_headaches_age_cate_table_lines_data <- 
  cate_sex_headaches_age |>
  mutate(
    age_grp = case_when(
      str_detect(sex_chronic_headache_age, "15-25") ~ 20,
      str_detect(sex_chronic_headache_age, "26-35") ~ 30.5,
      str_detect(sex_chronic_headache_age, "36-45") ~ 40.5,
      str_detect(sex_chronic_headache_age, "46-55") ~ 50.5,
      str_detect(sex_chronic_headache_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    sex_headaches = factor(case_when(
      str_detect(sex_chronic_headache_age, "^male_no_headaches") ~ "male_no_headaches",
      str_detect(sex_chronic_headache_age, "^male_headaches") ~ "male_headaches",
      str_detect(sex_chronic_headache_age, "^female_no_headaches") ~ "female_no_headaches",
      str_detect(sex_chronic_headache_age, "^female_headaches") ~ "female_headaches",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_headaches", "male_headaches", "female_no_headaches", "female_headaches"),
    labels = c("male without chronic or frequent headaches", 
               "male with chronic or frequent headaches", 
               "female without chronic or frequent headaches", 
               "female with chronic or frequent headaches")
    ),
    patchwork = factor(case_when(
      str_detect(sex_chronic_headache_age, "^male_no_headaches") ~ "male_no_headaches",
      str_detect(sex_chronic_headache_age, "^male_headaches") ~ "male_headaches",
      str_detect(sex_chronic_headache_age, "^female_no_headaches") ~ "female_no_headaches",
      str_detect(sex_chronic_headache_age, "^female_headaches") ~ "female_headaches",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_headaches", "male_headaches", "female_no_headaches", "female_headaches"),
    labels = c("male without health condition", 
               "male with health condition", 
               "female without health condition", 
               "female with health condition")
    )
  )

### additional data for result section
histogram_summary_sex_headaches_age <- 
  cate_histogram_data_sex_headaches_age |>
  group_by(sex_headaches_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - age, sex, and ptsd   ============================================
### table with effects across age, sex, and ptsd
cate_sex_ptsd_age <- CausalForestCATETable(
  cf,
  list(
    sex = rep(list("female" = 0, "male" = 1), each = 2, times = 5),
    ptsd = rep(list("ptsd" = 1, "no_ptsd" = 0), times = 10),
    age = c(rep(list("15-25" = 15:25), 4), 
            rep(list("26-35" = 26:35), 4), 
            rep(list("36-45" = 36:45), 4), 
            rep(list("46-55" = 46:55), 4),
            rep(list("56-65" = 56:65), 4))
  )
)

### t-tests for difference between ptsd and no ptsd
data <- bind_cols(oob, tibble(aipw = aipw_scores)) |>
  mutate(
    sex_ptsd_age = paste0(
      ifelse(sex == 0, "female_", "male_"),
      ifelse(ptsd == 0, "no_ptsd_", "ptsd_"),
      case_when(
        age %in% 15:25 ~ "15-25",
        age %in% 26:35 ~ "26-35",
        age %in% 36:45 ~ "36-45",
        age %in% 46:55 ~ "46-55",
        age %in% 56:65 ~ "56-65"
      )
    )
  ) |>
  select(sex_ptsd_age, aipw)
res <- tibble()
for(i in str_subset(unique(data$sex_ptsd_age), "no")) {
  tmp <- data |> 
    filter(
      str_detect(sex_ptsd_age, paste0("^",i)) | 
        str_detect(sex_ptsd_age, paste0("^", str_remove(i, "_no")))
    )
  ols <- lm(aipw ~ 1 + factor(sex_ptsd_age), data = tmp)
  res <- rbind(
    res, 
    as_tibble(
      coef(summary(ols))[2, c(1,2,4), drop=FALSE]
    ) %>%
      mutate(id = paste0(unique(tmp$sex_ptsd_age), collapse = " vs. "))  
  )
}
# Adjust p-values with Benjamini-Hockberg:
cate_sex_ptsd_age_ttest <- res %>%
  rename(`Orig. p-value` = `Pr(>|t|)`) %>%
  mutate(
    `95 % CI` = str_c("(",
                      round(Estimate - qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ", ",
                      round(Estimate + qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ")"),
    `Orig. p-value` = round(
      `Orig. p-value`,
      digits = 4),
    `Adj. p-value` = round(
      p.adjust(`Orig. p-value`, method = "BH"),
      digits = 4)) %>%
  select(id, Estimate, `Std. Error`, `95 % CI`, everything())

### data for plots with sex, ptsd, and age
cate_histogram_data_sex_ptsd_age <- oob |>
  mutate(
    sex_ptsd_age = factor(case_when(
      sex == 0 & ptsd == 1 & age < 50 ~ "female, <50 yrs, PTSD",
      sex == 0 & ptsd == 1 & age >= 50 ~ "female, \u226550 yrs, PTSD",
      sex == 0 & ptsd == 0 & age < 50 ~ "female, <50 yrs, no PTSD",
      sex == 0 & ptsd == 0 & age >= 50 ~ "female, \u226550 yrs, no PTSD",
      sex == 1 & ptsd == 1 & age < 50 ~ "male, < 50yrs, PTSD",
      sex == 1 & ptsd == 1 & age >= 50 ~ "male, \u226550 yrs, PTSD",
      sex == 1 & ptsd == 0 & age < 50 ~ "male, <50 yrs, no PTSD",
      sex == 1 & ptsd == 0 & age >= 50 ~ "male, \u226550 yrs, no PTSD",
      TRUE ~ "unknown"
    ),
    c(
      "male, <50 yrs, no PTSD", "female, <50 yrs, no PTSD", 
      "male, <50 yrs, PTSD", "female, <50 yrs, PTSD", 
      "male, \u226550 yrs, no PTSD", "female, \u226550 yrs, no PTSD", 
      "male, \u226550 yrs, PTSD", "female, \u226550 yrs, PTSD"
    )
    ),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_sex_ptsd_age <- 
  cate_histogram_data_sex_ptsd_age |>
  group_by(sex_ptsd_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.59,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

sex_ptsd_age_cate_table_lines_data <- 
  cate_sex_ptsd_age |>
  mutate(
    age_grp = case_when(
      str_detect(sex_ptsd_age, "15-25") ~ 20,
      str_detect(sex_ptsd_age, "26-35") ~ 30.5,
      str_detect(sex_ptsd_age, "36-45") ~ 40.5,
      str_detect(sex_ptsd_age, "46-55") ~ 50.5,
      str_detect(sex_ptsd_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    sex_ptsd = factor(case_when(
      str_detect(sex_ptsd_age, "^male_no_ptsd") ~ "male_no_ptsd",
      str_detect(sex_ptsd_age, "^male_ptsd") ~ "male_ptsd",
      str_detect(sex_ptsd_age, "^female_no_ptsd") ~ "female_no_ptsd",
      str_detect(sex_ptsd_age, "^female_ptsd") ~ "female_ptsd",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_ptsd", "male_ptsd", "female_no_ptsd", "female_ptsd"),
    labels = c("male without PTSD", 
               "male with PTSD", 
               "female without PTSD", 
               "female with PTSD")
    ),
    patchwork = factor(case_when(
      str_detect(sex_ptsd_age, "^male_no_ptsd") ~ "male_no_ptsd",
      str_detect(sex_ptsd_age, "^male_ptsd") ~ "male_ptsd",
      str_detect(sex_ptsd_age, "^female_no_ptsd") ~ "female_no_ptsd",
      str_detect(sex_ptsd_age, "^female_ptsd") ~ "female_ptsd",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_ptsd", "male_ptsd", "female_no_ptsd", "female_ptsd"),
    labels = c("male without health condition", 
               "male with health condition", 
               "female without health condition", 
               "female with health condition")
    )
  )

### additional data for result section
histogram_summary_sex_ptsd_age <- 
  cate_histogram_data_sex_ptsd_age |>
  group_by(sex_ptsd_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   CATE - age, sex, and diabetes   ========================================
### table with effects across age, sex, and diabetes
cate_sex_diabetes_age <- CausalForestCATETable(
  cf,
  list(
    sex = rep(list("female" = 0, "male" = 1), each = 2, times = 5),
    chronic_diabetes = rep(list("diabetes" = 1, "no_diabetes" = 0), times = 10),
    age = c(rep(list("15-25" = 15:25), 4), 
            rep(list("26-35" = 26:35), 4), 
            rep(list("36-45" = 36:45), 4), 
            rep(list("46-55" = 46:55), 4),
            rep(list("56-65" = 56:65), 4))
  )
)

### t-tests for difference between diabetes and no diabetes
data <- bind_cols(oob, tibble(aipw = aipw_scores)) |>
  mutate(
    sex_diabetes_age = paste0(
      ifelse(sex == 0, "female_", "male_"),
      ifelse(chronic_diabetes == 0, "no_diabetes_", "diabetes_"),
      case_when(
        age %in% 15:25 ~ "15-25",
        age %in% 26:35 ~ "26-35",
        age %in% 36:45 ~ "36-45",
        age %in% 46:55 ~ "46-55",
        age %in% 56:65 ~ "56-65"
      )
    )
  ) |>
  select(sex_diabetes_age, aipw)
res <- tibble()
for(i in str_subset(unique(data$sex_diabetes_age), "no")) {
  tmp <- data |> 
    filter(
      str_detect(sex_diabetes_age, paste0("^",i)) | 
        str_detect(sex_diabetes_age, paste0("^", str_remove(i, "_no")))
    )
  ols <- lm(aipw ~ 1 + factor(sex_diabetes_age), data = tmp)
  res <- rbind(
    res, 
    as_tibble(
      coef(summary(ols))[2, c(1,2,4), drop=FALSE]
    ) %>%
      mutate(id = paste0(unique(tmp$sex_diabetes_age), collapse = " vs. "))  
  )
}
# Adjust p-values with Benjamini-Hockberg:
cate_sex_diabetes_age_ttest <- res %>%
  rename(`Orig. p-value` = `Pr(>|t|)`) %>%
  mutate(
    `95 % CI` = str_c("(",
                      round(Estimate - qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ", ",
                      round(Estimate + qnorm(0.975) * `Std. Error`, 
                            digits = 4),
                      ")"),
    `Orig. p-value` = round(
      `Orig. p-value`,
      digits = 4),
    `Adj. p-value` = round(
      p.adjust(`Orig. p-value`, method = "BH"),
      digits = 4)) %>%
  select(id, Estimate, `Std. Error`, `95 % CI`, everything())

### data for plots with sex, diabetes and age
cate_histogram_data_sex_diabetes_age <- oob |>
  mutate(
    sex_diabetes_age = factor(case_when(
      sex == 0 & chronic_diabetes == 1 & age < 50 ~ "female, <50 yrs, diabetes",
      sex == 0 & chronic_diabetes == 1 & age >= 50 ~ "female, \u226550 yrs, diabetes",
      sex == 0 & chronic_diabetes == 0 & age < 50 ~ "female, <50 yrs, no diabetes",
      sex == 0 & chronic_diabetes == 0 & age >= 50 ~ "female, \u226550 yrs, no diabetes",
      sex == 1 & chronic_diabetes == 1 & age < 50 ~ "male, <50 yrs, diabetes",
      sex == 1 & chronic_diabetes == 1 & age >= 50 ~ "male, \u226550 yrs, diabetes",
      sex == 1 & chronic_diabetes == 0 & age < 50 ~ "male, <50 yrs, no diabetes",
      sex == 1 & chronic_diabetes == 0 & age >= 50 ~ "male, \u226550 yrs, no diabetes",
      TRUE ~ "unknown"
    ),
    c(
      "male, <50 yrs, no diabetes", "female, <50 yrs, no diabetes", 
      "male, <50 yrs, diabetes", "female, <50 yrs, diabetes", 
      "male, \u226550 yrs, no diabetes", "female, \u226550 yrs, no diabetes", 
      "male, \u226550 yrs, diabetes", "female, \u226550 yrs, diabetes"
    )
    ),
    ate_diff_signif = # TRUE is significantly different from ATE
      (ate_all$estimate < predictions - qnorm(0.975) * sqrt(variance.estimates)) |
      (ate_all$estimate > predictions + qnorm(0.975) * sqrt(variance.estimates))
  )
cate_histogram_label_sex_diabetes_age <- 
  cate_histogram_data_sex_diabetes_age |>
  group_by(sex_diabetes_age) |>
  summarise(
    count = n(),
    lower = sum(predictions < ate_all$estimate) / n(),
    higher = 1 - lower
  ) |>
  mutate(
    x_l = 100 * ate_all$estimate - 2.5,
    x_h = 100 * ate_all$estimate + 2.5,
    y = 0.59,
    lower = paste0(sprintf("%.1f", 100 * lower), "%"),
    higher = paste0(sprintf("%.1f", 100 * higher), "%")
  )

sex_diabetes_age_cate_table_lines_data <- 
  cate_sex_diabetes_age |>
  mutate(
    age_grp = case_when(
      str_detect(sex_chronic_diabetes_age, "15-25") ~ 20,
      str_detect(sex_chronic_diabetes_age, "26-35") ~ 30.5,
      str_detect(sex_chronic_diabetes_age, "36-45") ~ 40.5,
      str_detect(sex_chronic_diabetes_age, "46-55") ~ 50.5,
      str_detect(sex_chronic_diabetes_age, "56-65") ~ 60.5,
      TRUE ~ NA_real_
    ),
    sex_diabetes = factor(case_when(
      str_detect(sex_chronic_diabetes_age, "^male_no_diabetes") ~ "male_no_diabetes",
      str_detect(sex_chronic_diabetes_age, "^male_diabetes") ~ "male_diabetes",
      str_detect(sex_chronic_diabetes_age, "^female_no_diabetes") ~ "female_no_diabetes",
      str_detect(sex_chronic_diabetes_age, "^female_diabetes") ~ "female_diabetes",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_diabetes", "male_diabetes", "female_no_diabetes", "female_diabetes"),
    labels = c("male without diabetes", 
               "male with diabetes", 
               "female without diabetes", 
               "female with diabetes")
    ),
    patchwork = factor(case_when(
      str_detect(sex_chronic_diabetes_age, "^male_no_diabetes") ~ "male_no_diabetes",
      str_detect(sex_chronic_diabetes_age, "^male_diabetes") ~ "male_diabetes",
      str_detect(sex_chronic_diabetes_age, "^female_no_diabetes") ~ "female_no_diabetes",
      str_detect(sex_chronic_diabetes_age, "^female_diabetes") ~ "female_diabetes",
      TRUE ~ NA_character_
    ),
    levels = c("male_no_diabetes", "male_diabetes", "female_no_diabetes", "female_diabetes"),
    labels = c("male without health condition", 
               "male with health condition", 
               "female without health condition", 
               "female with health condition")
    )
  )

### additional data for result section
histogram_summary_sex_diabetes_age <- 
  cate_histogram_data_sex_diabetes_age |>
  group_by(sex_diabetes_age) |>
  summarise(
    count = n(),
    mean = mean(predictions),
    sd = sd(predictions),
    lower_signif = sum(predictions < ate_all$estimate & ate_diff_signif) / n(),
    higher_signif = sum(predictions > ate_all$estimate & ate_diff_signif) / n()
  )

# ===   data for plot with age, sex, and health condition interaction   ========
health_conditions_plot_data <- 
  sex_chronic_asthma_age_cate_table_lines_data |>
  select(2:4, 8, 10) |>
  structure(names = c("e1", "l1", "u1", "a", "c")) |>
  bind_cols(
    sex_headaches_age_cate_table_lines_data |>
      select(2:4) |>
      structure(names = c("e2", "l2", "u2"))
  ) |>
  bind_cols(
    sex_ptsd_age_cate_table_lines_data |>
      select(2:4) |>
      structure(names = c("e3", "l3", "u3"))
  ) |>
  bind_cols(
    sex_diabetes_age_cate_table_lines_data |>
      select(2:4) |>
      structure(names = c("e4", "l4", "u4"))
  )

# ===   deciles of conditional risk differences   ==============================
### risk difference subgroups
set.seed(seed)
dynamic_subgroups <- 
  CausalForestDynamicSubgroups(cf, 
                               num_rankings = 10, 
                               num_folds = 10, 
                               num.trees = 2000, 
                               num.threads = num_threads,
                               alpha = 0.01, 
                               min.node.size = 10,
                               seed = seed)

### forestplot with risk difference subgroups
forestplot_data <- dynamic_subgroups$cf_rank_ate |>
  mutate(
    lower = estimate - 1.96 * std_err,
    upper = estimate + 1.96 * std_err,
    `  ` = str_replace_all(
      sprintf("%.1f (%.1f to %.1f)", 100 * estimate, 100 * lower, 100 * upper),
      "\\.",
      "\\." 
    )
  )

# forest plot with deciles of conditional risk difference
forestplot_data <- forestplot_data |>
  mutate(
    `Predicted\nharm decile` = str_sub(ranking, 2, -1)
  ) |>
  inner_join(
    dynamic_subgroups$cf_subgroups$tau_hat |>
      as_tibble() |>
      rename(crd = value) |>
      arrange(crd) |>
      mutate(
        ranking = as.character(
          c(
            rep(1:8, each = ceiling(n() / 10)), 
            rep(9:10, each = floor(n() / 10))
          )
        )
      ) |>
      group_by(ranking) |>
      summarise(
        min = min(crd),
        max = max(crd),
        .groups = "drop"
      ),
    by = c("Predicted\nharm decile" = "ranking")
  ) |>
  inner_join(
    dynamic_subgroups$cf_subgroups$tau_hat |>
      as_tibble() |>
      mutate(
        test_result = as.numeric(data_retrospective_baseline_final_9m$test_result) - 1
      ) |>
      rename(crd = value) |>
      arrange(crd) |>
      mutate(
        ranking = as.character(
          c(
            rep(1:8, each = ceiling(n() / 10)), 
            rep(9:10, each = floor(n() / 10))
          )
        )
      ) |>
      count(ranking, test_result) |>
      mutate(
        test_result = ifelse(
          test_result, 
          "SARS-CoV-2\npositive", 
          "SARS-CoV-2\nnegative"
        )
      ) |>
      pivot_wider(names_from = test_result, values_from = n),
    by = c("Predicted\nharm decile" = "ranking")
  ) |>
  mutate(
    min = ifelse(min == min(min), -1, min),
    max = ifelse(max == max(max), 1, max),
    `Predicted CRD\ninterval` = str_replace_all(
      sprintf("%.3f to %.3f", min, max),
      "\\.",
      "\\." #\u00b7 to use middle dot
    ),
    ` ` = paste(rep(" ", 20), collapse = " ")
  ) |>
  select(
    estimate, std_err, lower, upper, 
    `Predicted\nharm decile`, `Predicted CRD\ninterval`, 
    `SARS-CoV-2\nnegative`, `SARS-CoV-2\npositive`, 
    `  `, ` `
  )

# ===   sensitivity analysis - RT-PCR test sensitivity   =======================
# time diff wrapper
time_diff <- function(t) {
  stopifnot("t must be a time or date-time object" = 
              inherits(t, c("POSIXct", "POSIXt")))
  time_diff <- Sys.time() - t
  time <- as.numeric(time_diff)
  unit <- attr(time_diff, "units")
  return(glue::glue("{sprintf('%.2f', time)} {unit}"))
}
test_sens <- 0.90
test_spec <- 0.99
# Numbers below determined by solving system of equations with number of test
# positive and negative as well as sensitivity and specificity
TN <- 47224
FP <- 477
FN <- 4112
TP <- 37005

pcr_sens_result <- list()
for(i in 1:20) {
  if (i == 1) t_t <- Sys.time()
  print(glue("Iteration {i}:"))
  t_i <- Sys.time()
  data_pcr_sens <- data_retrospective_baseline_final_9m
  ### sample false negatives and change them to positive
  data_pcr_sens$test_result[
    sample(
      which(data_retrospective_baseline_final_9m$test_result == "negative"),
      FN
    )
  ] <- "positive"
  ### sample false positives and change them to negative
  data_pcr_sens$test_result[
    sample(
      which(data_retrospective_baseline_final_9m$test_result == "positive"),
      FP
    )
  ] <- "negative"
  ### train causal forest using corrected exposure
  t_fy <- Sys.time()
  forest_Y_pcr_sens <- RegressionForestAnalysisWrapper(
    data = data_pcr_sens,
    cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", "anxiety", "depression", "ptsd"),
    binary = c(1, 3, 4, 5, 6),
    out = "Y", 
    num.trees = 2000, num.threads = 1L,
    alpha = 0.01, min.node.size = 10
  )
  print(glue("Total time: {time_diff(t_t)}, Y forest training time: {time_diff(t_fy)}"))
  t_fw <- Sys.time()
  forest_W_pcr_sens <- RegressionForestAnalysisWrapper(
    data = data_pcr_sens,
    cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", "anxiety", "depression", "ptsd"),
    binary = c(1, 3, 4, 5, 6),
    out = "W", 
    num.trees = 2000, num.threads = 1L,
    alpha = 0.01, min.node.size = 10
  )
  print(glue("Total time: {time_diff(t_t)}, W forest training time: {time_diff(t_fw)}"))
  Y_hat_pcr_sens <- predict(forest_Y_pcr_sens)$predictions
  W_hat_pcr_sens <- predict(forest_W_pcr_sens)$predictions
  t_fc <- Sys.time()
  cf_pcr_sens <- CausalForestAnalysisWrapper(
    data = data_pcr_sens,
    cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", "anxiety", "depression", "ptsd"),
    binary = c(1, 3, 4, 5, 6),
    Y.hat = Y_hat_pcr_sens, W.hat = W_hat_pcr_sens,
    num.trees = 2000, num.threads = 1L,
    alpha = 0.01, min.node.size = 10
  )
  print(glue("Total time: {time_diff(t_t)}, causal forest training time: {time_diff(t_fc)}"))
  ### variable importance
  variable_importance_pcr_sens <- tibble(
    variable_name = names(as_tibble(cf_pcr_sens$X.orig)),
    variable_importance = as.numeric(variable_importance(cf_pcr_sens))
  ) |>
    arrange(desc(variable_importance)) |>
    mutate(
      variable_importance_num = variable_importance,
      variable_importance = sprintf("%.3f", variable_importance)
    )
  
  ### CATE estimates by age, high BMI and depression
  cate_depression_high_bmi_age_pcr_sens <- CausalForestCATETable(
    cf_pcr_sens,
    list(
      depression = rep(list("no" = 0, "yes" = 1), each = 3, times = 5),
      high_bmi = rep(as.list(str_subset(names(cf_pcr_sens$X.orig), "^high_bmi")), times = 10),
      age = c(rep(list("15-25" = 15:25), 6), 
              rep(list("26-35" = 26:35), 6), 
              rep(list("36-45" = 36:45), 6), 
              rep(list("46-55" = 46:55), 6),
              rep(list("56-65" = 56:65), 6))
    )
  ) |>
    filter(!str_detect(depression_high_bmi_age, "unknown")) |>
    mutate(
      age_grp = case_when(
        str_detect(depression_high_bmi_age, "15-25") ~ 20,
        str_detect(depression_high_bmi_age, "26-35") ~ 30.5,
        str_detect(depression_high_bmi_age, "36-45") ~ 40.5,
        str_detect(depression_high_bmi_age, "46-55") ~ 50.5,
        str_detect(depression_high_bmi_age, "56-65") ~ 60.5,
        TRUE ~ NA_real_
      ),
      depression_high_bmi = factor(case_when(
        str_detect(depression_high_bmi_age, "no_no") ~ "no_no",
        str_detect(depression_high_bmi_age, "yes_no") ~ "yes_no",
        str_detect(depression_high_bmi_age, "no_yes") ~ "no_yes",
        str_detect(depression_high_bmi_age, "yes_yes") ~ "yes_yes",
        TRUE ~ NA_character_
      ),
      levels = c("no_no", "yes_no", "no_yes", "yes_yes"),
      labels = c("no depression, not high BMI", "depression, not high BMI", 
                 "no depression, high BMI", "depression, high BMI")
      )
    )
  
  ### ATE in full population and in subgroups along covariates
  ate_all_pcr_sens <- SubgroupATETable(NULL, NULL, cf_pcr_sens, 0.95)
  ate_age_pcr_sens <- c(map(seq(10, 60, 10), ~ .x + 0:9), list(14:49, 50:64)) |>
    map_dfr(\ (x) SubgroupATETable(x, "age", cf_pcr_sens, 0.95))
  ate_sex_pcr_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "sex", cf_pcr_sens, 0.95))
  ate_education_pcr_sens <- 
    stringr::str_subset(names(cf_pcr_sens$X.orig), "^status_education") |>
    map_dfr(
      \ (y) {
        SubgroupATETable(
          x = NULL,
          y = c(
            "subgroup", 
            paste0(
              "education - ", 
              stringr::str_replace_all(
                stringr::str_remove(y, "status_education_"), "_", " "
              )
            )
          ), 
          cf = cf_pcr_sens, 
          level = 0.95, 
          subset = cf_pcr_sens[["X.orig"]][[y]] == 1
        )
      }
    )
  ate_charlson_score_pcr_sens <- 0:3 |>
    map_dfr(\(x) SubgroupATETable(x, "charlson_score_cat", cf_pcr_sens, 0.95))
  ate_chronic_asthma_pcr_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_asthma", cf_pcr_sens, 0.95))
  ate_chronic_diabetes_pcr_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_diabetes", cf_pcr_sens, 0.95))
  ate_chronic_high_blood_pressure_pcr_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_high_blood_pressure", cf_pcr_sens, 0.95))
  ate_chronic_copd_lung_disease_pcr_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_copd_lung_disease", cf_pcr_sens, 0.95))
  ate_chronic_headache_pcr_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_headache", cf_pcr_sens, 0.95))
  ate_fibromyalgia_pcr_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "fibromyalgia", cf_pcr_sens, 0.95))
  ate_high_bmi_pcr_sens <- map2_dfr(
    rep(1, 3),
    str_subset(names(cf_pcr_sens$X.orig), "high_bmi"),
    \(x, y) SubgroupATETable(x, y, cf_pcr_sens, 0.95)
  )
  ate_chronic_fatigue_syndrome_pcr_sens <- 0:1 |>
    map_dfr(\(x) {
      SubgroupATETable(
        x, 
        "chronic_fatigue_syndrome", 
        cf_pcr_sens, 
        0.95)
    })
  ate_anxiety_pcr_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "anxiety", cf_pcr_sens, 0.95))
  ate_depression_pcr_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "depression", cf_pcr_sens, 0.95))
  ate_ptsd_pcr_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "ptsd", cf_pcr_sens, 0.95))
  ## combine ATE estimates
  ate_table_pcr_sens <- bind_rows(
    ate_all_pcr_sens,
    ate_age_pcr_sens,
    ate_sex_pcr_sens,
    ate_charlson_score_pcr_sens,
    ate_education_pcr_sens,
    ate_chronic_asthma_pcr_sens,
    ate_chronic_copd_lung_disease_pcr_sens,
    ate_chronic_diabetes_pcr_sens,
    ate_chronic_headache_pcr_sens,
    ate_chronic_high_blood_pressure_pcr_sens,
    ate_fibromyalgia_pcr_sens,
    ate_high_bmi_pcr_sens,
    ate_chronic_fatigue_syndrome_pcr_sens,
    ate_anxiety_pcr_sens,
    ate_depression_pcr_sens,
    ate_ptsd_pcr_sens
  ) |>
    mutate(
      `estimate (%)` = sprintf("%.1f", 100 * estimate),
      `95% CI - lower (%)` = sprintf("%.1f", 100 * `95% CI - lower`),
      `95% CI - upper (%)` = sprintf("%.1f", 100 * `95% CI - upper`)
    )
  
  ### Omnibus test based on linear model using demeaned covariates and outcomes
  test_calibration_pcr_sens <- test_calibration(cf_pcr_sens) |>
    structure(class = NULL) |>
    as_tibble() |>
    (\(x) bind_cols(
      tibble(` ` = c("mean forest prediction", "differential forest prediction")), 
      x
    ))()
  
  ### RATE based omnibus test
  t_ro <- Sys.time()
  rate_omnibus_test_pcr_sens <- 
    RATEOmnibusTest(cf_pcr_sens, alpha = 0.01, min.node.size = 10)
  print(glue("Total time: {time_diff(t_t)}, RATE omnibus test time: {time_diff(t_ro)}"))
  
  ### RATE based tests along covariates
  t_r <- Sys.time()
  rate_age_pcr_sens <- RATETest(cf_pcr_sens, "age")
  rate_sex_pcr_sens <- RATETest(cf_pcr_sens, "sex", "discrete")
  education_scores_pcr_sens <- tibble(
    education = case_when(
      cf_pcr_sens$X.orig$status_education_primary == 1 ~ "primary",
      cf_pcr_sens$X.orig$status_education_secondary == 1 ~ "secondary",
      cf_pcr_sens$X.orig$status_education_vocational_training == 1 ~ "vocational_training",
      cf_pcr_sens$X.orig$status_education_higher_short == 1 ~ "higher_short",
      cf_pcr_sens$X.orig$status_education_higher_medium == 1 ~ "higher_medium",
      cf_pcr_sens$X.orig$status_education_higher_long == 1 ~ "higher_long",
      cf_pcr_sens$X.orig$status_education_do_not_know == 1 ~ "do_not_know",
    ),
    dr_scores = get_scores(cf_pcr_sens),
  ) |>
    mutate(
      priority = case_when(
        education == "primary" ~ 6,
        education == "secondary" ~ 5,
        education == "vocational_training" ~ 4,
        education == "higher_short" ~ 3,
        education == "higher_medium" ~ 2,
        education == "higher_long" ~ 1,
        TRUE ~ NA_real_
      )
    ) |>
    filter(education != "do_not_know")
  education_q_pcr_sens <- c(
    0.001,
    cumsum(rev(table(education_scores_pcr_sens[["priority"]]))) / 
      length(education_scores_pcr_sens[["priority"]])
  )
  rate_education_pcr_sens <- list(
    rate = rank_average_treatment_effect.fit(
      DR.scores = education_scores_pcr_sens[["dr_scores"]],
      priorities = education_scores_pcr_sens[["priority"]],
      target = "AUTOC",
      q = education_q_pcr_sens,
      R = 500
    )) %>%
    c(list(
      confint = .$rate$estimate + 
        tibble(
          estimate = 0,
          lower = -qnorm(0.975) * .$rate$std.err,
          upper =  qnorm(0.975) * .$rate$std.err
        ),
      pval = 2 * pnorm(-abs(.$rate$estimate) / .$rate$std.err)
    ))
  rate_charlson_score_pcr_sens <- 
    RATETest(cf_pcr_sens, "charlson_score_cat", "discrete")
  rate_chronic_asthma_pcr_sens <- 
    RATETest(cf_pcr_sens, "chronic_asthma", "discrete")
  rate_chronic_diabetes_pcr_sens <- 
    RATETest(cf_pcr_sens, "chronic_diabetes", "discrete")
  rate_chronic_high_blood_pressure_pcr_sens <- 
    RATETest(cf_pcr_sens, "chronic_high_blood_pressure", "discrete")
  rate_chronic_copd_lung_disease_pcr_sens <- 
    RATETest(cf_pcr_sens, "chronic_copd_lung_disease", "discrete")
  rate_chronic_headache_pcr_sens <- 
    RATETest(cf_pcr_sens, "chronic_headache", "discrete")
  rate_fibromyalgia_pcr_sens <- 
    RATETest(cf_pcr_sens, "fibromyalgia", "discrete")
  high_bmi_scores_pcr_sens <- tibble(
    high_bmi = case_when(
      cf_pcr_sens$X.orig$high_bmi_yes == 1 ~ "yes",
      cf_pcr_sens$X.orig$high_bmi_no == 1 ~ "no",
      cf_pcr_sens$X.orig$high_bmi_unknown == 1 ~ "unknown",
    ),
    dr_scores = get_scores(cf_pcr_sens),
  ) |>
    mutate(
      priority = case_when(
        high_bmi == "yes" ~ 1,
        high_bmi == "no" ~ 0,
        TRUE ~ NA_real_
      )
    ) |>
    filter(high_bmi != "unknown")
  high_bmi_q_pcr_sens <- c(
    0.001,
    cumsum(rev(table(high_bmi_scores_pcr_sens[["priority"]]))) / 
      length(high_bmi_scores_pcr_sens[["priority"]])
  )
  rate_high_bmi_pcr_sens <- list(
    rate = rank_average_treatment_effect.fit(
      DR.scores = high_bmi_scores_pcr_sens[["dr_scores"]],
      priorities = high_bmi_scores_pcr_sens[["priority"]],
      target = "AUTOC",
      q = high_bmi_q_pcr_sens,
      R = 500
    )) %>%
    c(list(
      confint = .$rate$estimate + 
        tibble(
          estimate = 0,
          lower = -qnorm(0.975) * .$rate$std.err,
          upper =  qnorm(0.975) * .$rate$std.err
        ),
      pval = 2 * pnorm(-abs(.$rate$estimate) / .$rate$std.err)
    ))
  rate_chronic_fatigue_syndrome_pcr_sens <- 
    RATETest(cf_pcr_sens, "chronic_fatigue_syndrome", "discrete")
  rate_anxiety_pcr_sens <- RATETest(cf_pcr_sens, "anxiety", "discrete")
  rate_depression_pcr_sens <- RATETest(cf_pcr_sens, "depression", "discrete")
  rate_ptsd_pcr_sens <- RATETest(cf_pcr_sens, "ptsd", "discrete")
  
  ## combine RATE estimates
  rate_table_pcr_sens <- bind_rows(
    as_tibble(
      c(covariate = "omnibus", 
        rate_omnibus_test_pcr_sens$confint, 
        pval = rate_omnibus_test_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "age", 
        rate_age_pcr_sens$confint, 
        pval = rate_age_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "sex", 
        rate_sex_pcr_sens$confint, 
        pval = rate_sex_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "education", 
        rate_education_pcr_sens$confint, 
        pval = rate_education_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "charlson_score", 
        rate_charlson_score_pcr_sens$confint, 
        pval = rate_charlson_score_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "chronic_asthma", 
        rate_chronic_asthma_pcr_sens$confint, 
        pval = rate_chronic_asthma_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "chronic_copd_lung_disease", 
        rate_chronic_copd_lung_disease_pcr_sens$confint, 
        pval = rate_chronic_copd_lung_disease_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "chronic_diabetes", 
        rate_chronic_diabetes_pcr_sens$confint, 
        pval = rate_chronic_diabetes_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "chronic_headache", 
        rate_chronic_headache_pcr_sens$confint, 
        pval = rate_chronic_headache_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "chronic_high_blood_pressure", 
        rate_chronic_high_blood_pressure_pcr_sens$confint, 
        pval = rate_chronic_high_blood_pressure_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "fibromyalgia", 
        rate_fibromyalgia_pcr_sens$confint, 
        pval = rate_fibromyalgia_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "high_bmi", 
        rate_high_bmi_pcr_sens$confint, 
        pval = rate_high_bmi_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "chronic_fatigue_syndrome", 
        rate_chronic_fatigue_syndrome_pcr_sens$confint, 
        pval = rate_chronic_fatigue_syndrome_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "anxiety", 
        rate_anxiety_pcr_sens$confint, 
        pval = rate_anxiety_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "depression", 
        rate_depression_pcr_sens$confint, 
        pval = rate_depression_pcr_sens$pval)
    ),
    as_tibble(
      c(covariate = "ptsd", 
        rate_ptsd_pcr_sens$confint, 
        pval = rate_ptsd_pcr_sens$pval)
    )
  ) |>
    transmute(
      covariate = covariate,
      `RATE (AUTOC)` = sprintf("%.2E", estimate),
      `95% CI - lower` = sprintf("%.2E", lower),
      `95% CI - upper` = sprintf("%.2E", upper),
      `p-value` = sprintf("%.2E", pval)
    )
  print(glue("Total time: {time_diff(t_t)}, RATE test time: {time_diff(t_r)}"))
  ### save results in list
  pcr_sens_result[[i]] <- list(
    variable_importance = variable_importance_pcr_sens,
    cate_depression_high_bmi_age = cate_depression_high_bmi_age_pcr_sens,
    ate_table = ate_table_pcr_sens,
    test_calibration = test_calibration_pcr_sens,
    rate_table = rate_table_pcr_sens,
    rate_omnibus_test = rate_omnibus_test_pcr_sens
  )
  print(glue("Total time: {time_diff(t_t)}, iteration time: {time_diff(t_i)}"))
}

### prepare results
# calibration measures
pcr_sens_calibration <- map(pcr_sens_result, function(x) {
  tibble(
    mfp = x$test_calibration[1, 2, drop = TRUE],
    dfp = x$test_calibration[2, 2, drop = TRUE],
    dfp_pval = x$test_calibration[2, 5, drop = TRUE]
  )
}) |>
  list_rbind() |>
  summarise(
    mfp_mean = mean(mfp),
    mfp_median = median(mfp),
    mfp_sd = sd(mfp),
    mfp_min = min(mfp),
    mfp_max = max(mfp),
    dfp_mean = mean(dfp),
    dfp_median = median(dfp),
    dfp_sd = sd(dfp),
    dfp_min = min(dfp),
    dfp_max = max(dfp),
    dfp_pval_max = max(dfp_pval)
  )

# variable importance
pcr_sens_importance <- map(pcr_sens_result, function(x) {
  data <- x$variable_importance |>
    select(-variable_importance) |>
    bind_rows(
      x$variable_importance |>
        select(-variable_importance) |>
        filter(str_detect(variable_name, "high_bmi")) |>
        summarise(
          variable_name = "high_bmi",
          variable_importance_num = sum(variable_importance_num)
        )
    ) |>
    bind_rows(
      x$variable_importance |>
        select(-variable_importance) |>
        filter(str_detect(variable_name, "status_education")) |>
        summarise(
          variable_name = "education",
          variable_importance_num = sum(variable_importance_num)
        )
    ) |>
    filter(str_detect(variable_name, "status_education_", negate = TRUE)) |>
    filter(str_detect(variable_name, "high_bmi_", negate = TRUE)) |>
    mutate(
      variable_importance_ord = rank(-variable_importance_num, ties.method = "first")
    )
}) |>
  reduce(function(x, y) {
    inner_join(x, y, by = "variable_name")
  })

pcr_sens_importance_rank_plot_data <- pcr_sens_importance |>
  select(variable_name, starts_with("variable_importance_ord")) |>
  pivot_longer(cols = starts_with("variable_importance_ord"), values_to = "rank") |>
  select(variable_name, rank)  |>
  count(variable_name, rank) |>
  complete(variable_name, rank, fill = list(n = 0)) |>
  mutate(
    rank = factor(rank),
    variable_name = factor(
      variable_name,
      levels = c(
        "age", "high_bmi", "depression", "sex", "education", "chronic_asthma", 
        "charlson_score_cat", "ptsd", "chronic_headache", "chronic_diabetes", 
        "anxiety", "fibromyalgia", "chronic_high_blood_pressure", 
        "chronic_fatigue_syndrome", "chronic_copd_lung_disease"
      ),
      labels = c(
        "age", "high BMI", "depression", "sex", "education", "chronic asthma", 
        "Charlson Comorbidity Index", "post-traumatic stress disorder",
        "chronic or frequent headaches/migraines", "diabetes", "anxiety", 
        "fibromyalgia", "high blood pressure", 
        "chronic fatigue syndrome", "COPD or other chronic lung disease"
      )
    )
  )

# ATE for full population and single risk factors
pcr_sens_ate_summary <- map(
  pcr_sens_result,
  \(x) {
    x$ate_table |>
      select(subgroup, estimate) |>
      pivot_wider(names_from = "subgroup", values_from = "estimate")
  }
) |>
  list_rbind() |>
  summarise(
    across(
      everything(),
      list(
        mean = mean,
        median = median,
        sd = sd,
        min = min,
        max = max
      ),
      .names = "{.col} {.fn}"
    )
  ) |>
  pivot_longer(everything(), names_to = "name", values_to = "value")

# RATE for full population and single risk factors
pcr_sens_rate_summary <- map(
  pcr_sens_result,
  \(x) {
    x$rate_table |>
      select(covariate, `RATE (AUTOC)`) |>
      mutate(`RATE (AUTOC)` = 1000 * as.numeric(`RATE (AUTOC)`)) |>
      pivot_wider(names_from = "covariate", values_from = "RATE (AUTOC)")
  }
) |>
  list_rbind() |>
  summarise(
    across(
      everything(),
      list(
        mean = mean,
        median = median,
        sd = sd,
        min = min,
        max = max
      ),
      .names = "{.col} {.fn}"
    )
  ) |>
  pivot_longer(everything(), names_to = "name", values_to = "value")

pcr_sens_rate_pval_summary <- map(
  pcr_sens_result,
  \(x) {
    x$rate_table |>
      select(covariate, `p-value`) |>
      mutate(`p-value` = as.numeric(`p-value`)) |>
      pivot_wider(names_from = "covariate", values_from = "p-value")
  }
) |>
  list_rbind() |>
  summarise(
    across(
      everything(),
      list(
        mean = mean,
        median = median,
        sd = sd,
        min = min,
        max = max
      ),
      .names = "{.col} {.fn}"
    )
  ) |>
  pivot_longer(everything(), names_to = "name", values_to = "value")

# ===   sensitivity analysis - tuning parameters   =============================
data_tun_sens <- data_retrospective_baseline_final_9m
tuning_parameters <- expand_grid(
  sample.fraction = c(0.4, 0.5, 0.3),
  mtry = c(10, 15, 20),
  min.node.size = c(5, 10, 20),
  honesty.fraction = c(0.4, 0.5, 0.6),
)
tun_sens_result <- list()
for (i in 1:81) {
  if (i == 1) t_t <- Sys.time()
  print(glue("Iteration {i}:"))
  t_i <- Sys.time()
  t_fy <- Sys.time()
  forest_Y_tun_sens <- RegressionForestAnalysisWrapper(
    data = data_tun_sens,
    cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", 
            "anxiety", "depression", "ptsd"),
    binary = c(1, 3, 4, 5, 6),
    out = "Y", 
    num.trees = 2000, num.threads = 1L,
    alpha = 0.01, 
    sample.fraction = tuning_parameters[i, "sample.fraction", drop = TRUE],
    mtry = tuning_parameters[i, "mtry", drop = TRUE],
    min.node.size = tuning_parameters[i, "min.node.size", drop = TRUE],
    honesty.fraction = tuning_parameters[i, "honesty.fraction", drop = TRUE]
  )
  print(
    glue("Total time: {time_diff(t_t)}, Y forest training time: {time_diff(t_fy)}")
  )
  t_fw <- Sys.time()
  forest_W_tun_sens <- RegressionForestAnalysisWrapper(
    data = data_tun_sens,
    cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", 
            "anxiety", "depression", "ptsd"),
    binary = c(1, 3, 4, 5, 6),
    out = "W", 
    num.trees = 2000, num.threads = 1L,
    alpha = 0.01, 
    sample.fraction = tuning_parameters[i, "sample.fraction", drop = TRUE],
    mtry = tuning_parameters[i, "mtry", drop = TRUE],
    min.node.size = tuning_parameters[i, "min.node.size", drop = TRUE],
    honesty.fraction = tuning_parameters[i, "honesty.fraction", drop = TRUE]
  )
  print(
    glue("Total time: {time_diff(t_t)}, W forest training time: {time_diff(t_fw)}")
  )
  Y_hat_tun_sens <- predict(forest_Y_tun_sens)$predictions
  W_hat_tun_sens <- predict(forest_W_tun_sens)$predictions
  t_fc <- Sys.time()
  cf_tun_sens <- CausalForestAnalysisWrapper(
    data = data_tun_sens,
    cov = c("fibromyalgia", "high_bmi", "chronic_fatigue_syndrome", 
            "anxiety", "depression", "ptsd"),
    binary = c(1, 3, 4, 5, 6),
    Y.hat = Y_hat_tun_sens, W.hat = W_hat_tun_sens,
    num.trees = 2000, num.threads = 1L,
    alpha = 0.01, 
    sample.fraction = tuning_parameters[i, "sample.fraction", drop = TRUE],
    mtry = tuning_parameters[i, "mtry", drop = TRUE],
    min.node.size = tuning_parameters[i, "min.node.size", drop = TRUE],
    honesty.fraction = tuning_parameters[i, "honesty.fraction", drop = TRUE]
  )
  print(
    glue("Total time: {time_diff(t_t)}, causal forest training time: {time_diff(t_fc)}")
  )
  ### variable importance
  variable_importance_tun_sens <- tibble(
    variable_name = names(as_tibble(cf_tun_sens$X.orig)),
    variable_importance = as.numeric(variable_importance(cf_tun_sens))
  ) |>
    arrange(desc(variable_importance)) |>
    mutate(
      variable_importance_num = variable_importance,
      variable_importance = sprintf("%.3f", variable_importance)
    )
  ### CATE estimates by age, high BMI and depression
  cate_depression_high_bmi_age_tun_sens <- CausalForestCATETable(
    cf_tun_sens,
    list(
      depression = rep(list("no" = 0, "yes" = 1), each = 3, times = 5),
      high_bmi = rep(
        as.list(str_subset(names(cf_tun_sens$X.orig), "^high_bmi")), 
        times = 10
      ),
      age = c(rep(list("15-25" = 15:25), 6), 
              rep(list("26-35" = 26:35), 6), 
              rep(list("36-45" = 36:45), 6), 
              rep(list("46-55" = 46:55), 6),
              rep(list("56-65" = 56:65), 6))
    )
  ) |>
    filter(!str_detect(depression_high_bmi_age, "unknown")) |>
    mutate(
      age_grp = case_when(
        str_detect(depression_high_bmi_age, "15-25") ~ 20,
        str_detect(depression_high_bmi_age, "26-35") ~ 30.5,
        str_detect(depression_high_bmi_age, "36-45") ~ 40.5,
        str_detect(depression_high_bmi_age, "46-55") ~ 50.5,
        str_detect(depression_high_bmi_age, "56-65") ~ 60.5,
        TRUE ~ NA_real_
      ),
      depression_high_bmi = factor(case_when(
        str_detect(depression_high_bmi_age, "no_no") ~ "no_no",
        str_detect(depression_high_bmi_age, "yes_no") ~ "yes_no",
        str_detect(depression_high_bmi_age, "no_yes") ~ "no_yes",
        str_detect(depression_high_bmi_age, "yes_yes") ~ "yes_yes",
        TRUE ~ NA_character_
      ),
      levels = c("no_no", "yes_no", "no_yes", "yes_yes"),
      labels = c("no depression, not high BMI", "depression, not high BMI", 
                 "no depression, high BMI", "depression, high BMI")
      )
    )
  
  ### ATE in full population and in subgroups along covariates
  ate_all_tun_sens <- SubgroupATETable(NULL, NULL, cf_tun_sens, 0.95)
  ate_age_tun_sens <- c(map(seq(10, 60, 10), ~ .x + 0:9), list(14:49, 50:64)) |>
    map_dfr(\ (x) SubgroupATETable(x, "age", cf_tun_sens, 0.95))
  ate_sex_tun_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "sex", cf_tun_sens, 0.95))
  ate_education_tun_sens <- 
    stringr::str_subset(names(cf_tun_sens$X.orig), "^status_education") |>
    map_dfr(
      \ (y) {
        SubgroupATETable(
          x = NULL,
          y = c(
            "subgroup", 
            paste0(
              "education - ", 
              stringr::str_replace_all(
                stringr::str_remove(y, "status_education_"), "_", " "
              )
            )
          ), 
          cf = cf_tun_sens, 
          level = 0.95, 
          subset = cf_tun_sens[["X.orig"]][[y]] == 1
        )
      }
    )
  ate_charlson_score_tun_sens <- 0:3 |>
    map_dfr(\(x) SubgroupATETable(x, "charlson_score_cat", cf_tun_sens, 0.95))
  ate_chronic_asthma_tun_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_asthma", cf_tun_sens, 0.95))
  ate_chronic_diabetes_tun_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_diabetes", cf_tun_sens, 0.95))
  ate_chronic_high_blood_pressure_tun_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_high_blood_pressure", cf_tun_sens, 0.95))
  ate_chronic_copd_lung_disease_tun_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_copd_lung_disease", cf_tun_sens, 0.95))
  ate_chronic_headache_tun_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "chronic_headache", cf_tun_sens, 0.95))
  ate_fibromyalgia_tun_sens <- 0:1 |>
    map_dfr(\(x) SubgroupATETable(x, "fibromyalgia", cf_tun_sens, 0.95))
  ate_high_bmi_tun_sens <- map2_dfr(
    rep(1, 3),
    str_subset(names(cf_tun_sens$X.orig), "high_bmi"),
    \(x, y) SubgroupATETable(x, y, cf_tun_sens, 0.95)
  )
  ate_chronic_fatigue_syndrome_tun_sens <- 0:1 |>
    map_dfr(\(x) {
      SubgroupATETable(
        x, 
        "chronic_fatigue_syndrome", 
        cf_tun_sens, 
        0.95)
    })
  ate_anxiety_tun_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "anxiety", cf_tun_sens, 0.95))
  ate_depression_tun_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "depression", cf_tun_sens, 0.95))
  ate_ptsd_tun_sens <- 0:1 |>
    map_dfr(\ (x) SubgroupATETable(x, "ptsd", cf_tun_sens, 0.95))
  ## combine ATE estimates
  ate_table_tun_sens <- bind_rows(
    ate_all_tun_sens,
    ate_age_tun_sens,
    ate_sex_tun_sens,
    ate_charlson_score_tun_sens,
    ate_education_tun_sens,
    ate_chronic_asthma_tun_sens,
    ate_chronic_copd_lung_disease_tun_sens,
    ate_chronic_diabetes_tun_sens,
    ate_chronic_headache_tun_sens,
    ate_chronic_high_blood_pressure_tun_sens,
    ate_fibromyalgia_tun_sens,
    ate_high_bmi_tun_sens,
    ate_chronic_fatigue_syndrome_tun_sens,
    ate_anxiety_tun_sens,
    ate_depression_tun_sens,
    ate_ptsd_tun_sens
  ) |>
    mutate(
      `estimate (%)` = sprintf("%.1f", 100 * estimate),
      `95% CI - lower (%)` = sprintf("%.1f", 100 * `95% CI - lower`),
      `95% CI - upper (%)` = sprintf("%.1f", 100 * `95% CI - upper`)
    )
  
  ### Omnibus test based on linear model using demeaned covariates and outcomes
  test_calibration_tun_sens <- test_calibration(cf_tun_sens) |>
    structure(class = NULL) |>
    as_tibble() |>
    (\(x) bind_cols(
      tibble(` ` = c("mean forest prediction", "differential forest prediction")), 
      x
    ))()
  ### save results in list
  tun_sens_result[[i]] <- list(
    variable_importance = variable_importance_tun_sens,
    cate_depression_high_bmi_age = cate_depression_high_bmi_age_tun_sens,
    ate_table = ate_table_tun_sens,
    test_calibration = test_calibration_tun_sens
  )
  print(glue("Total time: {time_diff(t_t)}, iteration time: {time_diff(t_i)}"))
}

### prepare results
# calibration measures
tun_sens_calibration <- map(tun_sens_result, function(x) {
  tibble(
    mfp = x$test_calibration[1, 2, drop = TRUE],
    dfp = x$test_calibration[2, 2, drop = TRUE],
    dfp_pval = x$test_calibration[2, 5, drop = TRUE]
  )
}) |>
  list_rbind() |>
  summarise(
    mfp_mean = mean(mfp),
    mfp_median = median(mfp),
    mfp_sd = sd(mfp),
    mfp_min = min(mfp),
    mfp_max = max(mfp),
    dfp_mean = mean(dfp),
    dfp_median = median(dfp),
    dfp_sd = sd(dfp),
    dfp_min = min(dfp),
    dfp_max = max(dfp),
    dfp_pval_max = max(dfp_pval)
  )

tun_sens_importance <- map(tun_sens_result, function(x) {
  data <- x$variable_importance |>
    select(-variable_importance) |>
    bind_rows(
      x$variable_importance |>
        select(-variable_importance) |>
        filter(str_detect(variable_name, "high_bmi")) |>
        summarise(
          variable_name = "high_bmi",
          variable_importance_num = sum(variable_importance_num)
        )
    ) |>
    bind_rows(
      x$variable_importance |>
        select(-variable_importance) |>
        filter(str_detect(variable_name, "status_education")) |>
        summarise(
          variable_name = "education",
          variable_importance_num = sum(variable_importance_num)
        )
    ) |>
    filter(str_detect(variable_name, "status_education_", negate = TRUE)) |>
    filter(str_detect(variable_name, "high_bmi_", negate = TRUE)) |>
    mutate(
      variable_importance_ord = rank(-variable_importance_num, ties.method = "first")
    )
}) |>
  reduce(function(x, y) {
    inner_join(x, y, by = "variable_name")
  })

tun_sens_importance_rank_plot_data <- tun_sens_importance|>
  select(variable_name, starts_with("variable_importance_ord")) |>
  pivot_longer(cols = starts_with("variable_importance_ord"), values_to = "rank") |>
  select(variable_name, rank)  |>
  count(variable_name, rank) |>
  complete(variable_name, rank, fill = list(n = 0)) |>
  mutate(
    rank = factor(rank),
    variable_name = factor(
      variable_name,
      levels = c(
        "age", "high_bmi", "depression", "sex", "education", "chronic_asthma", 
        "charlson_score_cat", "ptsd", "chronic_headache", "chronic_diabetes", 
        "anxiety", "fibromyalgia", "chronic_high_blood_pressure", 
        "chronic_fatigue_syndrome", "chronic_copd_lung_disease"
      ),
      labels = c(
        "age", "high BMI", "depression", "sex", "education", "chronic asthma", 
        "Charlson Comorbidity Index", "post-traumatic stress disorder",
        "chronic or frequent headaches/migraines", "diabetes", "anxiety", 
        "fibromyalgia", "high blood pressure", 
        "chronic fatigue syndrome", "COPD or other chronic lung disease"
      )
    )
  )

# ATE for full population and single risk factors
tun_sens_ate_summary <- map(
  tun_sens_result,
  \(x) {
    x$ate_table |>
      select(subgroup, estimate) |>
      pivot_wider(names_from = "subgroup", values_from = "estimate")
  }
) |>
  list_rbind() |>
  summarise(
    across(
      everything(),
      list(
        mean = mean,
        median = median,
        sd = sd,
        min = min,
        max = max
      ),
      .names = "{.col} {.fn}"
    )
  ) |>
  pivot_longer(everything(), names_to = "name", values_to = "value") |>
  mutate(percent = sprintf("%.3f", 100 * value))

# ===   save data to reproduce plots   =========================================
saveRDS(
  variable_importance |>
    select(-variable_importance) |>
    bind_rows(
      variable_importance |>
        select(-variable_importance) |>
        filter(str_detect(variable_name, "high_bmi")) |>
        summarise(
          variable_name = "high_bmi",
          variable_importance_num = sum(variable_importance_num)
        )
    ) |>
    bind_rows(
      variable_importance |>
        select(-variable_importance) |>
        filter(str_detect(variable_name, "status_education")) |>
        summarise(
          variable_name = "education",
          variable_importance_num = sum(variable_importance_num)
        )
    ) |>
    filter(str_detect(variable_name, "status_education_", negate = TRUE)) |>
    filter(str_detect(variable_name, "high_bmi_", negate = TRUE)) |>
    mutate(
      variable_name = factor(
        variable_name,
        levels = rev(c(
          "age", "high_bmi", "depression", "sex", "education", "chronic_asthma", 
          "charlson_score_cat", "ptsd", "chronic_headache", "chronic_diabetes", 
          "anxiety", "fibromyalgia", "chronic_high_blood_pressure", 
          "chronic_fatigue_syndrome", "chronic_copd_lung_disease"
        )),
        labels = rev(c(
          "age", "high BMI", "depression", "sex", "education", "chronic asthma", 
          "Charlson Comorbidity Index", "post-traumatic stress disorder",
          "chronic or frequent headaches/migraines", "diabetes", "anxiety", 
          "fibromyalgia", "high blood pressure", 
          "chronic fatigue syndrome", "COPD or other chronic lung disease"
        )
        )) |>
        fct_reorder(variable_importance_num)
    ),
  file = "../variable_importance_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = tibble(
    W_hat = cf$W.hat, 
    W_orig = factor(
      cf$W.orig, 
      levels = c(0, 1), 
      labels = c("Negative", "Positive")
    ),
    IPW = ifelse(
      cf$W.orig == 1,
      1 / cf$W.hat,
      1 / (1 - cf$W.hat)
    )
  ),
  file = "../overlap_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = ate_table,
  file = "../ate_table.Rdata",
  compress = "xz"
)
saveRDS(
  object = oob$predictions,
  file = "../cate_predictions.Rdata",
  compress = "xz"
)
saveRDS(
  object = balance$love_data,
  file = "../love_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = balance_manu_plot$ecdf_data,
  file = "../ecdf_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_depression_high_bmi_age |>
      select(
        depression_high_bmi_age,
        predictions, 
        variance.estimates,
        ate_diff_signif
      ),
    label = cate_histogram_label_depression_high_bmi_age
  ),
  file = "../cate_histogram_data_depression_high_bmi_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_sex_high_bmi_age |>
      select(
        sex_high_bmi_age, 
        predictions, 
        variance.estimates, 
        ate_diff_signif
      ),
    label = cate_histogram_label_sex_high_bmi_age
  ),
  file = "../cate_histogram_data_sex_high_bmi_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_sex_depression_age |>
      select(
        sex_depression_age, 
        predictions, 
        variance.estimates, 
        ate_diff_signif
      ),
    label = cate_histogram_label_sex_depression_age
  ),
  file = "../cate_histogram_data_sex_depression_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_high_bmi_depression_sex |>
      select(
        high_bmi_depression_sex, 
        predictions, 
        variance.estimates, 
        ate_diff_signif
      ),
    label = cate_histogram_label_high_bmi_depression_sex
  ),
  file = "../cate_histogram_data_depression_high_bmi_sex.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_sex_asthma_age |>
      select(
        sex_asthma_age, 
        predictions, 
        variance.estimates, 
        ate_diff_signif
      ),
    label = cate_histogram_label_sex_asthma_age
  ),
  file = "../cate_histogram_data_sex_asthma_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_sex_headaches_age |>
      select(
        sex_headaches_age, 
        predictions, 
        variance.estimates, 
        ate_diff_signif
      ),
    label = cate_histogram_label_sex_headaches_age
  ),
  file = "../cate_histogram_data_sex_headaches_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = list(
    data = cate_histogram_data_sex_ptsd_age |>
      select(
        sex_ptsd_age, 
        predictions, 
        variance.estimates, 
        ate_diff_signif
      ),
    label = cate_histogram_label_sex_ptsd_age
  ),
  file = "../cate_histogram_data_sex_ptsd_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = depression_high_bmi_age_cate_table_lines_data,
  file = "../cate_lines_data_depression_high_bmi_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = sex_high_bmi_age_cate_table_lines_data,
  file = "../cate_lines_data_sex_high_bmi_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = sex_depression_age_cate_table_lines_data,
  file = "../cate_lines_data_sex_depression_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = high_bmi_depression_sex_cate_table_data,
  file = "../cate_lines_data_depression_high_bmi_sex.Rdata",
  compress = "xz"
)
saveRDS(
  object = sex_chronic_asthma_age_cate_table_lines_data,
  file = "../cate_lines_data_sex_asthma_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = sex_headaches_age_cate_table_lines_data,
  file = "../cate_lines_data_sex_headaches_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = sex_ptsd_age_cate_table_lines_data,
  file = "../cate_lines_data_sex_ptsd_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = sex_depression_age_cate_table_lines_data,
  file = "../cate_lines_data_sex_depression_age.Rdata",
  compress = "xz"
)
saveRDS(
  object = forestplot_data,
  file = "../forestplot_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = pcr_sens_importance,
  file = "../pcr_sens_importance.Rdata",
  compress = "xz"
)
saveRDS(
  object = tun_sens_importance,
  file = "../tun_sens_importance.Rdata",
  compress = "xz"
)
saveRDS(
  object = pcr_sens_importance_rank_plot_data,
  file = "../pcr_sens_importance_rank_plot_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = tun_sens_importance_rank_plot_data,
  file = "../tun_sens_importance_rank_plot_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = health_conditions_plot_data,
  file = "../health_conditions_plot_data.Rdata",
  compress = "xz"
)
saveRDS(
  object = pcr_sens_result,
  file = "../results/data/pcr_sens_result.Rdata",
  compress = "xz"
)
saveRDS(
  object = tun_sens_result,
  file = "../results/data/tun_sens_result.Rdata",
  compress = "xz"
)
