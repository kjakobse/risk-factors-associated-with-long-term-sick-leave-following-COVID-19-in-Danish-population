# ===   setup   ================================================================
### User inputs
setwd("path/to/plot/data")
output_path <- "path/to/output"
### required packages
library(cowplot)
library(DiagrammeR)
library(DiagrammeRsvg)
library(forestploter)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(grid)
library(gtable)
library(mmtable2)
library(patchwork)
library(ragg)
library(tidyverse)
library(writexl)
### plot theme for publication
theme_set(theme_pubr())
### helper functions
get_scale <- function(plot,
                      width_wanted,
                      height_wanted,
                      unit = "in"){
  h <- convertHeight(sum(plot$heights), unit, TRUE)
  w <- convertWidth(sum(plot$widths), unit, TRUE)
  max(c(w/width_wanted,  h/height_wanted))
}
floor_dec <- function(x, digits = 1) {
  round(x - 5 * 10^(-digits - 1), digits)
}
ceiling_dec <- function(x, digits = 1) {
  round(x + 5 * 10^(-digits - 1), digits)
}
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    0
  }
}
### read in data
variable_importance_data <- 
  readRDS("variable_importance_data.Rdata")
overlap_data <- 
  readRDS("overlap_data.Rdata")
love_data <- 
  readRDS("love_data.Rdata")
ecdf_data <- 
  readRDS("ecdf_data.Rdata")
ate_table <- 
  readRDS("ate_table.Rdata")
cate_predictions <- 
  readRDS("cate_predictions.Rdata")
forestplot_data <- 
  readRDS("forestplot_data.Rdata")
cate_histogram_data_depression_high_bmi_age <- 
  readRDS("cate_histogram_data_depression_high_bmi_age.Rdata")
cate_histogram_data_sex_high_bmi_age <- 
  readRDS("cate_histogram_data_sex_high_bmi_age.Rdata")
cate_histogram_data_sex_depression_age <- 
  readRDS("cate_histogram_data_sex_depression_age.Rdata")
cate_histogram_data_depression_high_bmi_sex <- 
  readRDS("cate_histogram_data_depression_high_bmi_sex.Rdata")
cate_histogram_data_sex_asthma_age <- 
  readRDS("cate_histogram_data_sex_asthma_age.Rdata")
cate_lines_data_depression_high_bmi_age <- 
  readRDS("cate_lines_data_depression_high_bmi_age.Rdata")
cate_lines_data_sex_high_bmi_age <- 
  readRDS("cate_lines_data_sex_high_bmi_age.Rdata")
cate_lines_data_sex_depression_age <- 
  readRDS("cate_lines_data_sex_depression_age.Rdata")
cate_lines_data_depression_high_bmi_sex <- 
  readRDS("cate_lines_data_depression_high_bmi_sex.Rdata")
cate_lines_data_sex_asthma_age <- 
  readRDS("cate_lines_data_sex_asthma_age.Rdata")
tree_plot <- 
  readRDS("tree_plot.Rdata")
pcr_sens_importance <- 
  readRDS("pcr_sens_importance.Rdata")
tun_sens_importance <- 
  readRDS("tun_sens_importance.Rdata")
pcr_sens_importance_rank_plot_data <- 
  readRDS("pcr_sens_importance_rank_plot_data.Rdata")
tun_sens_importance_rank_plot_data <- 
  readRDS("tun_sens_importance_rank_plot_data.Rdata")
health_conditions_plot_data <- 
  readRDS("health_conditions_plot_data.Rdata")
pcr_sens_result <- 
  readRDS("pcr_sens_result.Rdata")
tun_sens_result <- 
  readRDS("tun_sens_result.Rdata")

# ===   plots   ================================================================
### Figure 2 (forest plot + cate density)
Figure_2_b <- ggplot(tibble(tau_hat = cate_predictions)) + 
  geom_density(aes(x = 100 * tau_hat, y = after_stat(density)),
               color = "black", fill = "#eff3f2", linewidth = 0.25) +
  geom_vline(xintercept = 100 * ate_table[[2]][1], linetype = 2, linewidth = 0.3) +
  coord_cartesian(xlim = c(-5, 15)) +
  xlab("risk difference (percentage points)") + ylab("density") +
  theme(
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    text = element_text(family = "sans", size = 6.5)
  )

fp_theme <- forest_theme(
  base_size = 6.5,
  base_family = "sans",
  ci_pch = 20,
  ci_col = "black",
  ci_fill = "black",
  ci_lty = 1,
  ci_lwd = 1.5,
  ci_Theight = 0.2
)
Figure_2_a <- forest(
  forestplot_data[, 5:10],
  est = 100 * forestplot_data$estimate,
  lower = 100 * forestplot_data$lower,
  upper = 100 * forestplot_data$upper,
  sizes = 0.3,
  ci_column = 6,
  ref_line = 0,
  arrow_lab = c("Protective", "Harmful"),
  xlim = c(-5, 10),
  ticks_at = c(-3, 0, 3, 6, 9),
  theme = fp_theme
)
Figure_2_a <- add_text(
  Figure_2_a, 
  text = "Risk Difference, % (95% CI)",
  part = "header", 
  col = 5:6,
  gp = gpar(fontface = "bold", fontsize = 7, fontfamily = "sans")
)
Figure_2_a <- insert_text(
  Figure_2_a, 
  text = "No. of participants",
  part = "header",
  col = 3:4,
  gp = gpar(fontface = "bold", fontsize = 7, fontfamily = "sans"),
  just = "left"
)
Figure_2_a <- add_border(
  Figure_2_a, 
  part = "header",
  row = 1,
  col = 3:4,
  gp = gpar(lwd = .8)
)
Figure_2_a <- add_border(
  Figure_2_a, 
  part = "header",
  row = 2
)
forestplot_sc <- get_scale(Figure_2_a, 12, 8, unit = "cm")

Figure_2 <- cowplot::plot_grid(
  Figure_2_a, 
  Figure_2_b, 
  labels = c("a", "b"), 
  label_size = 8,
  rel_widths = c(2, 1),
  rel_heights = c(1, 1)
)

### Figure 3
cate_histogram_depression_high_bmi_age <- 
  ggplot(cate_histogram_data_depression_high_bmi_age$data) + 
  geom_histogram(aes(x = 100 * predictions,
                     y = after_stat(count) * 
                       5 / 
                       rep(
                         cate_histogram_data_depression_high_bmi_age$label$count, 
                         each = 2 * 87
                       ),
                     alpha = ate_diff_signif),
                 binwidth = 0.2,
                 fill = "red3") +
  geom_vline(xintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.25) +
  geom_label(aes(x = x_l, y = y, label = lower), vjust = 0.5, nudge_y = 0,
             data = cate_histogram_data_depression_high_bmi_age$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  geom_label(aes(x = x_h, y = y, label = higher), vjust = 0.5, nudge_y = 0,
             data = cate_histogram_data_depression_high_bmi_age$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  facet_wrap(~ depression_high_bmi_age, nrow = 4) +
  coord_cartesian(xlim = c(-5, 15)) + 
  xlab("conditional risk difference (percentage points)") + ylab("density") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.5)) +
  theme(legend.position = "none",
        text = element_text(family = "sans", size = 7.5),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
        strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

cate_lines_depression_high_bmi_age <- 
  ggplot(cate_lines_data_depression_high_bmi_age) +
  geom_line(aes(y = 100 * estimate, x = age_grp, color = depression_high_bmi),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * `95% CI - lower`, 
                    ymax = 100 * `95% CI - upper`,
                    x = age_grp, 
                    color = depression_high_bmi),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * estimate, x = age_grp, color = depression_high_bmi),
             position = position_dodge(width = 2.5), size = 1) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.3) +
  ylab("risk difference (percentage points)") + xlab("age (yrs)") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-5, 15)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        text = element_text(family = "sans", size = 7.5))

Figure_3 <- cowplot::plot_grid(
  cate_lines_depression_high_bmi_age, 
  cate_histogram_depression_high_bmi_age, 
  labels = c("a", "b"), 
  label_size = 8,
  rel_widths = c(1, 1),
  rel_heights = c(7, 9)
)

subgroup_counts_depression_high_bmi_age_data <- 
  cate_lines_data_depression_high_bmi_age |>
  select(depression_high_bmi_age, n) |>
  mutate(
    depression_high_bmi = case_when(
      str_detect(depression_high_bmi_age, "no_no") ~ "no depression, not high BMI",
      str_detect(depression_high_bmi_age, "yes_no") ~ "depression, not high BMI",
      str_detect(depression_high_bmi_age, "no_yes") ~ "no depression, high BMI",
      str_detect(depression_high_bmi_age, "yes_yes") ~ "depression, high BMI",
    ),
    age = paste(str_extract(depression_high_bmi_age, "\\d\\d-\\d\\d"), "yrs")
  ) |>
  select(-depression_high_bmi_age) |>
  pivot_wider(names_from = "depression_high_bmi", values_from = "n") |>
  bind_rows(
    cate_histogram_data_depression_high_bmi_age$data |>
      count(depression_high_bmi_age) |>
      mutate(
        depression_high_bmi = str_remove(depression_high_bmi_age, "^.{9}"),
        age = str_sub(depression_high_bmi_age, 1, 7)
      ) |>
      select(-depression_high_bmi_age) |>
      pivot_wider(names_from = "depression_high_bmi", values_from = "n")
  ) |>
  rename(" " = "age")
subgroup_counts_depression_high_bmi_age <- 
  subgroup_counts_depression_high_bmi_age_data |>
  tableGrob(
    theme = ttheme_minimal(
      base_size = 7, 
      base_family = "sans",
      core = list(
        bg_params = list(fill = c(rep(c("#eff3f2", "white"), 3), "#eff3f2")),
        fg_params = list(
          hjust = 0, 
          x = unit(0, "cm"),
          fontface = c(rep("bold", 7), rep("plain", 28))
        )
      ),
      colhead = list(fg_params = list(hjust = 0, x = unit(0, "cm")))
    )
  ) |>
  gtable_remove_grobs("rowhead-fg") |>
  gtable_add_grob(linesGrob(y = c(0, 0), gp = gpar(lwd = 2, col = "black")), t = 1, l = 2, r = 6) |>
  gtable_squash_cols(1)

Figure_3_t <- cowplot::plot_grid(
  Figure_3,
  subgroup_counts_depression_high_bmi_age,
  labels = c("", "c"), 
  label_size = 8,
  ncol = 1,
  rel_heights = c(9, 5.5)
)

### Figure 4
Figure_4_b <- 
  ggplot(cate_histogram_data_depression_high_bmi_sex$data) + 
  geom_histogram(aes(x = 100 * predictions,
                     y = after_stat(count) * 
                       5 / 
                       rep(
                         cate_histogram_data_depression_high_bmi_sex$label$count, 
                         each = 2 * 87
                       ),
                     alpha = ate_diff_signif), 
                 binwidth = 0.2,
                 fill = "red3") +
  geom_vline(xintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.25) +
  geom_label(aes(x = x_l, y = y, label = lower), vjust = 0.5,
             data = cate_histogram_data_depression_high_bmi_sex$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  geom_label(aes(x = x_h, y = y, label = higher), vjust = 0.5,
             data = cate_histogram_data_depression_high_bmi_sex$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  facet_wrap(~ high_bmi_depression_sex, nrow = 4) +
  coord_cartesian(xlim = c(-5, 15)) + 
  xlab("conditional risk difference (percentage points)") + ylab("density") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0, 0.43)) +
  theme(legend.position = "none",
        text = element_text(family = "sans", size = 7.5),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
        strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

Figure_4_a <- 
  cate_lines_data_depression_high_bmi_sex |>
  mutate(high_bmi_depression_sex = fct_rev(high_bmi_depression_sex)) |>
  ggplot() +
  geom_point(aes(x = 100 * estimate, y = high_bmi_depression_sex), 
             color = "red3", size = 1) +
  geom_errorbarh(aes(xmin = 100 * `95% CI - lower`, 
                     xmax = 100 * `95% CI - upper`,
                     y = high_bmi_depression_sex),
                 height = 0, color = "red3", linewidth = 0.7) + 
  geom_vline(xintercept = 100 * ate_table[1,2]$estimate, linetype = 2, linewidth = 0.3) +
  xlab("risk difference (percentage points)") + ylab("") +
  theme(axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        text = element_text(family = "sans", size = 7.5))

Figure_4 <- cowplot::plot_grid(
  Figure_4_a, 
  Figure_4_b, 
  labels = c("a", "b"), 
  label_size = 8,
  rel_widths = c(1, 1),
  rel_heights = c(7, 9)
)

subgroup_counts_depression_high_bmi_sex_data <- 
  cate_lines_data_depression_high_bmi_sex |>
  select(high_bmi_depression_sex, n) |>
  mutate(
    depression_high_bmi = case_when(
      str_detect(high_bmi_depression_sex, "no depression, not high BMI") ~ "no depression, not high BMI",
      str_detect(high_bmi_depression_sex, "depression, not high BMI") ~ "depression, not high BMI",
      str_detect(high_bmi_depression_sex, "no depression, high BMI") ~ "no depression, high BMI",
      str_detect(high_bmi_depression_sex, "depression, high BMI") ~ "depression, high BMI",
    ),
    sex = ifelse(str_detect(high_bmi_depression_sex, "^female"), "female", "male")
  ) |>
  select(-high_bmi_depression_sex) |>
  pivot_wider(names_from = "depression_high_bmi", values_from = "n") |>
  rename(" " = "sex") 

subgroup_counts_depression_high_bmi_sex <-
  subgroup_counts_depression_high_bmi_sex_data |>
  tableGrob(
    theme = ttheme_minimal(
      base_size = 7, 
      base_family = "sans",
      core = list(
        bg_params = list(fill = c("#eff3f2", "white")),
        fg_params = list(
          hjust = 0, 
          x = unit(0, "cm"),
          fontface = c(rep("bold", 2), rep("plain", 8))
        )
      ),
      colhead = list(fg_params = list(hjust = 0, x = unit(0, "cm")))
    )
  ) |>
  gtable_remove_grobs("rowhead-fg") |>
  gtable_add_grob(linesGrob(y = c(0, 0), gp = gpar(lwd = 2, col = "black")), t = 1, l = 2, r = 6) |>
  gtable_squash_cols(1)

Figure_4_t <- cowplot::plot_grid(
  Figure_4,
  subgroup_counts_depression_high_bmi_sex,
  labels = c("", "c"), 
  label_size = 8,
  ncol = 1,
  rel_heights = c(9, 2.2)
)

### Figure 5
Figure_5_a <- 
  ggplot(cate_lines_data_depression_high_bmi_age) +
  geom_line(aes(y = 100 * estimate, x = age_grp, color = depression_high_bmi),
            linewidth = 1) +
  ylab("risk difference (percentage points)") + xlab("age (yrs)") +
  coord_cartesian(xlim = c(18, 62), ylim = c(-1, 12)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(family = "sans", size = 7.5))
for(i in seq_along(pcr_sens_result)) {
  pcr_sens_depression_high_bmi_age_cate_table_lines_data <-
    pcr_sens_result[[i]]$cate_depression_high_bmi_age |>
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
  
  Figure_5_a <- 
    Figure_5_a +
    geom_line(aes(y = 100 * estimate, x = age_grp, color = depression_high_bmi),
              pcr_sens_depression_high_bmi_age_cate_table_lines_data,
              linewidth = 0.7, alpha = 0.3)
  
  if(i == 1) {
    figure_5a_data <- pcr_sens_depression_high_bmi_age_cate_table_lines_data |>
      mutate(iteration = i) |>
      select(
        iteration, 
        depression_high_bmi_age, 
        estimate, 
        `95% CI - lower`, 
        `95% CI - upper`
      ) |>
      rename("Risk difference" = "estimate")
  } else {
    figure_5a_data <- figure_5a_data |>
      bind_rows(
        pcr_sens_depression_high_bmi_age_cate_table_lines_data |>
          mutate(iteration = i) |>
          select(
            iteration, 
            depression_high_bmi_age, 
            estimate, 
            `95% CI - lower`, 
            `95% CI - upper`
          ) |>
          rename("Risk difference" = "estimate")
      )
  }
}

Figure_5_b <- 
  ggplot(cate_lines_data_depression_high_bmi_age) +
  geom_line(aes(y = 100 * estimate, x = age_grp, color = depression_high_bmi),
            linewidth = 1) +
  ylab("risk difference (percentage points)") + xlab("age (yrs)") +
  coord_cartesian(xlim = c(18, 62), ylim = c(-1, 12)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(family = "sans", size = 7.5))
for(i in seq_along(tun_sens_result)) {
  tun_sens_depression_high_bmi_age_cate_table_lines_data <-
    tun_sens_result[[i]]$cate_depression_high_bmi_age |>
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
  
  Figure_5_b <- 
    Figure_5_b +
    geom_line(aes(y = 100 * estimate, x = age_grp, color = depression_high_bmi),
              tun_sens_depression_high_bmi_age_cate_table_lines_data,
              linewidth = 0.7, alpha = 0.3)
  
  if(i == 1) {
    figure_5b_data <- tun_sens_depression_high_bmi_age_cate_table_lines_data |>
      mutate(iteration = i) |>
      select(
        iteration, 
        depression_high_bmi_age, 
        estimate, 
        `95% CI - lower`, 
        `95% CI - upper`
      ) |>
      rename("Risk difference" = "estimate")
  } else {
    figure_5b_data <- figure_5b_data |>
      bind_rows(
        tun_sens_depression_high_bmi_age_cate_table_lines_data |>
          mutate(iteration = i) |>
          select(
            iteration, 
            depression_high_bmi_age, 
            estimate, 
            `95% CI - lower`, 
            `95% CI - upper`
          ) |>
          rename("Risk difference" = "estimate")
      )
  }
}

Figure_5 <- cowplot::plot_grid(
  Figure_5_a, 
  Figure_5_b, 
  labels = c("a", "b"), 
  label_size = 8,
  rel_widths = c(1, 1)
)

subgroup_counts_depression_high_bmi_age_data_2 <- 
  cate_lines_data_depression_high_bmi_age |>
  select(depression_high_bmi_age, n) |>
  mutate(
    depression_high_bmi = case_when(
      str_detect(depression_high_bmi_age, "no_no") ~ "no depression, not high BMI",
      str_detect(depression_high_bmi_age, "yes_no") ~ "depression, not high BMI",
      str_detect(depression_high_bmi_age, "no_yes") ~ "no depression, high BMI",
      str_detect(depression_high_bmi_age, "yes_yes") ~ "depression, high BMI",
    ),
    age = paste(str_extract(depression_high_bmi_age, "\\d\\d-\\d\\d"), "yrs")
  ) |>
  select(-depression_high_bmi_age) |>
  pivot_wider(names_from = "depression_high_bmi", values_from = "n") |>
  rename(" " = "age")
subgroup_counts_depression_high_bmi_age_2 <- 
  subgroup_counts_depression_high_bmi_age_data_2 |>
  tableGrob(
    theme = ttheme_minimal(
      base_size = 7, 
      base_family = "sans",
      core = list(
        bg_params = list(fill = c(rep(c("#eff3f2", "white"), 2), "#eff3f2")),
        fg_params = list(
          hjust = 0, 
          x = unit(0, "cm"),
          fontface = c(rep("bold", 5), rep("plain", 20))
        )
      ),
      colhead = list(fg_params = list(hjust = 0, x = unit(0, "cm")))
    )
  ) |>
  gtable_remove_grobs("rowhead-fg") |>
  gtable_add_grob(linesGrob(y = c(0, 0), gp = gpar(lwd = 2, col = "black")), t = 1, l = 2, r = 6) |>
  gtable_squash_cols(1)

Figure_5_t <- cowplot::plot_grid(
  Figure_5,
  subgroup_counts_depression_high_bmi_age_2,
  labels = c("", "c"), 
  label_size = 8,
  ncol = 1,
  rel_heights = c(9, 4)
)

### Supplementary figure 2
SupplementaryFigure2 <- 
  variable_importance_data |>
  ggplot(aes(y = variable_name, x = variable_importance_num)) +
  geom_bar(stat = "identity", fill = "red3") +
  scale_x_continuous(breaks = seq(0, 0.6, 0.1)) +
  coord_cartesian(xlim = c(0, 0.6)) +
  xlab("variable importance") + ylab("") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        text = element_text(family = "sans", size = 11))

### Supplementary figure 3
SupplementaryFigure3a <-
  ggplot(cate_lines_data_sex_high_bmi_age) +
  geom_line(aes(y = 100 * estimate, x = age_grp, color = sex_high_bmi),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * `95% CI - lower`, 
                    ymax = 100 * `95% CI - upper`,
                    x = age_grp, 
                    color = sex_high_bmi),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * estimate, x = age_grp, color = sex_high_bmi),
             position = position_dodge(width = 2.5), size = 1) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.3) +
  ylab("risk difference (percentage points)") + xlab("age (yrs)") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-5, 15)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        text = element_text(family = "sans", size = 7.5))

SupplementaryFigure3b <- 
  ggplot(cate_histogram_data_sex_high_bmi_age$data) + 
  geom_histogram(aes(x = 100 * predictions,
                     y = after_stat(count) * 
                       5 / 
                       rep(
                         cate_histogram_data_sex_high_bmi_age$label$count, 
                         each = 2 * 87
                       ),
                     alpha = ate_diff_signif), 
                 binwidth = 0.2,
                 fill = "red3") +
  geom_vline(xintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.25) +
  geom_label(aes(x = x_l, y = y, label = lower), vjust = 0.5, nudge_y = 0,
             data = cate_histogram_data_sex_high_bmi_age$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  geom_label(aes(x = x_h, y = y, label = higher), vjust = 0.5, nudge_y = 0,
             data = cate_histogram_data_sex_high_bmi_age$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  facet_wrap(~ sex_high_bmi_age, nrow = 4) +
  coord_cartesian(xlim = c(-5, 15)) + 
  xlab("conditional risk difference (percentage points)") + ylab("density") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.52)) +
  theme(legend.position = "none",
        text = element_text(family = "sans", size = 7.5),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
        strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

SupplementaryFigure3_t <- cowplot::plot_grid(
  SupplementaryFigure3a, 
  SupplementaryFigure3b, 
  labels = c("a", "b"), 
  label_size = 8,
  rel_widths = c(1, 1),
  rel_heights = c(7, 9)
)

subgroup_counts_sex_high_bmi_age <- 
  cate_lines_data_sex_high_bmi_age |>
  select(sex_high_bmi_age, n) |>
  mutate(
    sex_high_bmi = case_when(
      str_detect(sex_high_bmi_age, "^male_no") ~ "male, not high BMI",
      str_detect(sex_high_bmi_age, "^female_no") ~ "female, not high BMI",
      str_detect(sex_high_bmi_age, "^male_yes") ~ "male, high BMI",
      str_detect(sex_high_bmi_age, "^female_yes") ~ "female, high BMI",
    ),
    age = paste(str_extract(sex_high_bmi_age, "\\d\\d-\\d\\d"), "yrs")
  ) |>
  select(-sex_high_bmi_age) |>
  pivot_wider(names_from = "sex_high_bmi", values_from = "n") |>
  bind_rows(
    cate_histogram_data_sex_high_bmi_age$data |>
      count(sex_high_bmi_age) |>
      mutate(
        sex_high_bmi = str_remove(sex_high_bmi_age, ", [<\u2265]50 yrs"),
        age = str_extract(sex_high_bmi_age, "[<\u2265]50 yrs")
      ) |>
      select(-sex_high_bmi_age) |>
      pivot_wider(names_from = "sex_high_bmi", values_from = "n")
  ) |>
  dplyr::select(age, starts_with("male"), starts_with("female")) |>
  rename(" " = "age") |>
  tableGrob(
    theme = ttheme_minimal(
      base_size = 7, 
      base_family = "sans",
      core = list(
        bg_params = list(fill = c(rep(c("#eff3f2", "white"), 3), "#eff3f2")),
        fg_params = list(
          hjust = 0, 
          x = unit(0, "cm"),
          fontface = c(rep("bold", 7), rep("plain", 28))
        )
      ),
      colhead = list(fg_params = list(hjust = 0, x = unit(0, "cm")))
    )
  ) |>
  gtable_remove_grobs("rowhead-fg") |>
  gtable_add_grob(linesGrob(y = c(0, 0), gp = gpar(lwd = 2, col = "black")), t = 1, l = 2, r = 6) |>
  gtable_squash_cols(1)

SupplementaryFigure3 <- cowplot::plot_grid(
  SupplementaryFigure3_t,
  subgroup_counts_sex_high_bmi_age,
  labels = c("", "c"), 
  label_size = 8,
  ncol = 1,
  rel_heights = c(9, 5.5)
)

### Supplementary figure 4
SupplementaryFigure4a <-
  ggplot(cate_lines_data_sex_depression_age) +
  geom_line(aes(y = 100 * estimate, x = age_grp, color = sex_depression),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * `95% CI - lower`, 
                    ymax = 100 * `95% CI - upper`,
                    x = age_grp, 
                    color = sex_depression),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * estimate, x = age_grp, color = sex_depression),
             position = position_dodge(width = 2.5), size = 1) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.3) +
  ylab("risk difference (percentage points)") + xlab("age (yrs)") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-5, 15)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        text = element_text(family = "sans", size = 7.5))

SupplementaryFigure4b <- 
  ggplot(cate_histogram_data_sex_depression_age$data) + 
  geom_histogram(aes(x = 100 * predictions,
                     y = after_stat(count) * 
                       5 / 
                       rep(
                         cate_histogram_data_sex_depression_age$label$count, 
                         each = 2 * 87
                       ),
                     alpha = ate_diff_signif), 
                 binwidth = 0.2,
                 fill = "red3") +
  geom_vline(xintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.25) +
  geom_label(aes(x = x_l, y = y, label = lower), vjust = 0.5, nudge_y = 0,
             data = cate_histogram_data_sex_depression_age$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  geom_label(aes(x = x_h, y = y, label = higher), vjust = 0.5, nudge_y = 0,
             data = cate_histogram_data_sex_depression_age$label,
             label.size = 0.12, size = 1.8,
             label.padding = unit(0.09, "lines"),
             label.r = unit(0.06, "lines")) +
  facet_wrap(~ sex_depression_age, nrow = 4) +
  coord_cartesian(xlim = c(-5, 15)) + 
  xlab("conditional risk difference (percentage points)") + ylab("density") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), limits = c(0, 0.48)) +
  theme(legend.position = "none",
        text = element_text(family = "sans", size = 7.5),
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
        strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm")))

SupplementaryFigure4_t <- cowplot::plot_grid(
  SupplementaryFigure4a, 
  SupplementaryFigure4b, 
  labels = c("a", "b"), 
  label_size = 8,
  rel_widths = c(1, 1),
  rel_heights = c(7, 9)
)

subgroup_counts_sex_depression_age <- 
  cate_lines_data_sex_depression_age |>
  select(sex_depression_age, n) |>
  mutate(
    sex_depression = case_when(
      str_detect(sex_depression_age, "^male_no") ~ "male, no depression",
      str_detect(sex_depression_age, "^female_no") ~ "female, no depression",
      str_detect(sex_depression_age, "^male_yes") ~ "male, depression",
      str_detect(sex_depression_age, "^female_yes") ~ "female, depression",
    ),
    age = paste(str_extract(sex_depression_age, "\\d\\d-\\d\\d"), "yrs")
  ) |>
  select(-sex_depression_age) |>
  pivot_wider(names_from = "sex_depression", values_from = "n") |>
  bind_rows(
    cate_histogram_data_sex_depression_age$data |>
      count(sex_depression_age) |>
      mutate(
        sex_depression = str_remove(sex_depression_age, ", [<\u2265]50 yrs"),
        age = str_extract(sex_depression_age, "[<\u2265]50 yrs")
      ) |>
      select(-sex_depression_age) |>
      pivot_wider(names_from = "sex_depression", values_from = "n")
  ) |>
  dplyr::select(age, starts_with("male"), starts_with("female")) |>
  rename(" " = "age") |>
  tableGrob(
    theme = ttheme_minimal(
      base_size = 7, 
      base_family = "sans",
      core = list(
        bg_params = list(fill = c(rep(c("#eff3f2", "white"), 3), "#eff3f2")),
        fg_params = list(
          hjust = 0, 
          x = unit(0, "cm"),
          fontface = c(rep("bold", 7), rep("plain", 28))
        )
      ),
      colhead = list(fg_params = list(hjust = 0, x = unit(0, "cm")))
    )
  ) |>
  gtable_remove_grobs("rowhead-fg") |>
  gtable_add_grob(linesGrob(y = c(0, 0), gp = gpar(lwd = 2, col = "black")), t = 1, l = 2, r = 6) |>
  gtable_squash_cols(1)

SupplementaryFigure4 <- cowplot::plot_grid(
  SupplementaryFigure4_t,
  subgroup_counts_sex_depression_age,
  labels = c("", "c"), 
  label_size = 8,
  ncol = 1,
  rel_heights = c(9, 5.5)
)

### Supplementary Figure 5
SupplementaryFigure5 <- 
  ggplot(health_conditions_plot_data) +
  geom_line(aes(y = 100 * e1, x = a, color = c),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * l1, 
                    ymax = 100 * u1,
                    x = a, 
                    color = c),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * e1, x = a, color = c),
             position = position_dodge(width = 2.5), size = 1.2) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.4) +
  ylab("risk difference (percentage points)") + xlab("") +
  ggtitle("chronic asthma") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-10, 20)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()) +
  ggplot(health_conditions_plot_data) +
  geom_line(aes(y = 100 * e2, x = a, color = c),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * l2, 
                    ymax = 100 * u2,
                    x = a, 
                    color = c),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * e2, x = a, color = c),
             position = position_dodge(width = 2.5), size = 1.2) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.4) +
  ylab("") + xlab("") + ggtitle("chronic or frequent headaches") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-10, 20)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()) +
  ggplot(health_conditions_plot_data) +
  geom_line(aes(y = 100 * e3, x = a, color = c),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * l3, 
                    ymax = 100 * u3,
                    x = a, 
                    color = c),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * e3, x = a, color = c),
             position = position_dodge(width = 2.5), size = 1.2) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.4) +
  ylab("risk difference (percentage points)") + xlab("age (yrs)") +
  ggtitle("PTSD") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-10, 20)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()) +
  ggplot(health_conditions_plot_data) +
  geom_line(aes(y = 100 * e4, x = a, color = c),
            position = position_dodge(width = 2.5), linewidth = 0.7) +
  geom_errorbar(aes(ymin = 100 * l4, 
                    ymax = 100 * u4,
                    x = a, 
                    color = c),
                position = position_dodge(width = 2.5),
                linewidth = 0.7,
                width = 0) + 
  geom_point(aes(y = 100 * e4, x = a, color = c),
             position = position_dodge(width = 2.5), size = 1.2) +
  geom_hline(yintercept = 100 * ate_table$estimate[1], linetype = 2, linewidth = 0.4) +
  ylab("") + xlab("age (yrs)") + ggtitle("diabetes") +
  coord_cartesian(xlim = c(14, 65), ylim = c(-10, 20)) +
  scale_color_jama() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        text = element_text(family = "sans", size = 7.5))

### Supplementary figure 6
SupplementaryFigure6_split <- ggplot(
  overlap_data
) + 
  geom_density(
    aes(x = W_hat, y = after_stat(density), fill = W_orig),
    alpha = 0.5
  ) +
  coord_cartesian(xlim = c(0, 1)) + 
  xlab("propensity score") + ylab("density") +
  labs(fill = "RT-PCR test:") +
  scale_y_continuous(breaks = 1:5) +
  scale_fill_jama() + 
  theme(
    text = element_text(family = "sans", size = 7.5)
  )
SupplementaryFigure6_split_wei <- ggplot(
  overlap_data
) + 
  geom_density(
    aes(x = W_hat, y = after_stat(density), fill = W_orig, weight = IPW),
    alpha = 0.5
  ) +
  coord_cartesian(xlim = c(0, 1)) + 
  xlab("propensity score") + ylab("density") +
  labs(fill = "RT-PCR test:") +
  scale_y_continuous(breaks = 1:5) +
  scale_fill_jama() + 
  theme(
    text = element_text(family = "sans", size = 7.5)
  )

# use ggplot_build to modify density to be negative
# add density of weighted pseudo population
SupplementaryFigure6_split_build <- ggplot_build(SupplementaryFigure6_split)
SupplementaryFigure6_split_wei_build <- ggplot_build(SupplementaryFigure6_split_wei)
SupplementaryFigure6_negative_orig <- SupplementaryFigure6_split_build$data[[1]] |>
  select(fill, y , x) |>
  filter(fill == "#374E55FF")
SupplementaryFigure6_positive_orig <- SupplementaryFigure6_split_build$data[[1]] |>
  select(fill, y , x) |>
  filter(fill != "#374E55FF")
SupplementaryFigure6_orig <- SupplementaryFigure6_split_build$data[[1]] |>
  as_tibble() |>
  select(fill, y , x) |>
  mutate(test = ifelse(fill == "#374E55FF", "Negative", "Positive"))
SupplementaryFigure6_negative_wei <- SupplementaryFigure6_split_wei_build$data[[1]] |>
  select(fill, y , x) |>
  filter(fill == "#374E55FF") |>
  mutate(y = -y)
SupplementaryFigure6_positive_wei <- SupplementaryFigure6_split_wei_build$data[[1]] |>
  select(fill, y , x) |>
  filter(fill != "#374E55FF") |>
  mutate(y = -y)
SupplementaryFigure6_wei <- SupplementaryFigure6_split_wei_build$data[[1]] |>
  as_tibble() |>
  select(fill, y , x) |>
  mutate(y = -y, test = ifelse(fill == "#374E55FF", "Negative", "Positive"))

SupplementaryFigure6 <- ggplot(SupplementaryFigure6_orig, aes(x = x, y = y)) +
  geom_area(aes(fill = test), position = "identity", alpha = 0.5) +
  geom_area(aes(fill = test), data = SupplementaryFigure6_wei, position = "identity", alpha = 0.5) +
  geom_segment(
    aes(x = x, xend = xend, y = y, yend = yend),
    data = tibble(
      x = min(SupplementaryFigure6_orig$x), xend = max(SupplementaryFigure6_orig$x),
      y = 0, yend = 0
    )
  ) +
  annotate("label", 0.4, -1, label = "weighted", size = 3.2, family = "sans") +
  annotate("label", 0.4, 1, label = "unweighted", size = 3.2, family = "sans") +
  coord_cartesian(xlim = c(0, 1)) + 
  xlab("propensity score") + ylab("density") +
  labs(fill = "RT-PCR test:") +
  scale_y_continuous(
    breaks = seq(-6, 6, 2), 
    labels = as.character(c(6, 4, 2, 0, 2, 4, 6))
  ) +
  scale_fill_jama() +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        text = element_text(family = "sans", size = 11))

### Supplementary Figure 7
SupplementaryFigure7 <- love_data |>
  ggplot(aes(x = value, y = covariate_name, colour = Cohort)) +
  geom_point(size = 2) +
  geom_line(aes(group = Cohort), orientation = "y") +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.4) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.4) +
  ggplot2::xlab("Absolute standardized mean difference") +
  ggplot2::ylab("") +
  scale_color_jama() + 
  ggplot2::scale_x_continuous(
    breaks = c(0, 0.1),
    limits = c(0, 0.12)
  ) +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        text = element_text(family = "sans", size = 11))

# Supplementary figure 8-10
X_type <- ecdf_data |>
  dplyr::summarise(
    dplyr::across(
      1:15, \(x) ifelse(is.factor(x), "dis", "con")
    )
  )
ec_x_scale_width <- c(10, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1)
# generate ecdf plots
ec_plots <- purrr::pmap(
  list(
    names = names(ecdf_data[1:15]),
    ec_x_scale_width = ec_x_scale_width,
    X_type = X_type
  ),
  \(names, ec_x_scale_width, X_type) {
    plot_data <- ecdf_data |>
      select("RT-PCR test result:",
             "IPW",
             "covariate_values" = dplyr::all_of(names)
      ) |>
      mutate("covariate_name" = names) |>
      group_by(`RT-PCR test result:`) |>
      arrange(`covariate_values`) |>
      mutate(
        cum_pct_ori = seq_len(n()) / n(),
        cum_pct_wei = cumsum(IPW) / sum(IPW)
      ) |>
      ungroup()
    ## plot eCDF
    if (X_type == "dis") {
      plot_data <- plot_data |>
        mutate(
          covariate_values_fct = covariate_values,
          covariate_values = as.numeric(covariate_values)
        )
      covariate_levels <- levels(plot_data[["covariate_values_fct"]])
    }
    covariate_values <- plot_data[["covariate_values"]]
    p_ori <- plot_data |>
      ggplot() +
      geom_step(
        aes(
          x = covariate_values,
          y = cum_pct_ori,
          color = `RT-PCR test result:`
        ),
        linewidth = 0.7
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, 1/3), 
        labels = c("0.00", "0.33", "0.66", "1.00")
      ) +
      xlab("") + ylab(names) +
      scale_color_jama() +
      theme(
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        text = element_text(family = "sans", size = 8)
      )
    
    p_wei <- ggplot(plot_data) +
      geom_step(
        aes(
          x = covariate_values,
          y = cum_pct_wei,
          color = `RT-PCR test result:`
        ),
        linewidth = 0.7
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, 1/3), 
        labels = c("0.00", "0.33", "0.66", "1.00")
      ) +
      xlab("") +
      ylab("") +
      scale_color_jama() +
      theme(
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        text = element_text(family = "sans", size = 8)
      )
    
    if (X_type == "con") {
      p_ori <- p_ori +
        scale_x_continuous(
          breaks = seq(
            floor_dec(
              min(covariate_values),
              digits = decimalplaces(ec_x_scale_width)
            ),
            ceiling_dec(
              max(covariate_values),
              digits = decimalplaces(ec_x_scale_width)
            ),
            ec_x_scale_width
          )
        )
      p_wei <- p_wei +
        scale_x_continuous(
          breaks = seq(
            floor_dec(
              min(covariate_values),
              digits = decimalplaces(ec_x_scale_width)
            ),
            ceiling_dec(
              max(covariate_values),
              digits = decimalplaces(ec_x_scale_width)
            ),
            ec_x_scale_width
          )
        )
    } else {
      p_ori <- p_ori +
        scale_x_continuous(
          breaks = seq_along(covariate_levels),
          labels = covariate_levels
        )
      p_wei <- p_wei +
        scale_x_continuous(
          breaks = seq_along(covariate_levels),
          labels = covariate_levels
        )
    }
    return(list(ori = p_ori, wei = p_wei))
  }
)
ec_plots_unadjusted <- purrr::list_transpose(ec_plots)[[1]]
ec_plots_adjusted <- purrr::list_transpose(ec_plots)[[2]]
# produce plot
SupplementaryFigure8_plots <- 
  SupplementaryFigure9_plots <- 
  SupplementaryFigure10_plots <- 
  vector("list", 10)
SupplementaryFigure8_plots[seq(1, 9, 2)] <- ec_plots_unadjusted[1:5]
SupplementaryFigure9_plots[seq(1, 9, 2)] <- ec_plots_unadjusted[6:10]
SupplementaryFigure10_plots[seq(1, 9, 2)] <- ec_plots_unadjusted[11:15]
SupplementaryFigure8_plots[seq(2, 10, 2)] <- ec_plots_adjusted[1:5]
SupplementaryFigure9_plots[seq(2, 10, 2)] <- ec_plots_adjusted[6:10]
SupplementaryFigure10_plots[seq(2, 10, 2)] <- ec_plots_adjusted[11:15]
SupplementaryFigure8 <- cowplot::ggdraw(
  gridExtra::arrangeGrob(
    patchwork::patchworkGrob(
      (
        Reduce("+", SupplementaryFigure8_plots) +
          patchwork::plot_layout(
            nrow = 5,
            ncol = 2,
            guides = "collect"
          ) &
          theme(legend.position = "top")
      )
    ),
    left = textGrob(
      label = "cumulative probability", 
      rot = 90, 
      gp = gpar(fontsize = 10, fontfamily = "sans")
    ), 
    bottom = textGrob(
      label = "covariate values", 
      rot = 0, 
      gp = gpar(fontsize = 10, fontfamily = "sans")
    )
  )
) +
  theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA)
  )
SupplementaryFigure9 <- cowplot::ggdraw(
  gridExtra::arrangeGrob(
    patchwork::patchworkGrob(
      (
        Reduce("+", SupplementaryFigure9_plots) +
          patchwork::plot_layout(
            nrow = 5,
            ncol = 2,
            guides = "collect"
          ) &
          ggplot2::theme(legend.position = "top")
      )
    ),
    left = textGrob(
      label = "cumulative probability", 
      rot = 90, 
      gp = gpar(fontsize = 10, fontfamily = "sans")
    ), 
    bottom = textGrob(
      label = "covariate values", 
      rot = 0, 
      gp = gpar(fontsize = 10, fontfamily = "sans")
    )
  )
) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA)
  )
SupplementaryFigure10 <- cowplot::ggdraw(
  gridExtra::arrangeGrob(
    patchwork::patchworkGrob(
      (
        Reduce("+", SupplementaryFigure10_plots) +
          patchwork::plot_layout(
            nrow = 5,
            ncol = 2,
            guides = "collect"
          ) &
          ggplot2::theme(legend.position = "top")
      )
    ),
    left = textGrob(
      label = "cumulative probability", 
      rot = 90, 
      gp = gpar(fontsize = 10, fontfamily = "sans")
    ), 
    bottom = textGrob(
      label = "covariate values", 
      rot = 0, 
      gp = gpar(fontsize = 10, fontfamily = "sans")
    )
  )
) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 12)
  )

### Supplementary figure 11
SupplementaryFigure11 <- pcr_sens_importance_rank_plot_data |>
  ggplot(aes(x = rank, y = n)) +
  geom_col(fill = "red3") +
  geom_vline(xintercept = 5.5, linetype = 2, linewidth = 0.4) +
  facet_wrap(~ variable_name, nrow = 5) +
  ylab("count") +
  theme(
    text = element_text(family = "sans", size = 7.5),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
    strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  )

### Supplementary figure 12
SupplementaryFigure12 <- tun_sens_importance_rank_plot_data |>
  ggplot(aes(x = rank, y = n)) +
  geom_col(fill = "red3") +
  geom_vline(xintercept = 5.5, linetype = 2, linewidth = 0.4) +
  facet_wrap(~ variable_name, nrow = 5) +
  ylab("count") +
  theme(
    text = element_text(family = "sans", size = 7.5),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
    strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
  )

# === plot data   ==============================================================
# figure 2
readr::write_excel_csv2(
  tibble(`conditional risk difference` = cate_predictions), 
  file = "../tables/figure2b_data.csv"
)

# figure 3
figure_3a_data <- cate_lines_data_depression_high_bmi_age |>
  select(depression_high_bmi_age, estimate, `95% CI - lower`, `95% CI - upper`) |>
  rename("Risk difference" = "estimate")
figure_3b_data <- cate_histogram_data_depression_high_bmi_age

readr::write_excel_csv2(
  figure_3a_data, 
  file = "../tables/figure3a_data.csv"
)
readr::write_excel_csv2(
  figure_3b_data, 
  file = "../tables/figure3b_data.csv"
)

# figure 4
figure_4a_data <- cate_lines_data_depression_high_bmi_sex |>
  select(depression_high_bmi_sex, estimate, `95% CI - lower`, `95% CI - upper`) |>
  rename("Risk difference" = "estimate")
figure_4b_data <- cate_histogram_data_depression_high_bmi_sex

readr::write_excel_csv2(
  figure_4a_data, 
  file = "../tables/figure4a_data.csv"
)
readr::write_excel_csv2(
  figure_4b_data, 
  file = "../tables/figure4b_data.csv"
)

# figure 5
readr::write_excel_csv2(
  figure_5a_data, 
  file = "../tables/figure5a_data.csv"
)
readr::write_excel_csv2(
  figure_5b_data, 
  file = "../tables/figure5b_data.csv"
)

# ===   save   =================================================================
### Figure 2
agg_jpeg(filename = "../figures/Figure2.jpg",
         width = 18, height = 8, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_2
dev.off()

### Figure 3
agg_jpeg(filename = "../figures/Figure3.jpg",
         width = 18, height = 9, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_3
dev.off()

agg_jpeg(filename = "../figures/Figure3_withtable.jpg",
         width = 18, height = 14.5, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_3_t
dev.off()

agg_jpeg(filename = "../figures/FeaturedImage.jpg",
         width = 1200, height = 675, units = "px", quality = 100,
         scaling = 2.3)
Figure_3_t
dev.off()

write_xlsx(
  x = subgroup_counts_depression_high_bmi_age_data,
  path = "../tables/table_for_Figure3.xlsx"
)

### Figure 4
agg_jpeg(filename = "../figures/Figure4.jpg",
         width = 18, height = 9, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_4
dev.off()

agg_jpeg(filename = "../figures/Figure4_t.jpg",
         width = 18, height = 11.2, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_4_t
dev.off()

write_xlsx(
  x = subgroup_counts_depression_high_bmi_sex_data,
  path = "../tables/table_for_Figure4.xlsx"
)

### Figure 5
agg_jpeg(filename = "../figures/Figure5.jpg",
         width = 18, height = 9, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_5
dev.off()

agg_jpeg(filename = "../figures/Figure5_t.jpg",
         width = 18, height = 13, units = "cm", res = 320, quality = 100,
         scaling = 1)
Figure_5_t
dev.off()

write_xlsx(
  x = subgroup_counts_depression_high_bmi_age_data_2,
  path = "../tables/table_for_Figure5.xlsx"
)


### Supplementary Figure 1
cat(export_svg(tree_plot), file = "../figures/SupplementaryFigure1.svg")

### Supplementary Figure 2
agg_jpeg(filename = "../figures/SupplementaryFigure2.jpg",
         width = 18, height = 12, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure2
dev.off()

### Supplementary Figure 3
agg_jpeg(filename = "../figures/SupplementaryFigure3.jpg",
         width = 18, height = 14.5, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure3
dev.off()

### Supplementary Figure 4
agg_jpeg(filename = "../figures/SupplementaryFigure4.jpg",
         width = 18, height = 14.5, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure4
dev.off()

### Supplementary Figure 5
agg_jpeg(filename = "../figures/SupplementaryFigure5.jpg",
         width = 18, height = 18, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure5
dev.off()

### Supplementary Figure 6
agg_jpeg(filename = "../figures/SupplementaryFigure6.jpg",
         width = 18, height = 9, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure6
dev.off()

### Supplementary Figure 7
agg_jpeg(filename = "../figures/SupplementaryFigure7.jpg",
         width = 18, height = 12, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure7
dev.off()

### Supplementary Figure 8
agg_jpeg(filename = "../figures/SupplementaryFigure8.jpg",
         width = 18, height = 20.5, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure8
dev.off()

### Supplementary Figure 9
agg_jpeg(filename = "../figures/SupplementaryFigure9.jpg",
         width = 18, height = 20.5, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure9
dev.off()

### Supplementary Figure 10
agg_jpeg(filename = "../figures/SupplementaryFigure10.jpg",
         width = 18, height = 20.5, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure10
dev.off()

### Supplementary Figure 11
agg_jpeg(filename = "../figures/SupplementaryFigure11.jpg",
         width = 18, height = 18, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure11
dev.off()

### Supplementary Figure 12
agg_jpeg(filename = "../figures/SupplementaryFigure12.jpg",
         width = 18, height = 18, units = "cm", res = 320, quality = 100,
         scaling = 1)
SupplementaryFigure12
dev.off()
