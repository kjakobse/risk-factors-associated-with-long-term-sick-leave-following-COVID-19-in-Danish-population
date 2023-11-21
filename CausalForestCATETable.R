CausalForestCATETable <- function (cf, 
                                   cov_list,
                                   subset = rep(TRUE, length(cf[["Y.orig"]])),
                                   level = 0.95) {
  if (!(is.list(subset) && length(subset) == 1 && is.logical(subset[[1]]) | 
       is.logical(subset) && length(subset) == length(cf[["Y.orig"]]))) {
    stop(
      glue::glue(
        "subset must be a named one element list with a logical vector of ",
        "length {length(cf[['Y.orig']])}, or a logical vector of ",
        "length {length(cf[['Y.orig']])}."
      )
    )
  }
  str_helper <- function (...) {
    string <- ""
    for(i in seq_along(cov_list)) {
      if (is.null(names(cov_list[[i]]))) {
        name <- str_remove(str_remove(...elt(i+1), names(cov_list)[i]), "^_")
      } else {
        name <- names(cov_list[[i]])[...elt(1)]
      }
      if(i == 1) {
        string <- paste0(string, name)
      } else {
        string <- paste0(string, "_", name)
      }
    }
    if(!is.null(subset_name)) {
      string <- paste0(string, "_", subset_name)
    }
    return(string)
  }
  subset_helper <- function (...) {
    for (i in seq_along(cov_list)) {
      name <- names(cov_list)[i]
      if (i == 1) {
        if (is.character(...elt(2))) {
          sub <- cf[["X.orig"]][[...elt(2)]] == 1
        } else {
          sub <- cf[["X.orig"]][[name]] %in% ...elt(2)
        }
      } else {
        if (is.character(...elt(i+1))) {
          sub <- sub & cf[["X.orig"]][[...elt(i+1)]] == 1
        } else {
          sub <- sub & (cf[["X.orig"]][[name]] %in% ...elt(i+1))
        }
      }
    }
    return(sub)
  }
  subset_name <- names(subset)
  if(is.list(subset)) subset <- subset[[1]]
  ci_names <- c(
    paste0(100 * level, "% CI - lower"),
    paste0(100 * level, "% CI - upper")
  )
  ci_names_pct <- c(
    paste0(100 * level, "% CI - lower (%)"),
    paste0(100 * level, "% CI - upper (%)")
  )
  cate_table <- pmap(
    c(
      list(id = seq_along(cov_list[[1]])),
      cov_list
    ),
    \(...) {
      string <- str_helper(...)
      sub <- subset_helper(...)
      SubgroupATETable(
        NULL, 
        c(paste0(names(cov_list), collapse = "_"), string), 
        cf, 
        level, 
        sub & subset
      )
    }
  ) |>
    list_rbind() |>
    mutate(
      `estimate (%)` = sprintf("%.1f", 100 * estimate),
      !!sym(ci_names_pct[1]) := sprintf("%.1f", 100 * !!sym(ci_names[1])),
      !!sym(ci_names_pct[2]) := sprintf("%.1f", 100 * !!sym(ci_names[2]))
    )
  return(cate_table)
}