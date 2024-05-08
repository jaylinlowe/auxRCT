

#' Mean Imputation
#' @description
#' Replaces all missing values with the mean of the column and creates new
#' columns to capture what values were initially missing.
#'
#'
#' @param df data frame to modify
#'
#' @return `mean_imputation` returns the original data frame, modified so that all missing values are
#' have been replaced with the mean of the column. All columns with missing values
#' now have an additional corresponding column with `_mis` tagged onto
#' the original column name. The new columns will have a 1 in any row where there
#' was originally a missing value in the corresponding column and a 0 otherwise.
#' @export
#'
mean_imputation <- function(df) {
  orig_names <- colnames(df) 
  df_new <- df %>%
    mutate(across(where(is.numeric), ~case_when(is.na(.x) ~ 1,
                                                TRUE ~ 0),
                  .names = "{.col}_mis")) %>%
    mutate(across(where(is.numeric), ~case_when(is.na(.x) ~ mean(.x, na.rm = TRUE),
                                                TRUE ~ .x)))

  #get table count of NAs per column, then pull any names with 0 and remove those columns
  na_count_df <- df_new %>% summarise(across(ends_with("_mis"), ~ sum(is.na(.))))
  cols_to_remove <- colnames(na_count_df)[which(na_count_df == 0)]
  final_cols_to_remove <- setdiff(cols_to_remove, orig_names) #solves problem of original columns ending in _mis
  
  
  df_new <- df_new %>%
    select(-all_of(final_cols_to_remove))

  return(df_new)
}
