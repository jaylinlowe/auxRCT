

#' Calculating the mode for a specified column of a data frame
#'
#' @param df a data frame to analyze
#' @param col the name of a column in `df` to calculate the mode for
#'
#' @return The value that appears most frequently in the `var` column of `df`
#' @export
#' @noRd
mode_calc <- function(df, col) {
  index <- which(colnames(df) == col)
  v <- df[, index]
  uniqv <- unique(v)
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

#' Get vector of variable names
#'
#' @description
#' An easy way to get the correct variable names depending on the Shiny inputs.
#' If "All Variables" is the input, all column names are returned. The function
#' also returns an error if `variables` is a list including variable names and
#' "All Variables".
#'
#'
#' @param variables user specified input from Shiny function of what variables to select
#' @param X data frame of covariates
#'
#' @return A vector of variable names
#' @export
#' @noRd
get_var_names <- function(variables, X) {
  if (length(variables) == 1) {
    if (variables == "All Variables") {
      variables <- colnames(X)
    }
  }

  if (length(variables) > 1) {
    if ("All Variables" %in% variables) {
      stop("If 'All Variables' is selected, no other columns should be selected")
    }
  }
  return(variables)
}

#' Create plot of correlations colored by strength of relationship
#'
#' @param df data frame of numeric variables used for correlations
#' @param variables vector of column names in df to include in plot
#'
#' @return A ggplot object with printed numerical correlations and a manual color scheme
#' @export
#' @noRd
correlation_plot <- function(df, variables) {
  p <- df %>%
    select(all_of(variables)) %>%
    select_if(is.numeric) %>%
    cor() %>%
    reshape2::melt() %>%
    ggplot(aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 2))) +
    labs(x = "", y = "", title = "Correlations for Numerical Variables") +
    scale_fill_gradient2(
      low = muted("orange"),
      mid = "white",
      high = muted("navyblue"),
      midpoint = 0,
      limits = c(-1, 1)) +
    scale_x_discrete(limits = rev)
  return(p)
}

#' Split table into two tables of equal size to be displayed side by side
#'
#' @param df data frame to split
#'
#' @return A list of 2 data frames with half of the rows each, with a column
#' labeled "Importance Ranking" capturing the original order of `df`
#' @export
#' @noRd
importance_split <- function(df) {
  full_table <- df %>%
    mutate(`Importance Ranking` = seq(1:nrow(df))) %>%
    select(`Importance Ranking`, Variable, Importance)
  half_size <- round(nrow(df)/2)
  table1 <- full_table[1:half_size,]
  table2 <- full_table[(1+half_size):nrow(full_table), ]
  return(list(table1, table2))
}

#' Create plot to compare subgroup to rest of population for a specific variable
#'
#' @param Y vector of responses
#' @param X data frame of covariates
#' @param subgroup_def vector with "Subgroup"/"General Population" to define observations
#' @param var variable to show in plot
#' @param out name of outcome variable
#'
#' @return A ggplot object of a boxplot for numeric variables or barplot for
#' categorical variables to compare the subgroup to the rest of the population.
#' @export
#' @noRd
subgroup_plot <- function(Y, X, subgroup_def, var, out) {
  df <- cbind(Y, subgroup_def, X)
  extra_subgroup <- filter(df, `subgroup_def` == "Subgroup")
  extra_subgroup$`subgroup_def` <- rep("General Population", nrow(extra_subgroup))
  df <- rbind(df, extra_subgroup)

  if (var == out) {
    var <- "Y"
    df[[var]] <- as.numeric(df[[var]])
  }
  if (class(df[[var]]) == "numeric") {
    p <- ggplot(df, aes(x = .data[[var]], fill = `subgroup_def`)) +
      geom_boxplot() +
      labs(fill = "") +
      theme(legend.position = "bottom") +
      scale_fill_manual(values = c(rgb(95, 92, 163, max = 255), rgb(36, 45, 49, max = 255)))
  }
  else {
    p <- df %>%
      group_by(`subgroup_def`, .data[[var]]) %>%
      summarize(n = n()) %>%
      group_by(`subgroup_def`) %>%
      mutate(prop = n/sum(n)) %>%
      ggplot(aes(y = prop, fill = `subgroup_def`, x = .data[[var]])) +
      geom_bar(position = "dodge", stat = "identity") +
      labs(fill = "", y = "Proportion of Group") +
      coord_flip() +
      theme(legend.position = "bottom") +
      scale_fill_manual(values = c(rgb(95, 92, 163, max = 255), rgb(36, 45, 49, max = 255)))
  }
  return(p)
}

#' Create plot of number of NA values for subgroup versus population across all variables
#'
#' @param df data frame of data
#' @param subgroup_def vector with "Subgroup"/"General Population" to define observations
#'
#' @return A ggplot density plot comparing the number of NAs for each observation
#' in the subgroup to the rest of the population
#' @export
#'
#' @noRd
na_plot <- function(df, subgroup_def) {
  df <- cbind(df, subgroup_def)
  extra_subgroup <- filter(df, subgroup_def == "Subgroup")
  extra_subgroup$`subgroup_def` <- rep("General Population", nrow(extra_subgroup))
  df <- rbind(df, extra_subgroup)

  na_count <- rowSums(is.na(df))
  if (sum(na_count) > 0) {
    na_df <- data.frame(cbind(na_count, subgroup_def))
    na_df$na_count <- as.numeric(na_df$na_count)
    colnames(na_df) <- c("na_count", "Group Definition")
    p <- ggplot(na_df, aes(x = na_count, fill = `Group Definition`)) +
      geom_density(stat = "count", alpha = 0.5) +
      scale_fill_manual(values = c(rgb(95, 92, 163, max = 255), rgb(36, 45, 49, max = 255))) +
      labs(fill = "", x = "Number of NA Values", y = "Count") +
      theme(legend.position = "bottom") +
      ggtitle("Distribution of Number of NA Values per Observation")
  }
  else {p <- NULL}
  return(p)
}
