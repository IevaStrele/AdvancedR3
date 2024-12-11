#' Summarize metabolites as means and SDs
#'
#' @param lipidomics data
#'
#' @return A data.frame/tibble

descriptive_stats <- function(data) {
  stats <- data %>%
    dplyr::group_by(metabolite) %>%
    dplyr::summarise(dplyr::across(
      value,
      list(
        mean = mean,
        sd = sd,
        iqr = IQR
      )
    )) %>%
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric),
      ~ round(.x, digits = 1)
    ))
  return(stats)
}



#' Distribution plots of metabolites
#'
#' @param data lipidomics data
#'
#' @return A plot object.

plot_distributions <- function(data) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free")
  return(plot)
}



#' Convert a column's char values to snake case
#'
#' @param data The lipidomics dataset
#' @param columns the column to convert
#'
#' @return A data frame
column_values_to_snake_case <- function(data, columns) {
  data %>%
    dplyr::mutate(dplyr::across(
      {{ columns }},
      snakecase::to_snake_case
    ))
}


#' Metabolite column to wider format
#'
#' @param data The lipidomics dataset
#'
#' @return A data frame

metabolites_to_wider <- function(data) {
    data %>%
        tidyr::pivot_wider(
            names_from = metabolite,
            values_from = value,
            values_fn = mean,
            names_prefix = "metabolite_"
        )
}

