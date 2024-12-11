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
        sd = sd
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
#' @param lipidomics data
#'
#' @return A plot object.

plot_distributions <- function(data) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram() +
        ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free")
    return(plot)
}
