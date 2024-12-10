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
