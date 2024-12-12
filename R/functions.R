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



#' A transformation recipe to pre-process the data
#'
#' @param data the lipidomics dataset
#' @param metabolite_variable the column of the metabolites
#'
#' @return recipe object
create_recipe_spec <- function(data, metabolite_variable) {
  recipes::recipe(data) %>%
    recipes::update_role(
      {{ metabolite_variable }},
      age,
      gender,
      new_role = "predictor"
    ) %>%
    recipes::update_role(class,
      new_role = "outcome"
    ) %>%
    recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}


#' Create a workflow object of the model and transformations
#'
#' @param model_specs The model specs
#' @param recipe_specs The recipe specs
#'
#' @return A workflow object
create_model_workflow <- function(model_specs, recipe_specs) {
  workflows::workflow() %>%
    workflows::add_model(model_specs) %>%
    workflows::add_recipe(recipe_specs)
}


#' Create a tidy output of model results
#'
#' @param workflow_fitted_model the model workflow object fitted
#'
#' @return A tidy dataframe
tidy_model_output <- function(workflow_fitted_model) {
  workflow_fitted_model %>%
    workflows::extract_fit_parsnip() %>%
    broom::tidy(exponentiate = TRUE)
}


#' Conver long form data to list of wide form data by metabolites
#'
#' @param data The lipidomics data set
#'
#' @return A list of data frames
split_by_metabolite <- function(data) {
  data %>%
    column_values_to_snake_case(metabolite) %>%
    dplyr::group_split(metabolite) %>%
    purrr::map(metabolites_to_wider)
}



#' Generate the results of a model
#'
#' @param data The lipidomics dataset
#'
#' @return A data frame
generate_model_results <- function(data) {
  create_model_workflow(
    parsnip::logistic_reg() %>%
      parsnip::set_engine("glm"),
    data %>% create_recipe_spec(tidyselect::starts_with("metabolite_"))
  ) %>%
    parsnip::fit(data) %>%
    tidy_model_output()
}



#' Model results for all metabolites
#'
#' @param data The lipidomics dataset
#'
#' @return A data frame
calculate_estimates <- function(data) {
    data %>%
        split_by_metabolite() %>%
        purrr::map(generate_model_results) %>%
        purrr::list_rbind() %>%
        dplyr::filter(stringr::str_detect(term, "metabolite_"))
}

