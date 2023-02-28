#' ncvi: a package for fitting Bayesian GLMM to RNA-Seq and Ribo-Seq data
#' using non-conjugate variational methods.
#'
#' The package provides functions to fit a Bayesian GLMM using
#' a variational algorithm.
#'
#'@section Prep functions
#'
#' @docType package
#' @name riboVI 
#'
#' @import dplyr
NULL

#' Run `riboVI` algorithm
#' @inheritParams prep_ribovi_data
#' @param options List of options settings for `fit_ribovi()` 
#' @export

ribovi <- function(counts, S, X, Z, normalize = TRUE, calc_S = FALSE,
                   priors, init, 
                   options = list(max_iter = 400,
                                  elbo_delta = 0.01,
                                  verbose = T,
                                  fixed_iter = F,
                                  short_out = TRUE)) {
    message("Preparing data...")
    ribovi_data <- prep_ribovi_data(counts = counts, X = X, Z = Z, 
                                    S = S, calc_S = calc_S,
                                    normalize = normalize)
    message("Running `riboVI`...")
    fit <- fit_ribovi(data = ribovi_data$data,
                      init = ribovi_data$init,
                      priors = ribovi_data$priors,
                      options = options)
}

