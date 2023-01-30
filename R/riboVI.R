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
#' @export

ribovi <- function(counts, S = NULL, X, Z, calc_S = FALSE,
                   normalize = TRUE) {
    ribovi_data <- prep_ribovi_data(counts = counts, X = X, Z = Z, 
                                    S = S, calc_S = calc_S,
                                    normalize = normalize)
}