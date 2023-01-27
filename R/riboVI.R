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

ribovi <- function(counts, S, X, Z, calc_S = FALSE) {
    
}