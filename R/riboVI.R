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
                   priors , init, 
                   options = list(max_iter = 400,
                                  elbo_delta = 0.01,
                                  verbose = T,
                                  fixed_iter = F,
                                  short_out = TRUE)) {
    message("Preparing data...")
    ribovi_data <- prep_ribovi_data(counts = counts, X = X, Z = Z, 
                                    S = S, calc_S = calc_S,
                                    normalize = normalize,
                                    priors = priors, init = init)
    message("Running `riboVI`...")
    fit <- fit_ribovi(data = ribovi_data$data,
                      init = ribovi_data$init,
                      priors = ribovi_data$priors,
                      options = options)
}

y5 <- edgeR::DGEList(counts = counts[, 1:18], group = group[1:18])
y5 <- edgeR::calcNormFactors(y5)
libs5 <- y5$samples$lib.size
norms5 <- y5$samples$norm.factors

y8 <- edgeR::DGEList(counts = counts[, 19:36], group = group[19:36])
y8 <- edgeR::calcNormFactors(y8)
libs8 <- y8$samples$lib.size
norms8 <- y8$samples$norm.factors

S5 <- log(libs5 * norms5) - median(log(libs5 * norms5))
S8 <- log(libs8 * norms8) - median(log(libs8 * norms8))
S <- c(S5, S8)
libs <- c(libs5, libs8)
norms <- c(norms5, norms8)


normed_counts <- round(t(t(real_data[, 2:37]) / exp(S)))

