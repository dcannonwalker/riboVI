## code to prepare `ribo_example` dataset goes here

geneid <- paste0("gene", 1:500)
devtools::load_all('../ncvi/')
N <- 8
settings = list(
    P = 4, # number fixed effects / treatments
    N = N, # number of samples / observations per group or gene
    U = N / 2, # number of random effects
    G = 500, # number of groups or genes
    X = cbind(1, rep(c(0, 1), each = N / 2),
              rep(rep(c(0, 1), each = N / 4), 2),
              c(rep(0, 3 * N / 4), rep(1, N / 4))), # f.e. design
    Z = rbind(diag(N / 2), diag(N / 2)), # r.e. design
    
    ## known parameters
    precision_mu0 = 1, # prior precision for mu0
    precision_u = 3, # prior precision for u
    priors = list(a_beta = c(1, 4, 4, 20), b_beta = c(3, 2, 2, 2),
                  a_u = 10, b_u = 2, pi0 = 0.8),
    seed = 33,
    mu0 = c(5, 0, 0, 0),
    adjust = list(decide = F, offset = 1000),
    family = "nbinom",
    nbinom = list(est = TRUE)
)
sim <- sim_data_mixture(settings = settings)
counts <- sim$other_data$counts
sim$other_data$design
ribo_example <- data.frame(geneid, counts)
colnames(ribo_example) <- 
    c('RiboCtrl1', 'RiboCtrl2', 'RiboTrt1', 'RiboTrt2',
      'RNACtrl1', 'RNACtrl2', 'RNATrt1', 'RNATrt2')
usethis::use_data(ribo_example, overwrite = TRUE)
