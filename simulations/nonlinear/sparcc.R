pkg_list = c("doParallel", "doRNG", "Hmisc", "readr", "compositions", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

library(doParallel)
library(doRNG)
library(Hmisc)
library(readr)
library(tidyverse)
library(SpiecEasi)

cor2cov = function(R, std) {
  Sigma = outer(std, std) * R
  return(Sigma)
}

hard_thresh = function(R, th){
  R_th = R
  R_th[abs(R) <= th] = 0
  return(R_th)
}

nonlinear_data_generation = function(n, d, d1, abn_mean, abn_prob, 
                                     dispersion, seed) {
  set.seed(seed)
  d2 = d - d1
  
  # ==========================Correlated pairs==================================
  # Log-Normal distribution to mimic NB distribution
  template_var = sapply(abn_mean, function(x) 
    log((x + dispersion * x^2)/x^2 + 1))
  template_mean = log(abn_mean) - template_var/2
  
  mu = sample(template_mean, d1, replace = TRUE, prob = abn_prob)
  std = rep(NA, d1)
  for (i in seq_along(template_mean)) {
    std[mu == template_mean[i]] = sqrt(template_var[i])
  }
  Sigma = cor2cov(R = diag(1, nrow = d1), std = std)
  
  # Absolute abundance in log scale
  a = t(MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma))
  for (i in seq(2, d1, 2)) {
    a[i, ] = poly(x = a[i - 1, ], degree = 2, raw = FALSE)[, 2]
    a[i, ] = mu[i] * a[i, ] + mu[i]
  }
  A1 = round(exp(a))
  
  # ========================Uncorrelated pairs==================================
  mu = sample(abn_mean, d2, replace = TRUE, prob = abn_prob)
  # Absolute abundances
  A2 = matrix(NA, ncol = n, nrow = d2)
  for (i in seq_len(d2)) {
    A2[i, ] = rnbinom(n = n, size = 1/dispersion, mu = mu[i])
  }
  
  A = rbind(A1, A2)
  
  # Sequencing efficiency
  C = exp(rnorm(d, mean = 0, sd = 1))
  
  # Microbial loads in the ecosystem
  A_prim = A * C
  A_dot = colSums(A_prim)
  
  # Relative abundances in the ecosystem
  R = A_prim/t(replicate(d, A_dot))
  
  # Sampling fractions
  S = rbeta(n = n, shape1 = 2, shape2 = 10)
  
  # Library sizes
  O_dot = round(S * A_dot)
  # Observed abundances
  O = matrix(NA, nrow = d, ncol = n)
  for (i in seq(n)) {
    O[, i] = rmultinom(1, size = O_dot[i], prob = R[, i])
  }
  
  R0 = matrix(0, ncol = d, nrow = d)
  lmat = replicate(d1/2, 
                   matrix(1, ncol = 2, nrow = 2), 
                   simplify = FALSE)
  R0_sub = as.matrix(Matrix::bdiag(lmat))
  R0[seq_len(d1), seq_len(d1)] = R0_sub
  
  res = list(O = O, R0 = R0)
  return(res)
}

n_d = c("50_100", "100_200")
d1 = 50
abn_mean = c(2000, 10000, 40000, 100000)
abn_prob = c(0.1, 0.4, 0.4, 0.1)
dispersion = c(0.5, 2)
iter_num = 100
seed = seq_len(iter_num)

simparams = data.frame(expand.grid(n_d, dispersion, seed)) %>%
  separate(col = Var1, into = c("n", "d"), sep = "_") %>%
  mutate(n = as.numeric(n),
         d = as.numeric(d))
colnames(simparams) = c("n", "d", "dispersion", "seed")
simparams = simparams %>%
  arrange(n, d, dispersion, seed)
simparams_list = apply(simparams, 1, paste0, collapse = "_")

cl = makeCluster(5)
registerDoParallel(cl)

res_sim = foreach(i = simparams_list, .combine = rbind, .verbose = TRUE, .packages = c("SpiecEasi")) %dorng% {
  params = strsplit(i, "_")[[1]]
  n = as.numeric(params[1])
  d = as.numeric(params[2])
  dispersion = as.numeric(params[3])
  seed = as.numeric(params[4])
  
  sim_data = nonlinear_data_generation(n, d, d1, abn_mean, abn_prob, dispersion, seed)
  O = sim_data$O
  R0 = sim_data$R0
  taxa_id = paste0("T", seq_len(d))
  sample_id = paste0("S", seq_len(n))
  dimnames(O) = list(taxa_id, sample_id)
  
  R_hat_sparcc = sparcc(t(O))$Cor
  R_hat_sparcc = hard_thresh(R_hat_sparcc, 0.3)
  dimnames(R_hat_sparcc) = list(taxa_id, taxa_id)
  
  # TPR
  true_ind = (R0[lower.tri(R0)] != 0)
  sparcc_ind = (R_hat_sparcc[lower.tri(R_hat_sparcc)] != 0)
  tpr_sparcc = sum(sparcc_ind * true_ind)/sum(true_ind)
  
  # FPR
  fpr_sparcc = sum(sparcc_ind * (!true_ind))/sum(!true_ind)
  
  c(tpr_sparcc, fpr_sparcc)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "nonlinear_sparcc.csv")