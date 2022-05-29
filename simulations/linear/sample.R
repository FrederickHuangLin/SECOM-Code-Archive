pkg_list = c("doParallel", "doRNG", "Hmisc", "readr", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

library(doParallel)
library(doRNG)
library(Hmisc)
library(readr)
library(tidyverse)

cor2cov = function(R, std) {
  Sigma = outer(std, std) * R
  return(Sigma)
}

p_filter = function(mat, mat_p, max_p){
  ind_p = mat_p
  ind_p[mat_p > max_p] = 0
  ind_p[mat_p <= max_p] = 1
  
  mat_filter = mat * ind_p
  return(mat_filter)
}

linear_data_generation = function(n, d, d1, corr_mu, corr_prob, 
                                  uncorr_mu, uncorr_prob, dispersion, seed) {
  set.seed(seed)
  d2 = d - d1
  
  mu = c(sample(corr_mu, d1, replace = TRUE, prob = corr_prob), # Correlated taxa
         sample(uncorr_mu, d2, replace = TRUE, prob = uncorr_prob)) # Uncorrelated taxa
  
  # Absolute abundances
  A = matrix(NA, ncol = n, nrow = d)
  for (i in seq_len(d)) {
    A[i, ] = rnbinom(n = n, size = 1/dispersion, mu = mu[i])
  }
  
  for (i in seq(2, d1, 2)) {
    A[i, ] = poly(x = A[i - 1, ], degree = 1, raw = FALSE)
    A[i, ] = round(mu[i] * A[i, ]) + mu[i]
  }
  
  R0 = matrix(0, ncol = d, nrow = d)
  lmat = replicate(d1/2, 
                   matrix(1, ncol = 2, nrow = 2), 
                   simplify = FALSE)
  R0_sub = as.matrix(Matrix::bdiag(lmat))
  R0_sub[R0_sub == 1] = cor(t(log(A[seq_len(d1), ] + 1)))[R0_sub == 1]
  R0[seq_len(d1), seq_len(d1)] = R0_sub
  
  # Sequencing efficiency
  C = rbeta(n = d, shape1 = 5, shape2 = 5)
  
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
  
  res = list(O = O, R0 = R0)
  return(res)
}

n_d = c("50_100", "100_200")
d1 = 50
corr_mu = c(2000, 10000, 40000, 100000)
corr_prob = c(0.1, 0.4, 0.4, 0.1)
uncorr_mu = c(2000, 10000, 40000, 100000)
uncorr_prob = c(0.1, 0.4, 0.4, 0.1)
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

res_sim = foreach(i = simparams_list, .verbose = TRUE, .combine = rbind) %dorng% {
  params = strsplit(i, "_")[[1]]
  n = as.numeric(params[1])
  d = as.numeric(params[2])
  dispersion = as.numeric(params[3])
  seed = as.numeric(params[4])
  
  sim_data = linear_data_generation(n, d, d1, corr_mu, corr_prob, 
                                    uncorr_mu, uncorr_prob, dispersion, seed)
  O = sim_data$O
  R0 = sim_data$R0
  taxa_id = paste0("T", seq_len(d))
  sample_id = paste0("S", seq_len(n))
  dimnames(O) = list(taxa_id, sample_id)
  o = log(O)
  o[is.infinite(o)] = NA
  
  res_sample = Hmisc::rcorr(t(O), type = "pearson")
  corr_sample = res_sample$r
  corr_sample_p = res_sample$P
  diag(corr_sample_p) = 0
  R_hat_sample = p_filter(corr_sample, corr_sample_p, max_p = 0.001)
  
  # Relative error of Frobenius norm
  rel_F_sample = norm(R_hat_sample - R0, type = "F")/norm(R0, type = "F")
  
  # Relative error of Spectral norm
  rel_S_sample = norm(R_hat_sample - R0, type = "F")/norm(R0, type = "2")
  
  # TPR
  true_ind = (R0[lower.tri(R0)] != 0)
  sample_ind = (R_hat_sample[lower.tri(R_hat_sample)] != 0)
  tpr_sample = sum(sample_ind * true_ind)/sum(true_ind)
  
  # FPR
  fpr_sample = sum(sample_ind * (!true_ind))/sum(!true_ind)
  
  c(rel_F_sample, rel_S_sample, tpr_sample, fpr_sample)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "linear_sample.csv")