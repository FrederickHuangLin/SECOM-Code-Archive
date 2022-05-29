pkg_list = c("doParallel", "doRNG", "Hmisc", "readr", "tidyverse", "Matrix")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

library(doParallel)
library(doRNG)
library(Hmisc)
library(readr)
library(tidyverse)
library(SpiecEasi)
library(Matrix)

cor2cov = function(R, std) {
  Sigma = outer(std, std) * R
  return(Sigma)
}

hard_thresh = function(R, th){
  R_th = R
  R_th[abs(R) <= th] = 0
  return(R_th)
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

cl = makeCluster(10)
registerDoParallel(cl)

res_sim = foreach(i = simparams_list, .combine = rbind, .verbose = TRUE, .packages = c("SpiecEasi", "Matrix")) %dorng% {
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
  
  se_mb = spiec.easi(t(O), method = "mb", lambda.min.ratio = 1e-2,
                     nlambda = 20, pulsar.params = list(rep.num = 50))
  se_gl = spiec.easi(t(O), method = "glasso", lambda.min.ratio = 1e-2,
                     nlambda = 20, pulsar.params = list(rep.num = 50))
  
  se_beta = symBeta(getOptBeta(se_mb), mode = "maxabs")
  se_cor  = cov2cor(getOptCov(se_gl))
  
  R_hat_se1 = as.matrix(se_beta)
  R_hat_se2 = as.matrix(se_cor * getRefit(se_gl))
  
  # Relative error of Frobenius norm
  rel_F_se1 = norm(R_hat_se1 - R0, type = "F")/norm(R0, type = "F")
  rel_F_se2 = norm(R_hat_se2 - R0, type = "F")/norm(R0, type = "F")
  
  # Relative error of Spectral norm
  rel_S_se1 = norm(R_hat_se1 - R0, type = "F")/norm(R0, type = "2")
  rel_S_se2 = norm(R_hat_se2 - R0, type = "F")/norm(R0, type = "2")
  
  # TPR
  true_ind = (R0[lower.tri(R0)] != 0)
  se_ind1 = (R_hat_se1[lower.tri(R_hat_se1)] != 0)
  se_ind2 = (R_hat_se2[lower.tri(R_hat_se2)] != 0)
  tpr_se1 = sum(se_ind1 * true_ind)/sum(true_ind)
  tpr_se2 = sum(se_ind2 * true_ind)/sum(true_ind)
  
  # FPR
  fpr_se1 = sum(se_ind1 * (!true_ind))/sum(!true_ind)
  fpr_se2 = sum(se_ind2 * (!true_ind))/sum(!true_ind)
  
  c(rel_F_se1, rel_S_se1, tpr_se1, fpr_se1, 
    rel_F_se2, rel_S_se2, tpr_se2, fpr_se2)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "linear_se.csv")