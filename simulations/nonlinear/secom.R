pkg_list = c("doParallel", "doRNG", "Hmisc", "DescTools", "readr", "energy", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkg_list = c("microbiome", "ANCOMBC")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) BiocManager::install(pkg_new)

library(doParallel)
library(doRNG)
library(DescTools)
library(Hmisc)
library(energy)
library(readr)
library(tidyverse)
library(microbiome)
library(ANCOMBC)

cor2cov = function(R, std) {
  Sigma = outer(std, std) * R
  return(Sigma)
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
  C = C = rbeta(n = d, shape1 = 5, shape2 = 5)
  
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

cl = makeCluster(10)
registerDoParallel(cl)

res_sim = foreach(i = simparams_list, .combine = rbind, .verbose = TRUE, .packages = c("microbiome", "tidyverse", "doParallel")) %dorng% {
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
  meta_data = data.frame(sample_id = sample_id)
  dimnames(O) = list(taxa_id, sample_id)
  OTU = otu_table(O, taxa_are_rows = TRUE)
  META = sample_data(meta_data)
  sample_names(META) = meta_data$sample_id
  otu_data = phyloseq(OTU, META)
  
  pseqs = list(c(otu_data, otu_data))
  pseudo = 0; prv_cut = 0.5; lib_cut = 1000; corr_cut = 0.5
  wins_quant = c(0, 1); method = "pearson"; soft = FALSE; thresh_len = 20
  n_cv = 10; thresh_hard = 0; max_p = 0.001; n_cl = 1
  
  set.seed(123)
  res_linear = secom_linear(pseqs, pseudo, prv_cut, lib_cut, corr_cut, 
                            wins_quant, method, soft, thresh_len, n_cv, 
                            thresh_hard, max_p, n_cl)
  R = 1000; max_p = 0.001
  set.seed(123)
  res_dist = secom_dist(pseqs, pseudo, prv_cut, lib_cut, corr_cut, 
                        wins_quant, R, thresh_hard, max_p, n_cl)
  
  taxa_keep = rownames(res_linear$corr)
  pos_idx = match(taxa_keep, taxa_id)
  R_hat_secom1 = matrix(0, ncol = d, nrow = d)
  R_hat_secom1[pos_idx, pos_idx] = res_linear$corr_th
  R_hat_secom2 = matrix(0, ncol = d, nrow = d)
  R_hat_secom2[pos_idx, pos_idx] = res_linear$corr_fl
  R_hat_secom3 = matrix(0, ncol = d, nrow = d)
  R_hat_secom3[pos_idx, pos_idx] = res_dist$dcorr_fl
  
  # TPR
  true_ind = (R0[lower.tri(R0)] != 0)
  secom_ind1 = (R_hat_secom1[lower.tri(R_hat_secom1)] != 0)
  secom_ind2 = (R_hat_secom2[lower.tri(R_hat_secom2)] != 0)
  secom_ind3 = (R_hat_secom3[lower.tri(R_hat_secom3)] != 0)
  tpr_secom1 = sum(secom_ind1 * true_ind)/sum(true_ind)
  tpr_secom2 = sum(secom_ind2 * true_ind)/sum(true_ind)
  tpr_secom3 = sum(secom_ind3 * true_ind)/sum(true_ind)
  
  # FPR
  fpr_secom1 = sum(secom_ind1 * (!true_ind))/sum(!true_ind)
  fpr_secom2 = sum(secom_ind2 * (!true_ind))/sum(!true_ind)
  fpr_secom3 = sum(secom_ind3 * (!true_ind))/sum(!true_ind) 
  
  c(tpr_secom1, fpr_secom1, 
    tpr_secom2, fpr_secom2,
    tpr_secom3, fpr_secom3)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "nonlinear_secom.csv")