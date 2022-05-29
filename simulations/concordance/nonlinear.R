pkg_list = c("doParallel", "doRNG", "Hmisc", "DescTools", "readr", "energy", "tidyverse")
pkg_new = pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if(length(pkg_new)) install.packages(pkg_new)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!"microbiome" %in% installed.packages()[, "Package"])
  BiocManager::install("microbiome")

library(doParallel)
library(doRNG)
library(DescTools)
library(Hmisc)
library(energy)
library(readr)
library(tidyverse)
library(microbiome)

source("programs/00_secom.R")

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

n_d = c("50_100")
d1 = 50
abn_mean = c(2000, 10000, 40000, 100000)
abn_prob = c(0.1, 0.4, 0.4, 0.1)
dispersion = 0.5
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
  pseudo = 0; zero_cut = 0.5; corr_cut = 0.5; lib_cut = 1000
  wins_quant = c(0, 1); method = "pearson"; soft = FALSE; thresh_len = 20
  n_cv = 10; seed = 123; thresh_hard = 0.3; max_p = 0.001; n_cl = 1
  
  res_linear = secom_linear(pseqs, pseudo, zero_cut, corr_cut, lib_cut, 
                            wins_quant, method, soft, thresh_len, n_cv, 
                            seed, thresh_hard, max_p, n_cl)
  R = 1000; max_p = 0.001
  res_dist = secom_dist(pseqs, pseudo, zero_cut, corr_cut, lib_cut, 
                        wins_quant, R, seed, max_p, n_cl)
  
  taxa_keep = rownames(res_linear$corr)
  pos_idx = match(taxa_keep, taxa_id)
  
  R_hat_secom1 = matrix(0, ncol = d, nrow = d)
  R_hat_secom1[pos_idx, pos_idx] = res_linear$corr_fl
  R_hat_secom2 = matrix(0, ncol = d, nrow = d)
  R_hat_secom2[pos_idx, pos_idx] = res_dist$dcorr_fl
  
  R11 = (R_hat_secom1 != 0) * (R_hat_secom2 != 0)
  R10 = (R_hat_secom1 != 0) * (R_hat_secom2 == 0)
  R01 = (R_hat_secom1 == 0) * (R_hat_secom2 != 0)
  R00 = (R_hat_secom1 == 0) * (R_hat_secom2 == 0)
  
  R11 = R11[lower.tri(R11)]
  R10 = R10[lower.tri(R10)]
  R01 = R01[lower.tri(R01)]
  R00 = R00[lower.tri(R00)]
  
  tp = (R0[lower.tri(R0)] != 0)
  tn = (R0[lower.tri(R0)] == 0)
  
  # TP
  tp11 = sum(R11 * tp)
  tp10 = sum(R10 * tp)
  tp01 = sum(R01 * tp)
  
  # FP
  fp11 = sum(R11 * tn)
  fp10 = sum(R10 * tn)
  fp01 = sum(R01 * tn)
  
  # TN/FN
  tn00 = sum(R00 * tn)
  fn00 = sum(R00 * tp)
  
  # Compare two measures
  count11 = sum(R11)
  count10 = sum(R10)
  count01 = sum(R01)
  count00 = sum(R00)
  
  c(count11, count10, count01, count00, 
    tp11, tp10, tp01, fp11, fp10, fp01, tn00, fn00)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "nonlinear_conc.csv")