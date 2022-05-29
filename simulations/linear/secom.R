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

res_sim = foreach(i = simparams_list, .combine = rbind, .verbose = TRUE, .packages = c("microbiome", "tidyverse", "doParallel")) %dorng% {
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
  R_hat_secom1[pos_idx, pos_idx] = res_linear$corr_th
  R_hat_secom2 = matrix(0, ncol = d, nrow = d)
  R_hat_secom2[pos_idx, pos_idx] = res_linear$corr_fl
  R_hat_secom3 = matrix(0, ncol = d, nrow = d)
  R_hat_secom3[pos_idx, pos_idx] = res_dist$dcorr_fl
  
  # Relative error of Frobenius norm
  rel_F_secom1 = norm(R_hat_secom1 - R0, type = "F")/norm(R0, type = "F")
  rel_F_secom2 = norm(R_hat_secom2 - R0, type = "F")/norm(R0, type = "F")
  rel_F_secom3 = norm(R_hat_secom3 - R0, type = "F")/norm(R0, type = "F")
  
  # Relative error of Spectral norm
  rel_S_secom1 = norm(R_hat_secom1 - R0, type = "F")/norm(R0, type = "2")
  rel_S_secom2 = norm(R_hat_secom2 - R0, type = "F")/norm(R0, type = "2")
  rel_S_secom3 = norm(R_hat_secom3 - R0, type = "F")/norm(R0, type = "2")
  
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
  
  # Compare two measures
  count11 = sum(secom_ind2 == 1 & secom_ind3 == 1)
  count10 = sum(secom_ind2 == 1 & secom_ind3 == 0)
  count01 = sum(secom_ind2 == 0 & secom_ind3 == 1)
  count00 = sum(secom_ind2 == 0 & secom_ind3 == 0)
  
  c(rel_F_secom1, rel_S_secom1, tpr_secom1, fpr_secom1, 
    rel_F_secom2, rel_S_secom2, tpr_secom2, fpr_secom2,
    rel_F_secom3, rel_S_secom3, tpr_secom3, fpr_secom3,
    count11, count10, count01, count00)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "linear_secom.csv")