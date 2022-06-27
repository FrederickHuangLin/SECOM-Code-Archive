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

complex_data_generation = function(n, d, d1, abn_mean1, abn_prob1, 
                                   abn_mean2, abn_prob2, dispersion, seed) {
  set.seed(seed)
  
  #=============================Correlated pairs================================
  mu1 = sample(abn_mean1, d, replace = TRUE, prob = abn_prob1)
  mu2 = sample(abn_mean2, d, replace = TRUE, prob = abn_prob2)
  
  # Absolute abundances
  A1_1 = matrix(NA, ncol = n, nrow = d1)
  for (i in seq_len(d1)) {
    A1_1[i, ] = rnbinom(n = n, size = 1/dispersion, mu = mu1[i])
  }
  
  A2_1 = matrix(NA, ncol = n, nrow = d1)
  for (i in seq_len(d1)) {
    A2_1[i, ] = poly(x = A1_1[i, ], degree = 1, raw = FALSE)
    A2_1[i, ] = round(mu2[i] * A2_1[i, ]) + mu2[i]
  }
  
  R0 = matrix(0, nrow = d, ncol = d)
  R0_sub = diag(1, nrow = d1, ncol = d1)
  for (i in seq_len(d1)) {
    x1 = log(A1_1[i, ])
    x2 = log(A2_1[i, ])
    x1[is.infinite(x1)] = NA
    x2[is.infinite(x2)] = NA
    R0_sub[i, i] = cor(x1, x2, use = "pairwise.complete.obs")
  }
  R0[seq_len(d1), seq_len(d1)] = R0_sub
  
  #==================================Data 1=====================================
  A1 = matrix(NA, ncol = n, nrow = d)
  for (i in seq(d1 + 1, d)) {
    A1[i, ] = rnbinom(n = n, size = 1/dispersion, mu = mu1[i])
  }
  A1[seq_len(d1), ] = A1_1
  
  # Sequencing efficiency
  C1 = rbeta(n = d, shape1 = 5, shape2 = 5)
  
  # Microbial loads in the ecosystem
  A1_prim = A1 * C1
  A1_dot = colSums(A1_prim)
  
  # Relative abundances in the ecosystem
  R1 = A1_prim/t(replicate(d, A1_dot))
  
  # Sampling fractions
  S1 = rbeta(n = n, shape1 = 2, shape2 = 10)
  
  # Library sizes
  O1_dot = round(S1 * A1_dot)
  # Observed abundances
  O1 = matrix(NA, nrow = d, ncol = n)
  for (i in seq(n)) {
    O1[, i] = rmultinom(1, size = O1_dot[i], prob = R1[, i])
  }
  
  #==================================Data 2=====================================
  A2 = matrix(NA, ncol = n, nrow = d)
  for (i in seq(d1 + 1, d)) {
    A2[i, ] = rnbinom(n = n, size = 1/dispersion, mu = mu2[i])
  }
  A2[seq_len(d1), ] = A2_1
  
  # Sequencing efficiency
  C2 = rbeta(n = d, shape1 = 5, shape2 = 5)
  
  # Microbial loads in the ecosystem
  A2_prim = A2 * C2
  A2_dot = colSums(A2_prim)
  
  # Relative abundances in the ecosystem
  R2 = A2_prim/t(replicate(d, A2_dot))
  
  # Sampling fractions
  S2 = rbeta(n = n, shape1 = 2, shape2 = 10)
  
  # Library sizes
  O2_dot = round(S2 * A2_dot)
  # Observed abundances
  O2 = matrix(NA, nrow = d, ncol = n)
  for (i in seq(n)) {
    O2[, i] = rmultinom(1, size = O2_dot[i], prob = R2[, i])
  }
  
  res = list(O1 = O1, O2 = O2, R0 = R0)
  return(res)
}

n_d = c("50_100", "100_200")
d1 = 50
abn_mean1 = c(2000, 10000, 40000, 100000)
abn_prob1 = c(0.1, 0.4, 0.4, 0.1)
abn_mean2 = c(2000, 10000, 40000, 100000)
abn_prob2 = c(0.1, 0.4, 0.4, 0.1)
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
  
  sim_data = complex_data_generation(n, d, d1, abn_mean1, abn_prob1, 
                                     abn_mean2, abn_prob2, dispersion, seed)
  O1 = sim_data$O1
  O2 = sim_data$O2
  R0 = sim_data$R0
  taxa_id1 = paste0("A", seq_len(d))
  taxa_id2 = paste0("B", seq_len(d))
  sample_id = paste0("S", seq_len(n))
  dimnames(O1) = list(taxa_id1, sample_id)
  dimnames(O2) = list(taxa_id2, sample_id)
  dimnames(R0) = list(taxa_id1, taxa_id2)
  
  OTU1 = otu_table(O1, taxa_are_rows = TRUE)
  OTU2 = otu_table(O2, taxa_are_rows = TRUE)
  meta_data = data.frame(sample_id = sample_id)
  META = sample_data(meta_data)
  sample_names(META) = meta_data$sample_id
  otu_data1 = phyloseq(OTU1, META)
  otu_data2 = phyloseq(OTU2, META)
  
  pseqs = list(c(otu_data1, otu_data1),
               c(otu_data2, otu_data2))
  pseudo = 0; prv_cut = 0.5; corr_cut = 0.5; lib_cut = 1000
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
  taxa_id = c(paste0("data1 - ", taxa_id1), paste0("data2 - ", taxa_id2))
  pos_idx = match(taxa_keep, taxa_id)
  R_hat_secom1 = matrix(0, ncol = 2 * d, nrow = 2 * d)
  R_hat_secom1[pos_idx, pos_idx] = res_linear$corr_th
  R_hat_secom2 = matrix(0, ncol = 2 * d, nrow = 2 * d)
  R_hat_secom2[pos_idx, pos_idx] = res_linear$corr_fl
  R_hat_secom3 = matrix(0, ncol = 2 * d, nrow = 2 * d)
  R_hat_secom3[pos_idx, pos_idx] = res_dist$dcorr_fl
  
  
  R_hat_secom1 = R_hat_secom1[seq_len(d), seq(d + 1, 2 * d)]
  dimnames(R_hat_secom1) = list(taxa_id1, taxa_id2)
  R_hat_secom2 = R_hat_secom2[seq_len(d), seq(d + 1, 2 * d)]
  dimnames(R_hat_secom2) = list(taxa_id1, taxa_id2)
  R_hat_secom3 = R_hat_secom3[seq_len(d), seq(d + 1, 2 * d)]
  dimnames(R_hat_secom3) = list(taxa_id1, taxa_id2)
  
  # Relative error of Frobenius norm
  rel_F_secom1 = norm(R_hat_secom1 - R0, type = "F")/norm(R0, type = "F")
  rel_F_secom2 = norm(R_hat_secom2 - R0, type = "F")/norm(R0, type = "F")
  rel_F_secom3 = norm(R_hat_secom3 - R0, type = "F")/norm(R0, type = "F")
  
  # Relative error of Spectral norm
  rel_S_secom1 = norm(R_hat_secom1 - R0, type = "F")/norm(R0, type = "2")
  rel_S_secom2 = norm(R_hat_secom2 - R0, type = "F")/norm(R0, type = "2")
  rel_S_secom3 = norm(R_hat_secom3 - R0, type = "F")/norm(R0, type = "2")
  
  # TPR
  true_ind = (R0 != 0)
  secom_ind1 = (R_hat_secom1 != 0)
  secom_ind2 = (R_hat_secom2 != 0)
  secom_ind3 = (R_hat_secom3 != 0)
  tpr_secom1 = sum(secom_ind1 * true_ind)/sum(true_ind)
  tpr_secom2 = sum(secom_ind2 * true_ind)/sum(true_ind)
  tpr_secom3 = sum(secom_ind3 * true_ind)/sum(true_ind)
  
  # FPR
  fpr_secom1 = sum(secom_ind1 * (!true_ind))/sum(!true_ind)
  fpr_secom2 = sum(secom_ind2 * (!true_ind))/sum(!true_ind)
  fpr_secom3 = sum(secom_ind3 * (!true_ind))/sum(!true_ind)
  
  c(rel_F_secom1, rel_S_secom1, tpr_secom1, fpr_secom1, 
    rel_F_secom2, rel_S_secom2, tpr_secom2, fpr_secom2,
    rel_F_secom3, rel_S_secom3, tpr_secom3, fpr_secom3)
}

stopCluster(cl)

write_csv(data.frame(res_sim), "complex_secom.csv")