library(doParallel)
library(doRNG)
library(tidyverse)
library(mice)
library(gbm)

cor2cov = function(R, std) {
  Sigma = outer(std, std) * R
  return(Sigma)
}

iter_num = 100
n = 100
d = 50
abn_mean = c(200, 1000, 10000)
dispersion = 2

cl = makeCluster(5)
registerDoParallel(cl)

res = foreach(iter = seq_len(iter_num), .combine = rbind, .verbose = TRUE, .packages = c("mice", "gbm")) %dorng% {
  set.seed(iter)
  
  template_var = sapply(abn_mean, function(x) 
    log((x + dispersion * x^2)/x^2 + 1))
  template_mean = log(abn_mean) - template_var/2
  
  mu = c(template_mean[1], template_mean[1],
         sample(template_mean, d - 2, replace = TRUE))
  std = rep(NA, length(mu))
  for (i in seq_along(template_mean)) {
    std[mu == template_mean[i]] = sqrt(template_var[i])
  }
  R = diag(1, nrow = length(mu))
  Sigma = cor2cov(R = R, std = std)
  
  a = t(MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma))
  A = round(exp(a))
  
  S_list = c(0.1, 0.05, 0.02, 0.01, 0.005)
  res1 = data.frame(S = S_list, zero_prop = NA, 
                    r_cca = NA, r_mice = NA, r_gbm = NA)
  
  for (i in seq_along(S_list)) {
    S = S_list[i]
    O = floor(A * S)
    res1$zero_prop[i] = sum(O == 0) / (2 * n)
    
    o = log(O); o[is.infinite(o)] = NA
    df = data.frame(t(o))
    colnames(df) = paste0("T", seq_len(d))
    
    # CCA
    res_cca = cor.test(df$T1, df$T2)
    res1$r_cca[i] = res_cca$estimate
    
    # Mice
    df_mice = mice(df[, c("T1", "T2")], m = 1, maxit = 50, 
                   meth = "pmm", seed = 123, print = FALSE)
    df_fill = complete(df_mice, 1)
    res_mice = cor.test(df_fill$T1, df_fill$T2)
    res1$r_mice[i] = res_mice$estimate
    
    # GBM
    df_gbm = df
    T1_good = which(!is.na(df_gbm$T1))
    T1_bad = which(is.na(df_gbm$T1))
    df_train = df_gbm[T1_good, ]
    df_test = df_gbm[T1_bad, ]
    gbm_fit1 = gbm(T1 ~ ., data = df_train, distribution="gaussian",
                   n.trees = 100, shrinkage = 0.1, interaction.depth = 1, 
                   bag.fraction = 0.5, train.fraction = 1, n.minobsinnode = 1,
                   cv.folds = 0, keep.data = TRUE, verbose = FALSE)
    T1_pred = predict(gbm_fit1, newdata = df_test[, -1], n.trees = 100)
    df_gbm$T1[T1_bad] = T1_pred
    
    T2_good = which(!is.na(df_gbm$T2))
    T2_bad = which(is.na(df_gbm$T2))
    df_train = df_gbm[T2_good, ]
    df_test = df_gbm[T2_bad, ]
    gbm_fiT2 = gbm(T2 ~ ., data = df_train, distribution = "gaussian",
                   n.trees = 100, shrinkage = 0.1, interaction.depth = 1, 
                   bag.fraction = 0.5, train.fraction = 1, n.minobsinnode = 1,
                   cv.folds = 0, keep.data = TRUE, verbose = FALSE)
    T2_pred = predict(gbm_fiT2, newdata = df_test[, -2], n.trees = 100)
    df_gbm$T2[T2_bad] = T2_pred
    
    res_gbm = cor.test(df_gbm$T1, df_gbm$T2)
    res1$r_gbm[i] = res_gbm$estimate
  }
  
  res1
}

stopCluster(cl)

write_csv(data.frame(res), "cca_vs_impute_uncorr.csv")



