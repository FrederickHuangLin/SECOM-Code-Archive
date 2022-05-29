library(doParallel)
library(doRNG)
library(tidyverse)
library(mice)

cor2cov = function(R, std) {
  Sigma = outer(std, std) * R
  return(Sigma)
}

iter_num = 100
n = 100
abn_mean = c(200, 200)
dispersion = 2

cl = makeCluster(5)
registerDoParallel(cl)

res = foreach(i = seq_len(iter_num), .combine = rbind, .verbose = TRUE, .packages = c("mice")) %dorng% {
  set.seed(i)
  
  template_var = sapply(abn_mean, function(x) 
    log((x + dispersion * x^2)/x^2 + 1))
  template_mean = log(abn_mean) - template_var/2
  mu = template_mean
  std = rep(NA, length(mu))
  for (i in seq_along(template_mean)) {
    std[mu == template_mean[i]] = sqrt(template_var[i])
  }
  R = matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  Sigma = cor2cov(R = R, std = std)
  
  a = t(MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma))
  A = round(exp(a))
  
  S_list = c(0.1, 0.05, 0.02, 0.01, 0.005)
  res2 = data.frame(S = S_list, zero_prop = NA, 
                    r_complete = NA, lo_complete = NA, up_complete = NA,
                    r_impute = NA, lo_impute = NA, up_impute = NA)
  
  for (i in seq_along(S_list)) {
    S = S_list[i]
    O = floor(A * S)
    res2$zero_prop[i] = sum(O == 0) / (2 * n)
    
    o = log(O); o[is.infinite(o)] = NA
    df = data.frame(t(o))
    colnames(df) = c("o1", "o2")
    
    df_impute = mice(df, m = 1, maxit = 50, meth = "pmm", seed = 123, print = FALSE)
    df_fill = complete(df_impute, 1)
    
    res_complete = cor.test(df$o1, df$o2)
    res2$r_complete[i] = res_complete$estimate
    res2$lo_complete[i] = res_complete$conf.int[1]
    res2$up_complete[i] = res_complete$conf.int[2]
    
    res_impute = cor.test(df_fill$o1, df_fill$o2)
    res2$r_impute[i] = res_impute$estimate
    res2$lo_impute[i] = res_impute$conf.int[1]
    res2$up_impute[i] = res_impute$conf.int[2]
  }
  
  res2
}

stopCluster(cl)

write_csv(data.frame(res), "cca_vs_mice_corr.csv")



