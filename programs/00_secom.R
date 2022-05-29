library(doParallel)
library(doRNG)
library(DescTools)
library(Hmisc)
library(energy)
library(microbiome)

# Sampling fraction difference estimation
s_diff_est = function(pseq, zero_cut, lib_cut) {
  # Data pre-processing
  O = abundances(pseq)
  n = dim(O)[2]
  d = dim(O)[1]
  
  # Filter taxa/samples with excess zeros
  tax_keep = which(rowSums(O == 0, na.rm = TRUE)/n <= zero_cut)
  samp_keep = which(colSums(O, na.rm = TRUE) > lib_cut)
  O_cut = O[tax_keep, samp_keep]
  
  o = log(O_cut)

  o[is.infinite(o)] = NA
  o_center = o - rowMeans(o, na.rm = TRUE)

  # Estimate weights
  wt = apply(o_center, 1, function(x) 1/var(x, na.rm = TRUE))
  o_center = o_center[is.finite(wt), ]
  wt = wt[is.finite(wt)]
  
  # Estimate sampling fraction difference
  s_diff_hat = apply(o_center, 2, function(x) {
    weighted.mean(x, wt, na.rm = TRUE)}
  )

  return(s_diff_hat)
}

# Absolute abundance estimation
abs_est = function(pseqs, pseudo, zero_cut, lib_cut) {
  pseq1 = pseqs[[1]]
  pseq2 = pseqs[[2]]
  
  if (! all(sample_names(pseq1) == sample_names(pseq2))) {
    stop("Sample names of two phyloseq objects does not match.")
  }
  
  # Sampling fraction difference estimation
  s_diff_hat = s_diff_est(pseq1, zero_cut, lib_cut)
  
  # Data pre-processing
  O = abundances(pseq2)
  n = dim(O)[2]
  d = dim(O)[1]
  tax_keep = which(rowSums(O == 0, na.rm = TRUE)/n <= zero_cut)
  samp_keep = names(s_diff_hat)
  O_cut = O[tax_keep, samp_keep] + pseudo
  
  o = log(O_cut)
  o[is.infinite(o)] = NA
  n = ncol(o)
  d = nrow(o)

  # Absolute abundance estimation
  y_hat = o - rowMeans(o, na.rm = TRUE) - t(replicate(d, s_diff_hat))
  
  abs_list = list(s_diff_hat = s_diff_hat, y_hat = y_hat)
  return(abs_list)
}

# Sparse estimation on linear correlations
sparse_linear = function(mat, wins_quant, method, soft, thresh_len, 
                         n_cv, seed, thresh_hard, max_p) {
  # Thresholding
  mat_thresh = function(mat, th, soft){
    mat_sign = sign(mat)
    mat_th = mat
    mat_th[abs(mat) <= th] = 0
    if (soft) {
      mat_th[abs(mat) > th] = abs(mat_th[abs(mat) > th]) - th
      mat_th = mat_th * mat_sign
    }
    return(mat_th)
  }

  # Threshold loss function
  thresh_loss = function(mat1, mat2, method, th, soft) {
    corr1 = cor(mat1, method = method, use = "pairwise.complete.obs")
    corr2 = cor(mat2, method = method, use = "pairwise.complete.obs")
    corr_diff = mat_thresh(corr1, th, soft) - corr2
    corr_diff[is.na(corr_diff)] = 0
    loss = norm(corr_diff, type = "F")
    return(loss)
  }

  # Filtering based on p-values
  p_filter = function(mat, mat_p, max_p){
    ind_p = mat_p
    ind_p[mat_p > max_p] = 0
    ind_p[mat_p <= max_p] = 1

    mat_filter = mat * ind_p
    return(mat_filter)
  }
  
  # Sort taxa
  sort_taxa = sort(colnames(mat))
  mat = mat[, sort_taxa]

  # Winsorization
  mat = apply(mat, 2, function(x) 
    DescTools::Winsorize(x, probs = wins_quant, na.rm = TRUE))
  
  # Co-occurrence matrix
  mat_occur = mat
  mat_occur[mat_occur != 0] = 1
  mat_occur[mat_occur == 0] = 0
  mat_occur[is.na(mat_occur)] = 0
  
  df_occur = as.data.frame(mat_occur) %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(cols = -sample_id, names_to = "taxon", values_to = "occur") %>%
    filter(occur == 1)
  
  mat_cooccur = crossprod(table(df_occur[, 1:2]))
  diag(mat_cooccur) = colSums(mat_occur)
  
  if (any(mat_cooccur < 10)) {
    warn_txt = sprintf(paste0("There are some pairs of taxa that have insufficient (< 10) overlapping samples.\n",
                              "Proceed with caution since the point estimates for these pairs are unstable.\n", 
                              "For pairs of taxa with no overlapping samples, the point estimates will be replaced with 0s,\n",
                              "and the corresponding p-values will be replaced with 1s.\n",
                              "Please check `mat_cooccur` for details about the co-occurrence pattern."))
    warning(warn_txt)
  }

  # Sample size for training and test sets
  n = dim(mat)[1]
  n1 = n - floor(n/log(n))
  n2 = n - n1
  d = dim(mat)[2]

  # Correlation matrix
  corr_list = suppressWarnings(Hmisc::rcorr(x = mat, type = method))
  corr = corr_list$r
  corr[mat_cooccur == 0] = 0

  # Cross-Validation
  max_thresh = max(abs(corr[corr != 1]), na.rm = TRUE)
  thresh_grid = seq(from = 0, to = max_thresh, length.out = thresh_len)
  
  set.seed(seed)
  loss_mat = foreach(i = seq_len(n_cv), .combine = rbind) %dorng% {
    index = sample(seq_len(n), size = n1, replace = FALSE)
    mat1 = mat[index,]
    mat2 = mat[-index,]
    loss = sapply(thresh_grid, FUN = thresh_loss,
                  mat1 = mat1, mat2 = mat2, 
                  method = method, soft = soft)
  }

  # Correlation matrix after thresholding
  loss_vec = colMeans(loss_mat)
  thresh_opt = thresh_grid[which.min(loss_vec)]
  if (thresh_opt >= thresh_hard) {
    corr_th = mat_thresh(mat = corr, th = thresh_opt, soft = soft)
  } else {
    corr_th = mat_thresh(mat = corr, th = thresh_hard, soft = soft)
  }

  # Correlation matrix after filtering
  corr_p = corr_list$P
  diag(corr_p) = 0
  corr_p[mat_cooccur == 0] = 1
  corr_p[is.na(corr_p)] = 1
  corr_fl = p_filter(mat = corr, mat_p = corr_p, max_p = max_p)

  # Output
  result = list(cv_error = loss_vec,
                n_cv = n_cv,
                thresh_grid = thresh_grid,
                thresh_opt = thresh_opt,
                mat_cooccur = mat_cooccur,
                corr = corr,
                corr_p = corr_p,
                corr_th = corr_th,
                corr_fl = corr_fl)
  return(result)
}

# Sparse estimation on distance correlations
sparse_dist = function(mat, wins_quant, R, seed, max_p) {
  # Filtering based on p-values
  p_filter = function(mat, mat_p, max_p){
    ind_p = mat_p
    ind_p[mat_p > max_p] = 0
    ind_p[mat_p <= max_p] = 1
    
    mat_filter = mat * ind_p
    return(mat_filter)
  }
  
  # Sort taxa
  sort_taxa = sort(colnames(mat))
  mat = mat[, sort_taxa]
  
  # Winsorization
  mat = apply(mat, 2, function(x) 
    DescTools::Winsorize(x, probs = wins_quant, na.rm = TRUE))
  
  # Co-occurrence matrix
  mat_occur = mat
  mat_occur[mat_occur != 0] = 1
  mat_occur[mat_occur == 0] = 0
  mat_occur[is.na(mat_occur)] = 0
  
  df_occur = as.data.frame(mat_occur) %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(cols = -sample_id, names_to = "taxon", values_to = "occur") %>%
    filter(occur == 1)
  
  mat_cooccur = crossprod(table(df_occur[, 1:2]))
  diag(mat_cooccur) = colSums(mat_occur)
  
  if (any(mat_cooccur < 10)) {
    warn_txt = sprintf(paste0("There are some pairs of taxa that have insufficient (< 10) overlapping samples.\n",
                              "Proceed with caution since the point estimates for these pairs are unstable.\n", 
                              "For pairs of taxa with no overlapping samples, the point estimates will be replaced with 0s,\n",
                              "and the corresponding p-values will be replaced with 1s.\n",
                              "Please check `mat_cooccur` for details about the co-occurrence pattern."))
    warning(warn_txt)
  }
  
  # Calculation
  d = dim(mat)[2]
  taxanames = colnames(mat)
  comb = function(...) {
    mapply('rbind', ..., SIMPLIFY = FALSE)
  }
  
  set.seed(seed)
  if (is.null(R)) {
    dcorr_list = foreach(i = seq_len(d - 1), .combine = 'comb', 
                         .multicombine = TRUE, .packages = "energy") %dorng% {
       dcorr_i = rep(NA, d)
       dcorr_p_i = rep(NA, d)
       
       mat_i = mat[!is.na(mat[, i]), ]
       x = mat_i[, i]
       
       # Distance correlation
       dcorr_i[(i + 1):d] = apply(mat_i[, (i + 1):d, drop = FALSE], 2, 
                                  function(y) {
                                    z = x[!is.na(y)]
                                    y = y[!is.na(y)]
                                    dcor(z, y, index = 1.0)
                                  })
       
       # P-values
       dcorr_p_i[(i + 1):d] = apply(mat_i[, (i + 1):d, drop = FALSE], 2, 
                                    function(y) {
                                      z = x[!is.na(y)]
                                      y = y[!is.na(y)]
                                      dcorT.test(z, y)$p.value
                                    })
       
       list(dcorr_i, dcorr_p_i)
     }
  } else {
    dcorr_list = foreach(i = seq_len(d - 1), .combine = 'comb', 
                         .multicombine = TRUE, .packages = "energy") %dorng% {
       dcorr_i = rep(NA, d)
       dcorr_p_i = rep(NA, d)
       
       mat_i = mat[!is.na(mat[, i]), ]
       x = mat_i[, i]
       
       # Distance correlation
       dcorr_i[(i + 1):d] = apply(mat_i[, (i + 1):d, drop = FALSE], 2, 
                                  function(y) {
                                    z = x[!is.na(y)]
                                    y = y[!is.na(y)]
                                    dcor(z, y, index = 1.0)
                                  })
       
       # P-values
       dcorr_p_i[(i + 1):d] = apply(mat_i[, (i + 1):d, drop = FALSE], 2, 
                                    function(y) {
                                      z = x[!is.na(y)]
                                      y = y[!is.na(y)]
                                      dcor.test(z, y, index = 1.0, R = R)$p.value
                                    })
       
       list(dcorr_i, dcorr_p_i)
     }
  }
  
  dcorr = rbind(dcorr_list[[1]], rep(NA, d))
  dcorr_p = rbind(dcorr_list[[2]], rep(NA, d))
  # Symmetrize the matrix
  dcorr[lower.tri(dcorr)] = t(dcorr)[lower.tri(dcorr)]
  diag(dcorr) = 1
  dcorr[mat_cooccur == 0] = 0
  dcorr_p[lower.tri(dcorr_p)] = t(dcorr_p)[lower.tri(dcorr_p)]
  diag(dcorr_p) = 0
  dcorr_p[mat_cooccur == 0] = 1
  dcorr_p[is.na(dcorr_p)] = 1
  dimnames(dcorr) = list(taxanames, taxanames)
  dimnames(dcorr_p) = list(taxanames, taxanames)
  
  dcorr_fl = p_filter(mat = dcorr, mat_p = dcorr_p, max_p = max_p)
  
  # Output
  result = list(mat_cooccur = mat_cooccur,
                dcorr = dcorr,
                dcorr_p = dcorr_p,
                dcorr_fl = dcorr_fl)
  return(result)
}

# Wrapper function
secom_linear = function(pseqs, pseudo = 0, zero_cut = 0.5, corr_cut = 0.5, 
                        lib_cut = 1000, wins_quant = c(0.05, 0.95), 
                        method = c("pearson", "kendall", "spearman"),
                        soft = FALSE, thresh_len = 100, n_cv = 10, seed = 123, 
                        thresh_hard = 0.3, max_p = 0.005, n_cl = 1) {
  # ===========Sampling fraction and absolute abundance estimation==============
  if (length(pseqs) == 1) {
    abs_list = abs_est(pseqs[[1]], pseudo, zero_cut, lib_cut)
    s_diff_hat = abs_list$s_diff_hat
    y_hat = abs_list$y_hat
  } else {
    if (is.null(names(pseqs))) names(pseqs) = paste0("data", seq_along(pseqs))
    
    # Check common samples
    samp_names = lapply(pseqs, function(x) sample_names(x[[2]]))
    samp_common = Reduce(intersect, samp_names)
    samp_txt = sprintf(paste0("The number of samples that are common across datasets: ",
                              length(samp_common)))
    message(samp_txt)
    if (length(samp_common) < 10) {
      stop("The number of common samples is too small. Across-dataset computation is not recommended.")
    }
    
    # Rename taxa
    for (i in seq_along(pseqs)) {
      taxa_names(pseqs[[i]][[1]]) = paste(names(pseqs)[i], 
                                          taxa_names(pseqs[[i]][[1]]), 
                                          sep = " - ")
      taxa_names(pseqs[[i]][[2]]) = paste(names(pseqs)[i], 
                                          taxa_names(pseqs[[i]][[2]]), 
                                          sep = " - ")
    }
    abs_list = lapply(pseqs, function(x) abs_est(x, pseudo, zero_cut, lib_cut))
    s_diff_hat = lapply(abs_list, function(x) x$s_diff_hat)
    y_hat = bind_rows(lapply(abs_list, function(x) as.data.frame(x$y_hat)))
    y_hat = as.matrix(y_hat)
  }
  
  # =================Sparse estimation on linear correlations===================
  cl = makeCluster(n_cl)
  registerDoParallel(cl)
  
  if (method %in% c("pearson", "kendall", "spearman")) {
    res_corr = sparse_linear(mat = t(y_hat), wins_quant, method, soft, 
                             thresh_len, n_cv, seed, thresh_hard, max_p)
  } else {
    stop("The correlation coefficien should be one of 'pearson', 'kendall', 'spearman'.", call. = FALSE)
  }
  
  stopCluster(cl)
  
  # To prevent FP from taxa with extremely small variances
  if (length(pseqs) == 1) {
    corr_s = cor(cbind(s_diff_hat, t(y_hat)), use = "pairwise.complete.obs")[1, -1]
    fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
    fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
    fp_ind = (fp_ind1 * fp_ind2 == 1)
    diag(fp_ind) = FALSE
    res_corr$corr[fp_ind] = 0
    res_corr$corr_th[fp_ind] = 0
    res_corr$corr_fl[fp_ind] = 0
    res_corr$corr_p[fp_ind] = 1
  } else {
    for (i in seq_along(pseqs)) {
      df_s = data.frame(s = s_diff_hat[[i]]) %>%
        rownames_to_column("sample_id")
      df_y = as.data.frame(t(y_hat)) %>%
        rownames_to_column("sample_id")
      df_merge = df_s %>%
        right_join(df_y, by = "sample_id") %>%
        dplyr::select(-sample_id)
      corr_s = cor(df_merge, use = "pairwise.complete.obs")[1, -1]
      fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
      fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
      fp_ind = (fp_ind1 * fp_ind2 == 1)
      diag(fp_ind) = FALSE
      res_corr$corr[fp_ind] = 0
      res_corr$corr_th[fp_ind] = 0
      res_corr$corr_fl[fp_ind] = 0
      res_corr$corr_p[fp_ind] = 1
    }
  }
  
  # ==================================Outputs===================================
  res = c(list(s_diff_hat = s_diff_hat, y_hat = y_hat), res_corr)
  return(res)
}

secom_dist = function(pseqs, pseudo = 0, zero_cut = 0.5, corr_cut = 0.5, 
                      lib_cut = 1000, wins_quant = c(0.05, 0.95), R = 1000, 
                      seed = 123, max_p = 0.005, n_cl = 1) {
  # ===========Sampling fraction and absolute abundance estimation==============
  if (length(pseqs) == 1) {
    abs_list = abs_est(pseqs[[1]], pseudo, zero_cut, lib_cut)
    s_diff_hat = abs_list$s_diff_hat
    y_hat = abs_list$y_hat
  } else {
    if (is.null(names(pseqs))) names(pseqs) = paste0("data", seq_along(pseqs))
    
    # Check common samples
    samp_names = lapply(pseqs, function(x) sample_names(x[[2]]))
    samp_common = Reduce(intersect, samp_names)
    samp_txt = sprintf(paste0("The number of samples that are common across datasets: ",
                              length(samp_common)))
    message(samp_txt)
    if (length(samp_common) < 10) {
      stop("The number of common samples is too small. Across-dataset computation is not recommended.")
    }
    
    # Rename taxa
    for (i in seq_along(pseqs)) {
      taxa_names(pseqs[[i]][[1]]) = paste(names(pseqs)[i], 
                                          taxa_names(pseqs[[i]][[1]]), 
                                          sep = " - ")
      taxa_names(pseqs[[i]][[2]]) = paste(names(pseqs)[i], 
                                          taxa_names(pseqs[[i]][[2]]), 
                                          sep = " - ")
    }
    abs_list = lapply(pseqs, function(x) abs_est(x, pseudo, zero_cut, lib_cut))
    s_diff_hat = lapply(abs_list, function(x) x$s_diff_hat)
    y_hat = bind_rows(lapply(abs_list, function(x) as.data.frame(x$y_hat)))
    y_hat = as.matrix(y_hat)
  }
  
  # ================Sparse estimation on distance correlations==================
  cl = makeCluster(n_cl)
  registerDoParallel(cl)
  
  res_corr = sparse_dist(mat = t(y_hat), wins_quant, R, seed, max_p)
  
  stopCluster(cl)
  
  # To prevent FP from taxa with extremely small variances
  if (length(pseqs) == 1) {
    corr_s = cor(cbind(s_diff_hat, t(y_hat)), use = "pairwise.complete.obs")[1, -1]
    fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
    fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
    fp_ind = (fp_ind1 * fp_ind2 == 1)
    diag(fp_ind) = FALSE
    res_corr$dcorr[fp_ind] = 0
    res_corr$dcorr_fl[fp_ind] = 0
    res_corr$dcorr_p[fp_ind] = 1
  } else {
    for (i in seq_along(pseqs)) {
      df_s = data.frame(s = s_diff_hat[[i]]) %>%
        rownames_to_column("sample_id")
      df_y = as.data.frame(t(y_hat)) %>%
        rownames_to_column("sample_id")
      df_merge = df_s %>%
        right_join(df_y, by = "sample_id") %>%
        dplyr::select(-sample_id)
      corr_s = cor(df_merge, use = "pairwise.complete.obs")[1, -1]
      fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
      fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
      fp_ind = (fp_ind1 * fp_ind2 == 1)
      diag(fp_ind) = FALSE
      res_corr$dcorr[fp_ind] = 0
      res_corr$dcorr_fl[fp_ind] = 0
      res_corr$dcorr_p[fp_ind] = 1
    }
  }
  
  # ==================================Outputs===================================
  res = c(list(s_diff_hat = s_diff_hat, y_hat = y_hat), res_corr)
  return(res)
}
