# Function for assessing the convergence of single-level FAST
RHat_FAST <- function(model, data, align_EF){
  n_chains = length(model@inits)  
  n_samp = model@stan_args[[1]]$iter - model@stan_args[[1]]$warmup
  
  # Extract by chain
  by_chain = FAST_byChain(model, data$N, data$M, data$K)
  
  aligned_objs = map(1:n_chains, function(c){
    aligned = align_weights(by_chain$Psi[[c]], 
                            by_chain$Score[[c]],
                            data$B_RE, 
                            anchor = align_EF)
    return(aligned)
  })
  
  score_chains = map(aligned_objs, function(x){
    out = x$Score %>% abind(along = 3)
    return(out)
  }) %>% abind(along = 4)
  
  lambda_chains = apply(score_chains, c(2,3,4), var)
  
  fpc_chains = map(aligned_objs, function(x){
    out = x$EF %>% abind(along = 3)
    return(out)
  }) %>% abind(along = 4)
  
  raw_samples = extract(model, permuted = F)
  dimen_names = dimnames(raw_samples)$parameters
  
  fpcs_idx = grepl("^H\\[[0-9]+\\]$", dimen_names)
  mus_idx = grepl("h_mu", dimen_names)
  smooth_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      fpc_smooth = raw_samples[i,j,fpcs_idx]
      mean_smooth = raw_samples[i,j,mus_idx]
      return(c(fpc_smooth, mean_smooth))
    }) %>% abind(along = 2)
  }) %>% abind(along = 3)
  
  sig_idx = grepl("sigma2", dimen_names)
  sigma_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      return(raw_samples[i,j,sig_idx])
    }) %>% abind(along = 1)
  }) %>% abind(along = 2)
  
  mu_idx = grepl("w_mu", dimen_names)
  mu_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      mean_func = (data$B_FE %*% raw_samples[i,j,mu_idx]) * data$sd_Y + data$mu_Y
      return(mean_func)
    }) %>% abind(along = 2)
  }) %>% abind(along = 3)
  
  # Calculate RHat statistics
  lambda_rhats = rep(0, length = data$K)
  score_rhats = matrix(0, nrow = data$N, ncol = data$K)
  fpc_rhats = matrix(0, nrow = data$M, ncol = data$K)
  sigma_rhat = Rhat(sigma_chains)
  mu_rhats = rep(0, length = data$M)
  smooth_rhats = rep(0, length = data$K + 1)

  for(m in 1:data$M){
    mu_rhats = Rhat(mu_chains[m,,])
  }
  for(k in 1:data$K){
    lambda_rhats[k] = Rhat(lambda_chains[k,,])
    smooth_rhats[k] = Rhat(smooth_chains[k,,])

    for(n in 1:data$N){
      score_rhats[n, k] = Rhat(score_chains[n,k,,])
    }
    
    for(m in 1:data$M){
      fpc_rhats[m, k] = Rhat(fpc_chains[m,k,,])
    }
  }
  smooth_rhats[data$K+1] = Rhat(smooth_chains[data$K,,])
  
  # Summarize and convert to dataframes
  fpc_names = paste0("FPC ", 1:data$K)
  
  score_rhats = data.frame(score_rhats)
  colnames(score_rhats) = fpc_names
  score_rhats = score_rhats %>%
    pivot_longer(cols = everything(), names_to = "FPC_Num", values_to = "Rhat") %>%
    group_by(FPC_Num) %>%
    summarize(Med_RHat = median(Rhat),
              Max_RHat = max(Rhat))
  
  func_rhats = data.frame(fpc_rhats)
  colnames(func_rhats) = fpc_names
  func_rhats = func_rhats %>%
    pivot_longer(cols = everything(), names_to = "Function", values_to = "Rhat") %>%
    group_by(Function) %>%
    summarize(Med_RHat = median(Rhat),
              Max_RHat = max(Rhat)) %>%
    rbind(data.frame(Function = "Mu", 
                     Med_RHat = median(mu_rhats),
                     Max_RHat = max(mu_rhats)))

  smooth_rhats = data.frame(Function = c(fpc_names, "Mu"), 
                            RHat = smooth_rhats)
  
  Variance_rhats = data.frame(Element = c(paste0("Lambda_", 1:data$K), "Sigma2"), 
                              RHat = c(lambda_rhats, sigma_rhat))
  
  outputs = list(Score = score_rhats, Func = func_rhats, 
                 Smoothing_Params = smooth_rhats, Variances = Variance_rhats)
  
  return(outputs)
}

# Function for assessing the convergence of single-level GFSR
RHat_GFSR <- function(model, data, align_EF){
  n_chains = length(model@inits)  
  n_samp = model@stan_args[[1]]$iter - model@stan_args[[1]]$warmup
  
  # Extract by chain
  by_chain = GFSR_byChain(model, data$I, data$D, data$Kp, data$BS)
  
  aligned_objs = map(1:n_chains, function(c){
    aligned = align_FPCs(by_chain$Phi[[c]], 
                         by_chain$Score[[c]],
                         anchor = align_EF)
    return(aligned)
  })
  
  score_chains = map(aligned_objs, function(x){
    out = x$Score %>% abind(along = 3)
    return(out)
  }) %>% abind(along = 4)
  
  lambda_chains = apply(score_chains, c(2,3,4), var)
  
  fpc_chains = map(aligned_objs, function(x){
    out = x$EF %>% abind(along = 3)
    return(out)
  }) %>% abind(along = 4)
  
  raw_samples = extract(model, permuted = F)
  dimen_names = dimnames(raw_samples)$parameters
  
  fpcs_idx = grepl("psi_sig", dimen_names)
  mus_idx = grepl("beta_sig", dimen_names)
  smooth_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      fpc_smooth = raw_samples[i,j,fpcs_idx]
      mean_smooth = raw_samples[i,j,mus_idx]
      return(c(fpc_smooth, mean_smooth))
    }) %>% abind(along = 2)
  }) %>% abind(along = 3)
  
  sig_idx = grepl("sigma2", dimen_names)
  sigma_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      return(raw_samples[i,j,sig_idx])
    }) %>% abind(along = 1)
  }) %>% abind(along = 2)
  
  mu_idx = grepl("^beta\\[[0-9]+,[0-9]+\\]$", dimen_names)
  mu_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      mean_func = data$BS %*% raw_samples[i,j,mu_idx]
      return(mean_func)
    }) %>% abind(along = 2)
  }) %>% abind(along = 3)
  
  # Calculate RHat statistics
  lambda_rhats = rep(0, length = data$Kp)
  score_rhats = matrix(0, nrow = data$I, ncol = data$Kp)
  fpc_rhats = matrix(0, nrow = data$D, ncol = data$Kp)
  sigma_rhat = Rhat(sigma_chains)
  mu_rhats = rep(0, length = data$D)
  smooth_rhats = rep(0, length = data$Kp + 1)
  
  for(m in 1:data$D){
    mu_rhats = Rhat(mu_chains[m,,])
  }
  for(k in 1:data$Kp){
    lambda_rhats[k] = Rhat(lambda_chains[k,,])
    smooth_rhats[k] = Rhat(smooth_chains[k,,])
    
    for(n in 1:data$I){
      score_rhats[n, k] = Rhat(score_chains[n,k,,])
    }
    
    for(m in 1:data$D){
      fpc_rhats[m, k] = Rhat(fpc_chains[m,k,,])
    }
  }
  smooth_rhats[data$Kp+1] = Rhat(smooth_chains[data$Kp,,])
  
  # Summarize and convert to dataframes
  fpc_names = paste0("FPC ", 1:data$Kp)
  
  score_rhats = data.frame(score_rhats)
  colnames(score_rhats) = fpc_names
  score_rhats = score_rhats %>%
    pivot_longer(cols = everything(), names_to = "FPC_Num", values_to = "Rhat") %>%
    group_by(FPC_Num) %>%
    summarize(Med_RHat = median(Rhat),
              Max_RHat = max(Rhat))
  
  func_rhats = data.frame(fpc_rhats)
  colnames(func_rhats) = fpc_names
  func_rhats = func_rhats %>%
    pivot_longer(cols = everything(), names_to = "Function", values_to = "Rhat") %>%
    group_by(Function) %>%
    summarize(Med_RHat = median(Rhat),
              Max_RHat = max(Rhat)) %>%
    rbind(data.frame(Function = "Mu", 
                     Med_RHat = median(mu_rhats),
                     Max_RHat = max(mu_rhats)))
  
  smooth_rhats = data.frame(Function = c(fpc_names, "Mu"), 
                            RHat = smooth_rhats)
  
  Variance_rhats = data.frame(Element = c(paste0("Lambda_", 1:data$Kp), "Sigma2"), 
                              RHat = c(lambda_rhats, sigma_rhat))
  
  outputs = list(Score = score_rhats, Func = func_rhats, 
                 Smoothing_Params = smooth_rhats, Variances = Variance_rhats)
  
  return(outputs)
}

# Function for assessing the convergence of single-level POLAR
RHat_POLAR <- function(model, data, align_EF){
  n_chains = length(model@inits)  
  n_samp = model@stan_args[[1]]$iter - model@stan_args[[1]]$warmup
  
  # Extract by chain
  by_chain = POLAR_byChain(jdh_mod, data_list$n, data_list$p, data_list$k)
  
  aligned_objs = map(1:n_chains, function(c){
    aligned = align_FPCs(by_chain$Phi[[c]], 
                         by_chain$Score[[c]],
                         anchor = align_EF)
    return(aligned)
  })
  
  score_chains = map(aligned_objs, function(x){
    out = x$Score %>% abind(along = 3)
    return(out)
  }) %>% abind(along = 4)
  
  lambda_chains = apply(score_chains, c(2,3,4), var)
  
  fpc_chains = map(aligned_objs, function(x){
    out = x$EF %>% abind(along = 3)
    return(out)
  }) %>% abind(along = 4)
  
  raw_samples = extract(model, permuted = F)
  dimen_names = dimnames(raw_samples)$parameters
  
  sig_idx = grepl("sig2", dimen_names)
  sigma_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      return(raw_samples[i,j,sig_idx])
    }) %>% abind(along = 1)
  }) %>% abind(along = 2)
  
  # Calculate RHat statistics
  lambda_rhats = rep(0, length = data$k)
  score_rhats = matrix(0, nrow = data$n, ncol = data$k)
  fpc_rhats = matrix(0, nrow = data$p, ncol = data$k)
  sigma_rhat = Rhat(sigma_chains)
  
  for(k in 1:data$k){
    lambda_rhats[k] = Rhat(lambda_chains[k,,])
    
    for(n in 1:data$n){
      score_rhats[n, k] = Rhat(score_chains[n,k,,])
    }
    
    for(m in 1:data$p){
      fpc_rhats[m, k] = Rhat(fpc_chains[m,k,,])
    }
  }
  
  # Summarize and convert to dataframes
  fpc_names = paste0("FPC ", 1:data$k)
  
  score_rhats = data.frame(score_rhats)
  colnames(score_rhats) = fpc_names
  score_rhats = score_rhats %>%
    pivot_longer(cols = everything(), names_to = "FPC_Num", values_to = "Rhat") %>%
    group_by(FPC_Num) %>%
    summarize(Med_RHat = median(Rhat),
              Max_RHat = max(Rhat))
  
  func_rhats = data.frame(fpc_rhats)
  colnames(func_rhats) = fpc_names
  func_rhats = func_rhats %>%
    pivot_longer(cols = everything(), names_to = "Function", values_to = "Rhat") %>%
    group_by(Function) %>%
    summarize(Med_RHat = median(Rhat),
              Max_RHat = max(Rhat))
  
  Variance_rhats = data.frame(Element = c(paste0("Lambda_", 1:data$k), "Sigma2"), 
                              RHat = c(lambda_rhats, sigma_rhat))
  
  outputs = list(Score = score_rhats, Func = func_rhats, Variances = Variance_rhats)
  
  return(outputs)
}

