# Function for assessing the convergence of single-level FAST
RHat_FAST <- function(model, data, B, align_EF){
  n_chains = length(model@inits)  
  n_samp = model@stan_args[[1]]$iter - model@stan_args[[1]]$warmup
  kB = kronecker(diag(data$P), B)
  
  # Extract by chain
  by_chain = FAST_byChain(model, data$N, data$M, data$Q, data$K, data$P)
  
  aligned_objs = map(1:n_chains, function(c){
    aligned = procrust_WEI(by_chain$Psi[[c]], B,
                           data$P, anchor = align_EF, 
                           Score_list = by_chain$Score[[c]])
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
  
  musmooth_idx = grepl("h_mu", dimen_names)
  mu_smooth_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      mean_smooth = raw_samples[i,j,musmooth_idx]
      return(mean_smooth)
    }) %>% abind(along = 2)
  }) %>% abind(along = 3)
  
  fpcsmooth_idx = grepl("^H", dimen_names)
  fpc_smooth_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      fpc_smooth = raw_samples[i,j,fpcsmooth_idx]
      dim(fpc_smooth) = c(data$P, data$K) # 1
      return(fpc_smooth)
    }) %>% abind(along = 3)
  }) %>% abind(along = 4)
  
  sig_idx = grepl("sigma2", dimen_names)
  sigma_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      return(raw_samples[i,j,sig_idx])
    }) %>% abind(along = 1)
  }) %>% abind(along = 2)
  
  mu_idx = grepl("w_mu", dimen_names)
  mu_chains = map(1:n_chains, function(j){
    output = map(1:n_samp, function(i){
      mean_func = kB %*% raw_samples[i,j,mu_idx]
      return(mean_func)
    }) %>% abind(along = 2)
  }) %>% abind(along = 3)
  
  # Calculate RHat statistics
  lambda_rhats = rep(0, length = data$K)
  score_rhats = matrix(0, nrow = data$N, ncol = data$K)
  fpc_rhats = matrix(0, nrow = data$P*data$M, ncol = data$K)
  sigma_rhat = Rhat(sigma_chains)
  mu_rhats = rep(0, length = data$P*data$M)
  mu_smooth_rhats = rep(0, length = data$P)
  fpc_smooth_rhats = matrix(0, nrow = data$P, ncol = data$K) # 1
  
  for(p in 1:data$P){
    mu_smooth_rhats[p] = Rhat(mu_smooth_chains[p,,])
  }
  for(m in 1:(data$P*data$M)){
    mu_rhats = Rhat(mu_chains[m,,])
  }
  for(k in 1:data$K){
    lambda_rhats[k] = Rhat(lambda_chains[k,,])
    
    for(p in 1:data$P){
      fpc_smooth_rhats[p,k] = Rhat(fpc_smooth_chains[p,k,,])
    }
    
    # fpc_smooth_rhats[1,k] = Rhat(fpc_smooth_chains[1,k,,])
    
    for(n in 1:data$N){
      score_rhats[n, k] = Rhat(score_chains[n,k,,])
    }
    
    for(m in 1:(data$P*data$M)){
      fpc_rhats[m, k] = Rhat(fpc_chains[m,k,,])
    }
  }
  
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
  
  mu_smooth_rhats = data.frame(Function = paste0("Mu", 1:data$P), 
                               RHat = mu_smooth_rhats)
  
  fpc_smooth_rhats = data.frame(fpc_smooth_rhats)
  colnames(fpc_smooth_rhats) = paste0("FPC ", 1:data$K)
  fpc_smooth_rhats$Var = paste0("Var ", 1:data$P)
  fpc_smooth_rhats = fpc_smooth_rhats %>%
    pivot_longer(-c(Var), names_to = "FPC", values_to = "RHat") # everything(),
  
  Variance_rhats = data.frame(Element = c(paste0("Lambda_", 1:data$K), "Sigma2"), 
                              RHat = c(lambda_rhats, sigma_rhat))
  
  outputs = list(Score = score_rhats, Func = func_rhats, Mu_Smoothing = mu_smooth_rhats, 
                 FPC_Smoothing = fpc_smooth_rhats, Variances = Variance_rhats)
  
  return(outputs)
}
