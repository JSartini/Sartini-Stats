# Convert FPC matrix to dataframe
FPC_df <- function(EF_mat, Domain, P = NULL){
  out_df = data.frame(EF_mat)
  colnames(out_df) = paste0("FPC ", 1:ncol(out_df))
  if(!is.null(P)){
    out_df$Arg = rep(Domain, P)
    out_df$Var = rep(1:P, each = length(Domain))
    out_df = out_df %>%
      pivot_longer(-c(Var, Arg), names_to = "FPC_Num", values_to = "FPC_Val") 
  }
  else{
    out_df$Arg = Domain
    out_df = out_df %>%
      pivot_longer(-c(Arg), names_to = "FPC_Num", values_to = "FPC_Val") 
  }
  return(out_df)
}

# Align FPC samples using Procrustes analysis
procrust_FPC <- function(EF_list, anchor = NULL){
  if(is.null(anchor)){
    anchor = EF_list[[1]]
  }
  
  out_EF = list()
  
  K = ncol(anchor)
  for(i in 1:length(EF_list)){
    Phi_s = EF_list[[i]]
    
    R_comps = svd(t(Phi_s) %*% anchor)
    R = R_comps$u %*% t(R_comps$v)
    
    out_EF[[i]] = Phi_s %*% R
  }
  
  return(out_EF)
}

# Align FPC weight samples using Procrustes analysis
procrust_WEI <- function(Weight_list, B, P, anchor = NULL, Score_list = NULL){
  kB = kronecker(diag(P), B)
  
  if(is.null(anchor)){
    anchor = kB %*% Weight_list[[1]]
  }
  
  EF = map(Weight_list, function(Psi_s){
    return(Phi_s = kB %*% Psi_s)
  })
  
  Rotations = map(EF, function(Phi_s){
    R_comps = svd(t(Phi_s) %*% anchor)
    R = R_comps$u %*% t(R_comps$v)
    return(R)
  })
  
  out_Weight = map2(Weight_list, Rotations, function(Psi_s, R){
    return(Psi_s %*% R)
  })
  out_EF = map2(EF, Rotations, function(Phi_s, R){
    return(Phi_s %*% R)
  })
  if(!is.null(Score_list)){
    out_Score = map2(Score_list, Rotations, function(Xi_s, R){
      return(Xi_s %*% R)
    }) 
  }
  else{
    out_Score = NULL
  }
  
  return(list(Weights = out_Weight, EF = out_EF, Score = out_Score))
}

# Calculate FPC estimate by taking the right SV of mean Phi - Goldsmith 2015
FPC_Est_RSV <- function(FPCs, M, P, Domain, anchor){
  phi_hat_uncon = abind(FPCs, along = 3) %>%
    apply(c(1,2), mean)
  phi_objs = svd(phi_hat_uncon)
  phi_est = phi_objs$u * sqrt(M) # Re-scale
  
  if(!is.null(anchor)){
    R_comps = svd(t(phi_est) %*% anchor)
    R = R_comps$u %*% t(R_comps$v)
    phi_est = phi_est %*% R
  }
  
  return(FPC_df(phi_est, Domain, P))
}

# Take mean of orthonormal weights, then orthogonalize
FPC_Est_WEI <- function(weights, B, P, Domain, anchor = NULL){
  kB = kronecker(diag(1, P), B)
  psi_hat = abind(weights, along = 3) %>%
    apply(c(1,2), mean)
  psi_svd = svd(psi_hat)
  psi_est = psi_svd$u %*% t(psi_svd$v)
  phi_est = kB %*% psi_est
  
  if(!is.null(anchor)){
    R_comps = svd(t(phi_est) %*% anchor)
    R = R_comps$u %*% t(R_comps$v)
    psi_fin = psi_est %*% R
    phi_est = kB %*% psi_fin 
  }
  
  return(FPC_df(phi_est, Domain, P))
}

# Calculate CI of FPC using vectorized operations
FPC_CI <- function(EF_list, Domain, P){
  n_samp = length(EF_list)
  
  EF_mat = EF_list %>% abind(along = 0)
  
  EF_LB = FPC_df(apply(EF_mat, c(2,3), quantile, probs = c(0.025)), Domain, P) %>%
    rename(LB = FPC_Val)
  EF_UB = FPC_df(apply(EF_mat, c(2,3), quantile, probs = c(0.975)), Domain, P) %>%
    rename(UB = FPC_Val)
  
  return(inner_join(EF_LB, EF_UB, by = c("Var", "Arg", "FPC_Num")))
}

# Calculate FE estimates with associated CI using vectorized operations
FE_Summary <- function(Mu_list, Domain, P){
  fe_df <- function(Mu_vec, Domain, P){
    return(data.frame(Mu = Mu_vec, 
                      Arg = rep(Domain, P),
                      Var = rep(1:P, each = length(Domain))))
  }
  
  fe_mat = Mu_list %>%
    abind(along = 0)
  
  fe_est = fe_df(apply(fe_mat, 2, mean), Domain, P) %>%
    rename(Est = Mu)
  fe_lb = fe_df(apply(fe_mat, 2, quantile, probs = c(0.025)), Domain, P) %>%
    rename(LB = Mu)
  fe_ub = fe_df(apply(fe_mat, 2, quantile, probs = c(0.975)), Domain, P) %>%
    rename(UB = Mu)
  
  fe_out = fe_est %>%
    inner_join(fe_lb, by = c("Var", "Arg")) %>%
    inner_join(fe_ub, by = c("Var", "Arg"))
  
  return(fe_out)
}

# Calculate smooth estimates with associated CI using vectorized operations
Smooth_Summary <- function(N, P, Mu_list, EF_list, Score_list, Domain){
  sm_df <- function(mat, ids, P, Domain){
    out_df = data.frame(mat)
    colnames(out_df) = paste0("Curve ", ids)
    out_df$Arg = rep(Domain, P)
    out_df$Var = rep(1:P, each = length(Domain))
    out_df = out_df %>%
      pivot_longer(-c(Var, Arg), names_to = "Curve", values_to = "Smooth")
  }
  
  n_samp = length(Mu_list)
  size_group = 25
  sub_groups = split(1:N, ceiling(seq_along(1:N)/size_group))
  
  plan(multisession, workers = 4)
  Smooth_DF = future_map(sub_groups, function(idxs){
    sg = length(idxs)
    
    smooths = map(1:n_samp, function(x){
      smooth = matrix(rep(Mu_list[[x]], sg), nrow = sg, byrow = T)
      smooth = smooth + Score_list[[x]][idxs,] %*% t(EF_list[[x]])
      return(smooth)
    }) %>% abind(along = 0)
    
    mean_sm = sm_df(t(apply(smooths, c(2,3), mean)), idxs, P, Domain) %>%
      rename(Est = Smooth)
    lb_sm = sm_df(t(apply(smooths, c(2,3), quantile, probs = c(0.025))), 
                         idxs, P, Domain) %>%
      rename(LB = Smooth)
    ub_sm = sm_df(t(apply(smooths, c(2,3), quantile, probs = c(0.975))), 
                         idxs, P, Domain) %>%
      rename(UB = Smooth)
    
    out_df = mean_sm %>%
      inner_join(lb_sm, by = c("Curve", "Var", "Arg")) %>%
      inner_join(ub_sm, by = c("Curve", "Var", "Arg"))
    
    return(out_df)}) %>% list_rbind()
  plan(sequential)
  return(Smooth_DF)
}

# Form smooth samples and return in matrix format
Smooth_Raw <- function(Mu_list, EF_list, Score_list, data_list){
  n_samp = length(Mu_list)
  N = data_list$N
  P = data_list$P
  M = data_list$M
  
  offsets_mu = data_list$consts %>%
    arrange(Var) %>%
    pull(mu_Y) %>%
    rep(each = M)
  offsets_mu = matrix(rep(offsets_mu, N), nrow = N, byrow = T)
  
  scales_sd = data_list$consts %>%
    arrange(Var) %>%
    pull(sd_Y) %>%
    rep(each = M)
  scales_sd = matrix(rep(scales_sd, N), nrow = N, byrow = T)
  
  plan(multisession, workers = 4)
  Smooth_mat = future_map(1:n_samp, function(x){
    smooth = matrix(rep(Mu_list[[x]], N), nrow = N, byrow = T)
    smooth = smooth + Score_list[[x]] %*% t(EF_list[[x]])
    smooth = smooth * scales_sd + offsets_mu
    return(smooth)
  }) %>% abind(along = 0)
  plan(sequential)
  
  return(Smooth_mat)
}
