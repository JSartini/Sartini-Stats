# Convert FPC matrix to dataframe
FPC_df <- function(EF_mat, Domain){
  out_df = data.frame(EF_mat)
  colnames(out_df) = paste0("FPC ", 1:ncol(out_df))
  out_df$Arg = Domain
  out_df = out_df %>%
    pivot_longer(-c(Arg), names_to = "FPC_Num", values_to = "FPC_Val")
  return(out_df)
}

# Align FPC samples using sign flips and re-ordering
align_weights <- function(Weight_list, Score_list, B, anchor = NULL){
  if(is.null(anchor)){
    evals = apply(Score_list[[1]], 2, var)
    fpc_order = sort(evals, decreasing = T, index.return = T)$ix
    anchor = (B %*% Weight_list[[1]])[,fpc_order]
  }
  
  out_Weight = list()
  out_EF = list()
  out_Score = list()
  
  K = ncol(anchor)
  for(i in 1:length(Weight_list)){
    evals = apply(Score_list[[i]], 2, var)
    fpc_order = sort(evals, decreasing = T, index.return = T)$ix
    
    Psi = Weight_list[[i]][,fpc_order]
    Xi = Score_list[[i]][,fpc_order]
    Phi = B %*% Psi
    
    for(k in 1:K){
      if(sum(Phi[,k]*anchor[,k]) < 0){
        Psi[,k] = -Psi[,k]
        Xi[,k] = -Xi[,k]
      }
    }
    
    out_Weight[[i]] = Psi
    out_EF[[i]] = B %*% Psi
    out_Score[[i]] = Xi
  }
  
  return(list(Weights = out_Weight, EF = out_EF, Score = out_Score))
}

align_FPCs <- function(EF_list, Score_list, anchor = NULL){
  if(is.null(anchor)){
    evals = apply(Score_list[[1]], 2, var)
    fpc_order = sort(evals, decreasing = T, index.return = T)$ix
    anchor = EF_list[[1]][,fpc_order]
  }
  
  out_EF = list()
  out_Score = list()
  
  K = ncol(anchor)
  for(i in 1:length(EF_list)){
    evals = apply(Score_list[[i]], 2, var)
    fpc_order = sort(evals, decreasing = T, index.return = T)$ix
    
    Phi = EF_list[[i]][,fpc_order]
    Xi = Score_list[[i]][,fpc_order]
    for(k in 1:K){
      if(sum(Phi[,k]*anchor[,k]) < 0){
        Phi[,k] = -Phi[,k]
        Xi[,k] = -Xi[,k]
      }
    }
    
    out_EF[[i]] = Phi
    out_Score[[i]] = Xi
  }
  
  return(list(EF = out_EF, Score = out_Score))
}

# Align chains separately, then collate - simple method
chain_align <- function(EF_list, Score_list, anchor = NULL){
  if(is.null(anchor)){
    evals = apply(Score_list[[1]][[1]], 2, var)
    fpc_order = sort(evals, decreasing = T, index.return = T)$ix
    anchor = EF_list[[1]][[1]][,fpc_order]
  }
  
  out_EF = list()
  out_Score = list()
  n_chain = length(EF_list)
  for(i in 1:n_chain){
    align_chain = full_align_simple(EF_list[[i]], Score_list[[i]], anchor)
    out_EF[[i]] = align_chain$EF
    out_Score[[i]] = align_chain$Score
  }
  
  return(list(EF = out_EF, Score = out_Score))
}

# Calculate CI of FPC in dataframe format
CI_EF <- function(EF_list, Domain){
  n_samp = length(EF_list)
  
  EF_df = map(1:n_samp, function(idx){
    fpc_mat = EF_list[[idx]]
    out_df = FPC_df(fpc_mat, Domain)
    out_df$Sample = idx
    return(out_df)
  }) %>% 
    list_rbind() %>%
    group_by(Arg, FPC_Num) %>%
    summarize(LB = quantile(FPC_Val, probs = c(0.025)), 
              UB = quantile(FPC_Val, probs = c(0.975)))
  
  return(EF_df)
}

# Calculate CI and estimates of scores in dataframe format
out_Score <- function(Score_list){
  n_samp = length(Score_list)
  
  score_df = map(1:n_samp, function(idx){
    smat = Score_list[[idx]]
    out_df = data.frame(t(smat))
    colnames(out_df) = paste0("Curve ", 1:ncol(out_df))
    out_df$FPC = paste0("FPC ", 1:nrow(out_df))
    out_df = out_df %>%
      pivot_longer(-c(FPC), names_to = "Curve", values_to = "Score")
    out_df$Sample = idx
    return(out_df)
  }) %>% 
    list_rbind() %>%
    group_by(FPC, Curve) %>%
    summarize(Est = mean(Score), 
              LB = quantile(Score, probs = c(0.025)), 
              UB = quantile(Score, probs = c(0.975)))
  
  return(score_df)
}

# Calculate FE dataframe
out_FE <- function(Mu_list, Domain){
  n_sample = length(Mu_list)
  
  fe_df = map(1:n_sample, function(x){
    fe_sample = data.frame(Mu = Mu_list[[x]], 
                           Arg = Domain, Sample = x)
    return(fe_sample)
  }) %>% list_rbind()
  
  return(fe_df)
}

# Calculate GP deviations using FPCs and associated scores
out_Devs <- function(EF_list, Score_list){
  n_samp = length(EF_list)
  
  smooths = map(1:n_samp, function(x){
    delta = Score_list[[x]] %*% t(EF_list[[x]])
    return(delta)
  })
  
  return(smooths)
}

# Calculate smooths for (M)FPCA
out_Smooths <- function(N, Mu_list, EFL_list, ScoreL_list, Domain){
  n_samp = length(Mu_list)
  n_levels = length(EFL_list)
  
  devs = list()
  for(i in 1:n_levels){
    devs[[i]] = out_Devs(EFL_list[[i]], ScoreL_list[[i]])
  }
  
  smooths = map(1:n_samp, function(x){
    smooth = matrix(rep(Mu_list[[x]], N), nrow = N, byrow = T)
    for(i in 1:n_levels){
      smooth = smooth + devs[[i]][[x]]
    }
    out_df = data.frame(t(smooth))
    colnames(out_df) = paste0("Curve ", 1:ncol(out_df))
    out_df$Arg = Domain
    out_df$Sample = x
    out_df = out_df %>%
      pivot_longer(-c(Arg, Sample), names_to = "Curve", values_to = "Smooth")
    
    return(out_df)
  }) %>% list_rbind() %>%
    group_by(Arg, Curve) %>%
    summarize(Est = mean(Smooth), 
              LB = quantile(Smooth, probs = c(0.025)), 
              UB = quantile(Smooth, probs = c(0.975)))
  
  return(smooths)
}

# Calculate FPC estimate by taking right SV of mean smooth SVD - Jauch
Smooth_SVD_FPC <- function(devs, Domain, K, anchor){
  M = length(Domain)
  mean_dev = abind(devs, along = 3) %>%
    apply(c(1,2), mean)
  
  FPC_est = svd(mean_dev, nv = K)$v*sqrt(M)
  for(k in 1:K){
    if(sum(FPC_est[,k]*anchor[,k]) < 0){
      FPC_est[,k] = -FPC_est[,k]
    }
  }
  return(FPC_df(FPC_est, Domain))
}

# Calculate FPC estimate by taking the right SV of mean Phi - Goldsmith
Phi_SVD_FPC <- function(FPCs, Domain, K, anchor){
  M = length(Domain)
  phi_hat_uncon = abind(FPCs, along = 3) %>%
    apply(c(1,2), mean)
  phi_objs = svd(phi_hat_uncon)
  FPC_est = phi_objs$u %*% t(phi_objs$v) * sqrt(M)
  for(k in 1:K){
    if(sum(FPC_est[,k]*anchor[,k]) < 0){
      FPC_est[,k] = -FPC_est[,k]
    }
  }
  return(FPC_df(FPC_est, Domain))
}

# Take mean of orthonormal weights, then orthogonalize - Novel
Psi_SVD_FPC <- function(weights, B, Domain, K, anchor = NULL){
  psi_hat = abind(weights, along = 3) %>%
    apply(c(1,2), mean)
  psi_svd = svd(psi_hat)
  FPC_est = B %*% psi_svd$u %*% t(psi_svd$v)
  if(!is.null(anchor)){
    for(k in 1:K){
      if(sum(FPC_est[,k]*anchor[,k]) < 0){
        FPC_est[,k] = -FPC_est[,k]
      }
    } 
  }
  return(FPC_df(FPC_est, Domain))
}

# Form smooth samples and return in matrix format
Smooth_Raw <- function(Mu_list, EF_list, Score_list, data_list){
  n_samp = length(Mu_list)
  N = data_list$N
  
  plan(multisession, workers = 4)
  Smooth_mat = future_map(1:n_samp, function(x){
    smooth = matrix(rep(Mu_list[[x]], N), nrow = N, byrow = T)
    smooth = smooth + Score_list[[x]] %*% t(EF_list[[x]])
    return(smooth)
  }) %>% abind(along = 0)
  plan(sequential)
  
  return(Smooth_mat)
}

