sample_scores <- function(Yi, Bi, samples){
  P = length(Yi)
  Q = ncol(Bi[[1]])
  
  plan(multisession, workers = 4)
  samples_df = future_map(1:dim(samples$sigma2)[1], function(x){
    M = 0
    V_comp = 0
    for(p in 1:P){
      sdx = (p-1)*Q+1
      edx = p*Q
        
      Ri = Yi[[p]] - Bi[[p]] %*% samples$w_mu[x,sdx:edx]
      DH = Bi[[p]] %*% samples$Psi[x,sdx:edx,]
      
      M = M + t(Ri) %*% DH / samples$sigma2[x,p]
      V_comp = V_comp + t(DH) %*% DH / samples$sigma2[x,p]
    }
    
    V_inv = V_comp + diag(1/rev(samples$lambda[x,]))
    V = chol2inv(chol(V_inv))
    
    score_samples = mvrnorm(n = 1, mu = as.vector(V %*% t(M)), Sigma = V)
    
    output = data.frame(Scores = score_samples, 
                        SNames = paste0("xi", 1:length(score_samples)), 
                        Sample = x) %>%
      pivot_wider(names_from = SNames, values_from = Scores)
    
    return(output)
  }, .options = furrr_options(seed = T)) %>% list_rbind()
  plan(sequential)
  
  return(samples_df)
}

