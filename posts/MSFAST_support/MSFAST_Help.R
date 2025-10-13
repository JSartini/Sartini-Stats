# Produce P_alpha using functional basis - numerical integration
FAST_P <- function(basis, derivs, Q){
  P0 = matrix(0, ncol = Q, nrow = Q)
  P2 = matrix(0, ncol = Q, nrow = Q)
  # Adjust desired accuracy as problem becomes more complex due to large # bases
  if(Q <= 30){
    des.tol = .Machine$double.eps^0.30
  }
  else{
    des.tol = .Machine$double.eps^0.40
  }
  for(i in 1:Q){
    for(j in 1:i){
      # Rescale for numerical precision purposes
      f0 <- function(x){ return(basis[[i]](x) * basis[[j]](x)) } 
      f2 <- function(x){ return(derivs[[i]](x) * derivs[[j]](x)) }
      P0[i,j] = stats::integrate(f0, lower = 0, upper = 1, subdivisions = 10000, 
                                 rel.tol = des.tol)$value
      P2[i,j] = stats::integrate(f2, lower = 0, upper = 1, subdivisions = 10000, 
                                 rel.tol = des.tol)$value
      if(i != j){
        P0[j,i] = P0[i,j]
        P2[j,i] = P2[i,j]
      }
    }
  } 
  return(list(P0 = P0, P2 = P2))
}

# Evaluate the orthogonal spline basis at specified time points
FAST_B <- function(basis_type = "B", Q, Domain){
  if(basis_type == "B"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = Q-2), rep(1, 3)))
    B = evaluate(B_f, Domain)
  }
  else{
    if(basis_type == "Fourier"){
      basis = Fourier_bases(Q)
    }
    else if(basis_type == "Legendre"){
      basis = Legendre_bases(Q)
    }
    else if(basis_type == "Splinet"){
      basis = Splinet_bases(Q)$B
    }
    else{
      stop("Basis not supported")
    }
    B = map(1:Q, function(q){
      return(basis[[q]](Domain))
    }) %>% abind(along = 2)
  }
  return(B)
}

# Produce inputs for FAST
MSFAST_datalist <- function(df, N, K, Q, Domain, basis_type = "B",
                            alpha = 0.1, scale = T){
  
  # Generate spline basis with derivatives
  if(basis_type == "B"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = Q-2), rep(1, 3)))
    B = evaluate(B_f, Domain)
    P0 = GramMatrix(B_f)
    P2 = OuterProdSecondDerivative(B_f)
    P_alpha = alpha * P0 + (1-alpha) * P2
  }
  else{
    B = FAST_B(basis_type, Q, Domain)
    if(basis_type == "Fourier"){
      basis = Fourier_bases(Q)
      derivs = Fourier_d2(Q)
    }
    else if(basis_type == "Splinet"){
      bobj = Splinet_bases(Q)
      basis = bobj$B
      derivs = Splinet_d2(Q, bobj$cInt, bobj$cSlo)
    }
    else if(basis_type == "Legendre"){
      basis = Legendre_bases(Q)
      derivs = Legendre_d2(Q)
    }
    else{
      stop("Basis not supported")
    }
    P_est = FAST_P(basis, derivs, Q)
    P_alpha = alpha * P_est$P0 + (1-alpha) * P_est$P2
  }
  
  # Scale the Y-values, storing the location and scale
  if(scale){
    const_df = df %>%
      group_by(Var) %>%
      summarize(mu_Y = mean(Y), 
                sd_Y = sd(Y))
  }
  else{
    const_df = df %>%
      group_by(Var) %>%
      summarize(mu_Y = 0, 
                sd_Y = 1)
  }
  arranged_df = df %>%
    arrange(Var, Arg) %>%
    left_join(const_df, by = "Var") %>%
    mutate(Y = (Y - mu_Y)/sd_Y) %>%
    select(-c(mu_Y, sd_Y))  
  
  TPC = arranged_df %>%
    group_by(Var) %>%
    summarize(m_c = n()) %>%
    pull(m_c)
  dim(TPC) = c(P)
  
  M_val = length(Domain)
  
  # Format for FAST
  fast_list = list(N = N, M = M_val, L = nrow(arranged_df), Q = Q, K = K, P = length(TPC), 
                   Y = arranged_df$Y, Subj = arranged_df$Subj, S = arranged_df$S,
                   Tp_card = TPC, B = B, P_alpha = P_alpha, consts = const_df)
  return(fast_list)
}

# Place EF and scores in standard form for post-processing
FAST_extract <- function(mod_fit, B, DL){
  samples = extract(mod_fit)
  n_samp = dim(samples$sigma2)[1]
  kB = kronecker(diag(DL$P), B)
  
  Weight_list = map(1:n_samp, function(x){
    return(samples$Psi[x,,])
  })
  
  EF_list = map(Weight_list, function(weights){
    EF_sample = kB %*% weights
    return(EF_sample)
  })
  
  Score_list = map(1:n_samp, function(x){
    Score_sample = samples$Scores[x,,]
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return((kB %*% samples$w_mu[x,]))
  })
  
  return(list(Weights = Weight_list, EF = EF_list, 
              Score = Score_list, Mu = Mu_list))
}

# Extract EF and Scores by-chain
FAST_byChain <- function(mod_fit, N, M, Q, K, P){
  chain_samples = extract(mod_fit, permuted = F)
  dimen_names = dimnames(chain_samples)$parameters
  
  n_samp = dim(chain_samples)[1]
  n_chain = dim(chain_samples)[2]
  
  Psi_idx = grepl("Psi", dimen_names)
  Score_idx = grepl("Scores", dimen_names)
  
  Psi_byChain = list()
  Score_byChain = list()
  
  for(j in 1:n_chain){
    Psi_byChain[[j]] = map(1:n_samp, function(i){
      Psi = chain_samples[i,j,Psi_idx]
      dim(Psi) = c(P*Q, K)
      return(Psi)
    })
    Score_byChain[[j]] = map(1:n_samp, function(i){
      Xi = chain_samples[i,j,Score_idx]
      dim(Xi) = c(N, K)
      return(Xi)
    })
  }
  return(list(Psi = Psi_byChain, Score = Score_byChain))
}
