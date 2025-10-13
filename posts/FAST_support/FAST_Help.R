# Produce P_alpha using functional basis - numerical integration
FAST_P <- function(derivs, Q, upp){
  P2 = matrix(0, ncol = Q, nrow = Q)
  for(i in 1:Q){
    for(j in 1:i){
      # Rescale for numerical precision purposes
      f2 <- function(x){ return(derivs[[i]](x) * derivs[[j]](x)) }
      P2[i,j] = stats::integrate(f2, lower = 0, upper = upp, 
                                 subdivisions = 10000)$value
      if(i != j){
        P2[j,i] = P2[i,j]
      }
    }
  } 
  return(P2)
}

FAST_B <- function(basis_type = "Fourier", Q, Domain){
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

# Produce data/hyper-priors for FAST
FAST_datalist <- function(Y, Q, K, Domain, basis_type = "Splinet", 
                          alpha = 0.1, scale = T){
  
  # Generate FE spline basis with penalty
  spline_knots = c(rep(0, 3), seq(0, 1, length.out = Q-2), rep(1, 3))
  B_f = SplineBasis(spline_knots)
  B_FE = evaluate(B_f, Domain)
  P0 = GramMatrix(B_f)
  P2 = OuterProdSecondDerivative(B_f)
  P_FE = alpha * P0 + (1-alpha) * P2
  
  # Generate RE spline basis with penalty
  P0 = diag(Q)
  if(basis_type == "B"){
    B_f = OBasis(spline_knots)
    B_RE = evaluate(B_f, Domain)
    P2 = OuterProdSecondDerivative(B_f)
    P_RE = alpha * P0 + (1-alpha) * P2
  }
  else{
    B_RE = FAST_B(basis_type, Q, Domain)
    upp = 1
    if(basis_type == "Fourier"){
      derivs = Fourier_d2(Q)
    }
    else if(basis_type == "Splinet"){
      bobj = Splinet_bases(Q)
      derivs = Splinet_d2(Q, bobj$cInt, bobj$cSlo)
      upp = 10
    }
    else if(basis_type == "Legendre"){
      derivs = Legendre_d2(Q)
    }
    else{
      stop("Basis not supported")
    }
    P2 = FAST_P(derivs, Q, upp)
    P2 = P2/eigen(P2)$values[1]
    P_RE = alpha * P0 + (1-alpha) * P2
  }
  
  # Scale the Y-values, storing the location and scale
  if(scale){
    sd_Y = sd(Y)
    mu_Y = mean(Y) 
  }
  else{
    sd_Y = 1
    mu_Y = 0
  }
  scaled_Y = (Y - mu_Y)/sd_Y
  M_val = length(Domain)
  
  # Format for FAST
  fast_list = list(N = nrow(Y), M = M_val, Q = Q, K = K, Y = scaled_Y, 
                   B_FE = B_FE, P_FE = P_FE, B_RE = B_RE, P_RE = P_RE,
                   sd_Y = sd_Y, mu_Y = mu_Y)
  return(fast_list)
}

FAST_DL_ML <- function(Y, N, IDs, K1, K2, Q, Domain, basis_type = "Fourier", 
                        alpha = 0.1, scale = T){
  
  # Generate spline basis with derivatives
  if(basis_type == "B"){
    B_f = OBasis(c(rep(0, 3), seq(0, 1, length.out = Q-2), rep(1, 3)))
    B = evaluate(B_f, Domain)
    P0 = GramMatrix(B_f)
    P2 = OuterProdSecondDerivative(B_f)
    P_alpha = alpha * P0 + (1-alpha) * P2
  }
  else if(basis_type == "Approx"){ 
    M = length(Domain)
    
    B_f = bSpline(Domain, df = Q, intercept = T)
    B = eigen(B_f %*% t(B_f))$vectors[,1:Q]
    sizes = diag(t(B) %*% diag(W) %*% B)
    B = B %*% diag(sqrt(1/sizes))
    
    # Calculate penalty matrix
    diff0 = diag(1, M, M)
    diff2 = matrix(rep(c(1,-2,1, rep(0, M-2)), M-2)[1:((M-2)*M)], M-2, M, byrow = TRUE)
    P0 = t(B) %*% t(diff0) %*% diff0 %*% B
    P2 = t(B) %*% t(diff2) %*% diff2 %*% B
    P_alpha = alpha * P0 + (1-alpha) * P2
  }
  else{
    B = FAST_B(basis_type, Q, Domain)
    P0 = diag(Q)
    upp = 1
    if(basis_type == "Fourier"){
      derivs = Fourier_d2(Q)
    }
    else if(basis_type == "Legendre"){
      derivs = Legendre_d2(Q)
    }
    else if(basis_type == "Splinet"){
      bobj = Splinet_bases(Q)
      derivs = Splinet_d2(Q, bobj$cInt, bobj$cSlo)
      upp = 10
    }
    else{
      stop("Basis not supported")
    }
    P2 = FAST_P(derivs, Q, upp)
    P2 = P2/eigen(P2)$values[1]
    P_alpha = alpha * P0 + (1-alpha) * P2
  }
  
  # Scale the Y-values, storing the location and scale
  if(scale){
    sd_Y = sd(Y)
    mu_Y = mean(Y) 
  }
  else{
    sd_Y = 1
    mu_Y = 0
  }
  scaled_Y = (Y - mu_Y)/sd_Y
  M_val = length(Domain)
  
  # Format for FAST
  fast_list = list(N = N, I = n_distinct(IDs), ID = IDs, M = M_val, Q = Q, 
                   K1 = K1, K2 = K2, Y = scaled_Y, B = B, P_alpha = P_alpha, 
                   sd_Y = sd_Y, mu_Y = mu_Y)
  return(fast_list)
}

# Place EF and scores in standard form for post-processing
FAST_extract <- function(mod_fit, B_FE, B_RE, Domain, DL){
  samples = extract(mod_fit)
  n_samp = length(samples$sigma2)
  # M = length(Domain)
  
  Weight_list = map(1:n_samp, function(x){
    return(samples$Psi[x,,])
  })
  
  EF_list = map(Weight_list, function(weights){
    EF_sample = B_RE %*% weights
    return(EF_sample)
  })
  
  Score_list = map(1:n_samp, function(x){
    Score_sample = samples$Scores[x,,] * DL$sd_Y
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return((B_FE %*% samples$w_mu[x,]) * DL$sd_Y + DL$mu_Y)
  })
  
  return(list(Weights = Weight_list, EF = EF_list, 
              Score = Score_list, Mu = Mu_list))
}

FAST_ML_extract <- function(mod_fit, B, Domain, DL){
  samples = extract(mod_fit)
  n_samp = length(samples$sigma2)
  # M = length(Domain)
  
  Weight1_list = map(1:n_samp, function(x){
    return(samples$Psi_1[x,,])
  })
  
  EF1_list = map(Weight1_list, function(weights){
    EF_sample = B %*% weights
    return(EF_sample)
  })
  
  Score1_list = map(1:n_samp, function(x){
    Score_sample = samples$Scores1[x,,] * DL$sd_Y
    return(Score_sample)
  })
  
  Weight2_list = map(1:n_samp, function(x){
    return(samples$Psi_2[x,,])
  })
  
  EF2_list = map(Weight2_list, function(weights){
    EF_sample = B %*% weights
    return(EF_sample)
  })
  
  Score2_list = map(1:n_samp, function(x){
    Score_sample = samples$Scores2[x,,] * DL$sd_Y
    return(Score_sample)
  })
  
  Mu_list = map(1:n_samp, function(x){
    return((B %*% samples$w_mu[x,]) * DL$sd_Y + DL$mu_Y)
  })
  
  return(list(W1 = Weight1_list, EF1 = EF1_list, 
              W2 = Weight2_list, EF2 = EF2_list,
              S1 = Score1_list, S2 = Score2_list,
              Mu = Mu_list))
}

# Extract EF and Scores by-chain
FAST_byChain <- function(mod_fit, N, M, K){
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
      dim(Psi) = c(Q, K)
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

FAST_ML_byChain <- function(mod_fit, N, I, M, K1, K2){
  chain_samples = extract(mod_fit, permuted = F)
  dimen_names = dimnames(chain_samples)$parameters
  
  n_samp = dim(chain_samples)[1]
  n_chain = dim(chain_samples)[2]
  
  Psi1_idx = grepl("Psi_1", dimen_names)
  Score1_idx = grepl("Scores1", dimen_names)
  Psi1_byChain = list()
  Score1_byChain = list()
  
  Psi2_idx = grepl("Psi_2", dimen_names)
  Score2_idx = grepl("Scores2", dimen_names)
  Psi2_byChain = list()
  Score2_byChain = list()
  
  for(j in 1:n_chain){
    Psi1_byChain[[j]] = map(1:n_samp, function(i){
      Psi = chain_samples[i,j,Psi1_idx]
      dim(Psi) = c(Q, K1)
      return(Psi)
    })
    Score1_byChain[[j]] = map(1:n_samp, function(i){
      Xi = chain_samples[i,j,Score1_idx]
      dim(Xi) = c(I, K1)
      return(Xi)
    })
    Psi2_byChain[[j]] = map(1:n_samp, function(i){
      Psi = chain_samples[i,j,Psi2_idx]
      dim(Psi) = c(Q, K2)
      return(Psi)
    })
    Score2_byChain[[j]] = map(1:n_samp, function(i){
      Xi = chain_samples[i,j,Score2_idx]
      dim(Xi) = c(N, K2)
      return(Xi)
    })
  }
  return(list(P1 = Psi1_byChain, S1 = Score1_byChain, 
              P2 = Psi2_byChain, S2 = Score2_byChain))
}
