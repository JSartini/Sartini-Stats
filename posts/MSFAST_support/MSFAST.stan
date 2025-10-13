data {
  int N;   // Number of time series/individuals
  int M;   // Cardinality of observation time set
  int L;   // Total number of observed data points
  int Q;   // Number of spline bases
  int K;   // Number of eigenfunctions
  int P;   // Number of covariates
  
  vector[L] Y;                         // Observed data for all variables stacked
  array[L] int<upper=N> Subj;          // Time series indices all variables stacked
  array[L] int<upper=M> S;             // Time point indices all variables stacked
  array[P] int Tp_card;                // Number of observed data points for each variable
  
  matrix[M, Q] B;               // Orthogonalized basis over T
  matrix[Q, Q] P_alpha;         // Penalty matrix for splines
}

parameters {
  vector<lower=0>[P] sigma2; // Error in observation
  
  // Fixed-effect components
  vector[P*Q] w_mu;         // Population mean parameters
  vector<lower=0>[P] h_mu;  // Population mean smoothing parameter
  
  // Components/weights
  positive_ordered[K] lambda;        // Eigenvalues
  matrix<lower=0>[P,K] H;            // EF Smoothing parameters
  matrix[P*Q, K] X;                  // Unconstrained EF matrix
  matrix[N, K] Xi_Raw;               // EF scores unscaled
}

transformed parameters{
  // Orthogonal basis weights
  matrix[P*Q, K] Psi;
  matrix[N, K] Scores;
  
  // Polar decomposition
  {
    matrix[K,K] evec_XtX = eigenvectors_sym(crossprod(X)); 
    vector[K] eval_XtX = eigenvalues_sym(crossprod(X));
    Psi = X*evec_XtX*diag_matrix(1/sqrt(eval_XtX))*evec_XtX'; 
  }
  
  // Scaled scores
  Scores = Xi_Raw * diag_matrix(sqrt(lambda));
}

model {
  // Variance component priors
  lambda ~ inv_gamma(0.01, 0.01); 
  sigma2 ~ inv_gamma(0.01, 0.01);
  
  // Smoothing priors
  h_mu ~ gamma(0.01, 0.01); 
  to_vector(H) ~ gamma(0.01, 0.01); 
  
  int sdx;
  int edx;
  for(p in 1:P){
    sdx = (p-1)*Q+1;
    edx = p*Q;
    
    target += Q / 2.0 * log(h_mu[p]) - h_mu[p] / 2.0 * quad_form(P_alpha, w_mu[sdx:edx]);
    
    for(k in 1:K){
      target += Q / 2.0 * log(H[p,k]) -  H[p,k] / 2.0 * quad_form(P_alpha, Psi[sdx:edx,k]);
    }
  }
  
  // Uniform priors through matrix normals
  to_vector(X) ~ std_normal();
  
  // Score priors 
  to_vector(Xi_Raw) ~ std_normal();
  
  // Model likelihood
  int pos = 1;
  for(p in 1:P){
    array[Tp_card[p]] int Subj_p;
    array[Tp_card[p]] int S_p;
    matrix[M, K] Phi_mat;
    vector[Tp_card[p]] Theta;
    vector[M] mu;
    
    // Constants/subsets
    sdx = (p-1)*Q+1;
    edx = p*Q;
    Subj_p = segment(Subj, pos, Tp_card[p]);
    S_p = segment(S, pos, Tp_card[p]);
    
    // Smooth estimates for component
    Phi_mat = B * Psi[sdx:edx,];   // Eigenfunctions
    Theta = rows_dot_product(Scores[Subj_p, ], Phi_mat[S_p, ]); // Smooth deviations
    mu = B * w_mu[sdx:edx];             // Population mean
    
    // Likelihood for this variable
    segment(Y, pos, Tp_card[p]) ~ normal(mu[S_p] + Theta, sqrt(sigma2[p]));
    
    // Increment position
    pos = pos + Tp_card[p];
  }
}
