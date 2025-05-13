data {
  int N;  // Number of observations
  int P;  // Number of fixed effect covariates
  
  array[N] int<lower=0, upper=1> Y;  // Binary outcomes
  matrix[N, P] X;                    // Fixed effects design matrix
}
transformed data {
  matrix[N, P] Q_coef = qr_thin_Q(X) * sqrt(N-1);
  matrix[P, P] R_coef = qr_thin_R(X) / sqrt(N-1);
  matrix[P, P] R_coef_inv = inverse(R_coef);
}
parameters {
  vector[P] theta;  // Coefficients
}
model {
  profile("Priors") {
    target += normal_lpdf(theta| 0, 100);
  }
  profile("Likelihood") {
    target += bernoulli_logit_lpmf(Y| Q_coef * theta);
  }
}
generated quantities {
  profile("Generated") {
    vector[P] beta = R_coef_inv * theta;
  }
}