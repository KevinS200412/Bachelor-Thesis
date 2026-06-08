// Bayesian Sigmoidal Quantile Regression (BSQR)
// Implements the GSCQF likelihood of Hutson (2024, 2026)
// embedded in a Bayesian framework following Yu & Moyeed (2001).
//
// Prior types:
//   1 = Normal       (ridge / L2)
//   2 = Student-t(3) (heavy-tailed)
//   3 = Cauchy       (heavy-tailed)
//   4 = Laplace      (Bayesian lasso / L1)

data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;

  real<lower=0, upper=1> q;
  real<lower=0> bw;

  int<lower=1, upper=4> prior_type;
  vector[K] beta_loc;
  vector<lower=0>[K] beta_scale;
}

parameters {
  vector[K] beta;
}

model {
  // --- Prior ---
  if (prior_type == 1) {
    beta ~ normal(beta_loc, beta_scale);
  } else if (prior_type == 2) {
    beta ~ student_t(3, beta_loc, beta_scale);
  } else if (prior_type == 3) {
    beta ~ cauchy(beta_loc, beta_scale);
  } else {
    // Laplace (double-exponential) — Bayesian lasso L1
    beta ~ double_exponential(beta_loc, beta_scale);
  }

  // --- Sigmoidal log-likelihood (GSCQF) ---
  vector[N] r = y - X * beta;
  target += -sum(bw * (log1p_exp(-r / bw) + log1p_exp(r / bw))
                 + (2.0 * q - 1.0) * r);
}
