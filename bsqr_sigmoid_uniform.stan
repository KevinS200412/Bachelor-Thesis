// Bayesian Sigmoidal Quantile Regression — Uniform Prior
// Parameters are bounded to [-bound, bound] per the Uniform(-c, c) prior.
// The bound c is passed as data so one model handles all three Uniform specs.

data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;

  real<lower=0, upper=1> q;
  real<lower=0> bw;
  real<lower=0> bound;              // c in Uniform(-c, c)
}

parameters {
  vector<lower=-bound, upper=bound>[K] beta;
}

model {
  // Uniform prior is implicit from the bounded parameter declaration.
  // No explicit prior statement needed — Stan uses uniform over the support.

  // --- Sigmoidal log-likelihood (GSCQF) ---
  vector[N] r = y - X * beta;
  target += -sum(bw * (log1p_exp(-r / bw) + log1p_exp(r / bw))
                 + (2.0 * q - 1.0) * r);
}
