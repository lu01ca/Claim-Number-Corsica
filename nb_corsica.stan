data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> y[n];
  matrix[n,p] X;
  vector<lower=0>[n] exposure;
  vector[p] beta0;
  matrix[p,p] Sigma0;
}
parameters {
  vector[p] beta;
  real<lower=0> phi;          // dispersione: var = mu + mu^2/phi
}
transformed parameters {
  vector[n] log_mu;
  for (i in 1:n)
    log_mu[i] = log(exposure[i]) + X[i] * beta;
}
model {
  beta ~ multi_normal(beta0, Sigma0);
  phi  ~ exponential(1);
  y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  vector[n] log_lik;
  for (j in 1:n)
    log_lik[j] = neg_binomial_2_log_lpmf(y[j] | log_mu[j], phi);
}

