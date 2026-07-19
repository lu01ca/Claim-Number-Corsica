data {
  int<lower=0> N;
  int<lower=0> p_fix;
  int<lower=0> y[N];
  matrix[N, p_fix] X;
  vector<lower=0>[N] exposure;
  int<lower=1> ngr;
  int<lower=1, upper=ngr> group[N];
  real<lower=0> sigma20;
}
parameters {
  vector[p_fix] beta;
  vector[ngr] a_raw;
  real<lower=0> sigma;
  real<lower=0> phi;          // <-- dispersione NB
}
transformed parameters {
  vector[ngr] gamma;
  {
    vector[ngr] a_tmp = a_raw * sigma;
    gamma = a_tmp - mean(a_tmp);   // sum-to-zero
  }
  vector[N] log_mu;
  for (i in 1:N)
    log_mu[i] = log(exposure[i]) + X[i] * beta + gamma[group[i]];
}
model {
  beta  ~ normal(0, pow(sigma20, 0.5));
  sigma ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  phi   ~ exponential(1);
  y ~ neg_binomial_2_log(log_mu, phi);   // <-- NB, non Poisson
}
generated quantities {
  vector[N] log_lik;
  for (j in 1:N)
    log_lik[j] = neg_binomial_2_log_lpmf(y[j] | log_mu[j], phi);
}
