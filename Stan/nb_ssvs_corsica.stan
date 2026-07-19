data {
  int<lower=0> N;
  int<lower=0> p_fix;                 // SEMPRE inclusi (intercetta + fascia_eta)
  int<lower=0> p_sel;                 // soggetti a selezione
  int<lower=0> y[N];
  matrix[N, p_fix] X_fix;
  matrix[N, p_sel] X_sel;
  vector<lower=0>[N] exposure;
  int<lower=1> ngr;
  int<lower=1, upper=ngr> group[N];
  real<lower=0> c;
  real<lower=0> tau2;
  real<lower=0, upper=1> theta;
}
parameters {
  vector[p_fix] beta_fix;
  vector[p_sel] beta_sel;
  vector[ngr] a_raw;
  real<lower=0> sigma;
  real<lower=0> phi;                  // <-- dispersione NB
}
transformed parameters {
  vector[ngr] gamma;
  {
    vector[ngr] a_tmp = a_raw * sigma;
    gamma = a_tmp - mean(a_tmp);
  }
}
model {
  beta_fix ~ normal(0, 5);

  for (j in 1:p_sel) {
    real log_spike = log(1 - theta) + normal_lpdf(beta_sel[j] | 0, sqrt(c * tau2));
    real log_slab  = log(theta)     + normal_lpdf(beta_sel[j] | 0, sqrt(tau2));
    target += log_sum_exp(log_spike, log_slab);
  }

  sigma ~ normal(0, 1);
  a_raw   ~ normal(0, 1);
  phi     ~ exponential(1);           // <-- prior dispersione

  for (i in 1:N) {
    real log_lambda = log(exposure[i]) + X_fix[i] * beta_fix
                      + X_sel[i] * beta_sel + gamma[group[i]];
    y[i] ~ neg_binomial_2_log(log_lambda, phi);   // <-- NB
  }
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real log_lambda = log(exposure[i]) + X_fix[i] * beta_fix
                      + X_sel[i] * beta_sel + gamma[group[i]];
    log_lik[i] = neg_binomial_2_log_lpmf(y[i] | log_lambda, phi);
  }
}


