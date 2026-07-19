

##############################
#                            #
#     Negative Binomial      #
#                            #
##############################

# matrice disegno modello 3

stan_data_nb <- list(
  n = nrow(df_stand),
  p = ncol(matrice_disegno3),
  y = df_stand$num_sinistri,
  X = matrice_disegno3,
  exposure = df_stand$esposizione,
  beta0 = rep(0, ncol(matrice_disegno3)),
  Sigma0 = diag(1, ncol(matrice_disegno3))
)

fit_nb <- stan(file = "nb_corsica.stan",
               data = stan_data_nb,
               chains = 4,
               warmup = 1000,
               iter = 2000,
               seed = 1234,
               control = list(adapt_delta = 0.95))

# esplicitare quale parametri visualizzare
rstan::traceplot(fit_nb, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(fit_nb, pars = "phi", inc_warmup = TRUE)

# Risultato catene per fare diagnostica
param_nb <- As.mcmc.list(fit_nb, pars = c("beta", "phi"))
summary(param_nb)
# controllo se ce autocorrelazione e se la catena è stabile
geweke.diag(param_nb)
# tolto il warmup controlla se la catena è stazionaria
# calcola i zeta score, maggiore di |2| indica non stazionarietà


# Intervallo di credibilità
ci(param_nb, method = "ETI")
ci(param_nb, method = "HDI")

# Valutazione del modello
log_lik_nb <- extract(fit_nb, pars = c("log_lik"))$log_lik
out_dens_nb <- exp(log_lik_nb)

# Calcolo WAIC, LPML e LOO
WAIC_nb <- my_waic(out_dens_nb)
LPML_nb <- my_lpml(out_dens_nb)
LOO_nb <- loo(log_lik_nb)$estimates["looic", "Estimate"]

# Calcolo del Bayes Factor (Necessario girare file STAN in questa sessione)
bs_nb <- bridge_sampler(fit_nb, silent = TRUE)
Log_Marginale_nb <- bs_nb$logml

# Salvataggio risultati
metriche_nb <- list(
  Modello = "Negative_Binomial",
  WAIC = WAIC_nb,
  LPML = LPML_nb,
  LOO_IC = LOO_nb,
  Log_Marginale = Log_Marginale_nb
)

print(metriche_nb)

# Confronto modello NB vs modello 3
metriche_m3$WAIC; metriche_nb$WAIC  # scelgo il più piccolo
metriche_m3$LPML; metriche_nb$LPML  # scelgo il più grande
metriche_m3$LOO_IC; metriche_nb$LOO_IC # scelgo il più piccolo
exp(metriche_nb$Log_Marginale - metriche_m3$Log_Marginale) #bf a favore del modello nb

# il giudizio è unanime, il modello NB con i stessi predittori è migliore del modello poisson

# Ri-check dispersione sul NB
post_nb <- extract(fit_nb)
beta_nb <- post_nb$beta; phi_nb <- post_nb$phi
n_sim <- 1000; idx <- sample(nrow(beta_nb), n_sim); N <- nrow(df_stand)

y_rep_nb <- matrix(NA, n_sim, N)
for (s in seq_len(n_sim)) {
  mu <- df_stand$esposizione * exp(matrice_disegno3 %*% beta_nb[idx[s], ])
  y_rep_nb[s, ] <- rnbinom(N, mu = mu, size = phi_nb[idx[s]])
}
disp_obs <- var(df_stand$num_sinistri) / mean(df_stand$num_sinistri)
disp_rep_nb <- apply(y_rep_nb, 1, function(x) var(x)/mean(x))
cat("p-value dispersione (NB):", round(mean(disp_rep_nb >= disp_obs), 3), "\n")

# il modello NB cattura correttamente la sovradispersione a discapito del Poisson


####################################
#                                  #
#     Negative Binomial Mixed      #
#                                  #
####################################

#rendo integer per fare funzionare stan
df_stand$zona_id <- as.integer(as.factor(df_stand$tipo_zona))
n_zone <- max(df_stand$zona_id)

# uso la matrice disegno del modello 3
colnames(matrice_disegno3)

stan_data_glmm <- list(
  N        = nrow(df_stand),
  p_fix    = ncol(matrice_disegno3),
  y        = df_stand$num_sinistri,
  X        = matrice_disegno3,
  exposure = df_stand$esposizione,
  ngr      = n_zone,
  group    = df_stand$zona_id,
  sigma20  = 1
)

fit_glmm <- stan(file = "nb_mixed_corsica.stan",
                 data = stan_data_glmm,
                 chains = 4,
                 warmup = 1000,
                 iter = 2000,
                 seed = 1234,
                 control = list(adapt_delta = 0.95))

print(fit_glmm, pars = c("beta", "gamma", "sigma"))

rstan::traceplot(fit_glmm, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(fit_glmm, pars= "gamma", inc_warmup = TRUE)
rstan::traceplot(fit_glmm, pars = "sigma", inc_warmup = TRUE)

param_GLMM <- As.mcmc.list(fit_glmm, pars = c("beta", "gamma", "sigma"))
summary(param_GLMM)
geweke.diag(param_GLMM)

ci(param_GLMM, method = "ETI")
ci(param_GLMM, method = "HDI")

# Valutazione del modello
log_lik_glmm <- extract(fit_glmm, pars = c("log_lik"))$log_lik
out_dens_glmm <- exp(log_lik_glmm)

# Calcolo WAIC, LPML e LOO
WAIC_glmm <- my_waic(out_dens_glmm)
LPML_glmm <- my_lpml(out_dens_glmm)
LOO_glmm <- loo(log_lik_glmm)$estimates["looic", "Estimate"]

# Calcolo del Bayes Factor (Necessario girare file STAN in questa sessione)
bs_glmm <- bridge_sampler(fit_glmm, silent = TRUE)
Log_Marginale_glmm <- bs_glmm$logml

# Salvataggio risultati
metriche_glmm <- list(
  Modello = "Negative_Binomial_Mixed",
  WAIC = WAIC_glmm,
  LPML = LPML_glmm,
  LOO_IC = LOO_glmm,
  Log_Marginale = Log_Marginale_glmm
)

print(metriche_glmm)

# Confronto modello NB vs modello 3
metriche_nb$WAIC; metriche_glmm$WAIC  # scelgo il più piccolo
metriche_nb$LPML; metriche_glmm$LPML  # scelgo il più grande
metriche_nb$LOO_IC; metriche_glmm$LOO_IC # scelgo il più piccolo
exp(metriche_glmm$Log_Marginale - metriche_nb$Log_Marginale) #bf a favore del modello glmm

# i due modelli sono simili, l'effetto territoriale non aggiunge nulla, meglio modello NB base per parsimonia
