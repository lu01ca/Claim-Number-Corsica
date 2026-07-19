
######################################
#                                    #
#     Negative Binomial finale       #
#                                    #
######################################

matrice_disegno_finale <- model.matrix(~ fascia_eta +
                                         log_densita +
                                         potenza_auto, data = df_stand)
colnames(matrice_disegno_finale)

stan_data_finale <- list(
  n = nrow(df_stand),
  p = ncol(matrice_disegno_finale),
  y = df_stand$num_sinistri,
  X = matrice_disegno_finale,
  exposure = df_stand$esposizione,
  beta0 = rep(0, ncol(matrice_disegno_finale)),
  Sigma0 = diag(1, ncol(matrice_disegno_finale))
)

fit_finale <- stan(file = "nb_corsica.stan",
                   data = stan_data_finale,
                   chains = 4,
                   warmup = 1000,
                   iter = 3000,
                   seed = 1234,
                   control = list(adapt_delta = 0.95))

# esplicitare quale parametri visualizzare
rstan::traceplot(fit_finale, pars = "beta", inc_warmup = TRUE)
rstan::traceplot(fit_finale, pars = "phi", inc_warmup = TRUE)

print(fit_finale, pars = c("beta", "phi"))
# Risultato catene per fare diagnostica
param_finale <- As.mcmc.list(fit_finale, pars = c("beta", "phi"))
summary(param_finale)
# controllo se ce autocorrelazione e se la catena è stabile
geweke.diag(param_finale)
# tolto il warmup controlla se la catena è stazionaria
# calcola i zeta score, maggiore di |2| indica non stazionarietà


# Intervallo di credibilità
ci(param_finale, method = "ETI")
ci(param_finale, method = "HDI")
mcmc_intervals(fit_finale, pars = c("beta[1]", "beta[2]",
                                    "beta[3]", "beta[4]",
                                    "beta[5]", "phi"))


# Valutazione del modello
log_lik_finale <- extract(fit_finale, pars = c("log_lik"))$log_lik
out_dens_finale <- exp(log_lik_finale)

# Calcolo WAIC, LPML e LOO
WAIC_finale <- my_waic(out_dens_finale)
LPML_finale <- my_lpml(out_dens_finale)
LOO_finale <- loo(log_lik_finale)$estimates["looic", "Estimate"]

# Calcolo del Bayes Factor (Necessario girare file STAN in questa sessione)
bs_finale <- bridge_sampler(fit_finale, silent = TRUE)
Log_Marginale_finale <- bs_finale$logml

# Salvataggio risultati
metriche_finale <- list(
  Modello = "Modello_Finale",
  WAIC = WAIC_finale,
  LPML = LPML_finale,
  LOO_IC = LOO_finale,
  Log_Marginale = Log_Marginale_finale
)

print(metriche_finale)

# Confronto modello NB completo vs modello NB ridotto
metriche_nb$WAIC; metriche_finale$WAIC  # scelgo il più piccolo
metriche_nb$LPML; metriche_finale$LPML  # scelgo il più grande
metriche_nb$LOO_IC; metriche_finale$LOO_IC # scelgo il più piccolo
exp(metriche_finale$Log_Marginale - metriche_nb$Log_Marginale) #bf a favore del modello finale

# confermato il modello ridotto per parsimonia, il modello finale esclude le variabili irrilevanti del NB completo

################################
#                              #
#     Diagnostica finale (NB)  #
#                              #
################################
post_finale <- extract(fit_finale)
beta_finale <- post_finale$beta; phi_finale <- post_finale$phi
n_sim <- 1000; idx <- sample(nrow(beta_finale), n_sim); N <- nrow(df_stand)

# Genera replicazioni dal modello FINALE
y_rep_finale <- matrix(NA, n_sim, N)
for (s in seq_len(n_sim)) {
  mu <- df_stand$esposizione * exp(matrice_disegno_finale %*% beta_finale[idx[s], ])
  y_rep_finale[s, ] <- rnbinom(N, mu = mu, size = phi_finale[idx[s]])
}

# Check dispersione
disp_obs_finale <- var(df_stand$num_sinistri) / mean(df_stand$num_sinistri)
disp_rep_finale <- apply(y_rep_finale, 1, function(x) var(x) / mean(x))
p_disp_finale <- mean(disp_rep_finale >= disp_obs_finale)
cat("Dispersione osservata:", round(disp_obs_finale, 3), "\n")
cat("p-value dispersione (NB):", round(p_disp_finale, 3), "\n")
# Sovradispersione catturata

# Check eccesso di zeri
prop_zero_obs_finale <- mean(df_stand$num_sinistri == 0)
prop_zero_rep_finale <- apply(y_rep_finale, 1, function(x) mean(x == 0))
p_zero_finale <- mean(prop_zero_rep_finale >= prop_zero_obs_finale)
cat("Prop. zeri osservata:", round(prop_zero_obs_finale, 4), "\n")
cat("Prop. zeri media simulata:", round(mean(prop_zero_rep_finale), 4), "\n")
cat("p-value zeri (NB):", round(p_zero_finale, 3), "\n")
# niente problema di eccesso di zeri

# Check simulazione del campione
tab_finale <- rbind(
  reale  = as.integer(table(factor(df_stand$num_sinistri, levels=0:3))),
  sim_NB = round(rowMeans(sapply(1:n_sim, function(s)
    table(factor(y_rep_finale[s,], levels=0:3)))), 1)   # y_rep_finale
)
colnames(tab_finale) <- 0:3; print(tab_finale)
# il modello riesce a simulare adeguatamente il campione


################################
#                              #
#   Salvataggio risultati      #
#                              #
################################


metriche <- bind_rows(
  metriche_base,
  metriche_m2,
  metriche_m3,
  metriche_nb,
  metriche_glmm,
  metriche_finale
)
print(metriche)
