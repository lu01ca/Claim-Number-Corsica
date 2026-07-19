
################
#              #
#     GLM      #
#              #
################
# uso tutti i predittori classici (no fasce d'età conducente e auto)
matrice_disegno = model.matrix(~ eta_conducente +
                                 eta_auto +
                                 tipo_alimentazione +
                                 log_densita +
                                 potenza_auto + 
                                 bonus_malus, data = df_stand)
colnames(matrice_disegno)

stan_data <- list(
  n = nrow(df_stand),                     
  p = ncol(matrice_disegno),                  
  y = df_stand$num_sinistri,              
  X = matrice_disegno,                        
  exposure = df_stand$esposizione,        
  
  # medie 0 e varianze 1 per prior non informative
  beta0 = rep(0, ncol(matrice_disegno)),      
  Sigma0 = diag(1, ncol(matrice_disegno))    
)

#primo modello
fit_base <- stan(
  file = "poisson_corsica.stan", 
  data = stan_data,                
  chains = 4,     # più catene più ergodicità                 
  iter = 2000,                     
  warmup = 1000,                   
  cores = 4,
  seed = 104
)

# esplicitare quale parametri visualizzare
rstan::traceplot(fit_base, pars = "beta", inc_warmup = TRUE)

# Risultato catene per fare diagnostica
param_GLM <- As.mcmc.list(fit_base, pars = c("beta"))
summary(param_GLM)
# controllo se ce autocorrelazione e se la catena è stabile
geweke.diag(param_GLM)
# tolto il warmup controlla se la catena è stazionaria
# calcola i zeta score, maggiore di |2| indica non stazionarietà


# Intervallo di credibilità
ci(param_GLM, method = "ETI")
ci(param_GLM, method = "HDI")
mcmc_intervals(fit_base, pars = c("beta[1]", "beta[2]",
                                  "beta[3]", "beta[4]",
                                  "beta[5]", "beta[6]",
                                  "beta[7]"))

# Valutazione del modello
log_lik_base <- extract(fit_base, pars = c("log_lik"))$log_lik
out_dens_base <- exp(log_lik_base)

# Calcolo WAIC, LPML e LOO
WAIC_base <- my_waic(out_dens_base)
LPML_base <- my_lpml(out_dens_base)
LOO_base <- loo(log_lik_base)$estimates["looic", "Estimate"]

# Calcolo del Bayes Factor (Necessario girare file STAN in questa sessione)
bs_base <- bridge_sampler(fit_base, silent = TRUE)
Log_Marginale_base <- bs_base$logml

# Salvataggio risultati
metriche_base <- list(
  Modello = "Poisson_base",
  WAIC = WAIC_base,
  LPML = LPML_base,
  LOO_IC = LOO_base,
  Log_Marginale = Log_Marginale_base
)

print(metriche_base)



################
#              #
#     GLM2     #
#              #
################

#uso fasce età del conducende e auto

matrice_disegno2 <- model.matrix(~ fascia_eta +
                                   fascia_eta_auto +
                                   tipo_alimentazione +
                                   log_densita +
                                   potenza_auto + 
                                   bonus_malus, data = df_stand)
colnames(matrice_disegno2)

stan_data_2 <- list(
  n = nrow(df_stand),
  p = ncol(matrice_disegno2),
  y = df_stand$num_sinistri,
  X = matrice_disegno2,
  exposure = df_stand$esposizione,
  
  # Prior quasi non informative adattate al nuovo numero di colonne
  beta0 = rep(0, ncol(matrice_disegno2)),
  Sigma0 = diag(1, ncol(matrice_disegno2))
)

#secondo modello
fit_modello2 <- stan(
  file = "poisson_corsica.stan", 
  data = stan_data_2,                
  chains = 4,                      
  iter = 2000,                     
  warmup = 1000,                   
  cores = 4,
  seed = 104
)
# esplicitare quale parametri visualizzare
rstan::traceplot(fit_modello2, pars = "beta", inc_warmup = TRUE)

# Risultato catene per fare diagnostica
param_GLM2 <- As.mcmc.list(fit_modello2, pars = c("beta"))
summary(param_GLM2)
# controllo se ce autocorrelazione e se la catena è stabile
geweke.diag(param_GLM2)
# tolto il warmup controlla se la catena è stazionaria
# calcola i zeta score, maggiore di |2| indica non stazionarietà


# Intervallo di credibilità
ci(param_GLM2, method = "ETI")
ci(param_GLM2, method = "HDI")

# Valutazione del modello
log_lik_m2 <- extract(fit_modello2, pars = c("log_lik"))$log_lik
out_dens_m2 <- exp(log_lik_m2)

# Calcolo WAIC, LPML e LOO
WAIC_m2 <- my_waic(out_dens_m2)
LPML_m2 <- my_lpml(out_dens_m2)
LOO_m2 <- loo(log_lik_m2)$estimates["looic", "Estimate"]

# Calcolo del Bayes Factor (Necessario girare file STAN in questa sessione)
bs_m2 <- bridge_sampler(fit_modello2, silent = TRUE)
Log_Marginale_m2 <- bs_m2$logml

# Salvataggio risultati
metriche_m2 <- list(
  Modello = "Poisson_Modello_2",
  WAIC = WAIC_m2,
  LPML = LPML_m2,
  LOO_IC = LOO_m2,
  Log_Marginale = Log_Marginale_m2
)

print(metriche_m2)

# Confronto modello base vs modello 2
metriche_base$WAIC; metriche_m2$WAIC  # scelgo il più piccolo
metriche_base$LPML; metriche_m2$LPML  # scelgo il più grande
metriche_base$LOO_IC; metriche_m2$LOO_IC # scelgo il più piccolo
exp(metriche_m2$Log_Marginale - metriche_base$Log_Marginale) #bf a favore del modello 2

# verdetto unanime, Modello 2 migliore ma va capito se categorizzare auto aiuta



################
#              #
#     GLM3     #
#              #
################

#uso fasce età del conducende e auto

matrice_disegno3 <- model.matrix(~ fascia_eta +
                                   eta_auto +
                                   tipo_alimentazione +
                                   log_densita +
                                   potenza_auto + 
                                   bonus_malus, data = df_stand)
colnames(matrice_disegno3)

stan_data_3 <- list(
  n = nrow(df_stand),
  p = ncol(matrice_disegno3),
  y = df_stand$num_sinistri,
  X = matrice_disegno3,
  exposure = df_stand$esposizione,
  
  # Prior non informative adattate al nuovo numero di colonne
  beta0 = rep(0, ncol(matrice_disegno3)),
  Sigma0 = diag(1, ncol(matrice_disegno3))
)

#terzo modello
fit_modello3 <- stan(
  file = "poisson_corsica.stan", 
  data = stan_data_3,                
  chains = 4,                      
  iter = 2000,                     
  warmup = 1000,                   
  cores = 4,
  seed = 104
)
# esplicitare quale parametri visualizzare
rstan::traceplot(fit_modello3, pars = "beta", inc_warmup = TRUE)

# Risultato catene per fare diagnostica
param_GLM3 <- As.mcmc.list(fit_modello3, pars = c("beta"))
summary(param_GLM3)
# controllo se ce autocorrelazione e se la catena è stabile
geweke.diag(param_GLM3)
# tolto il warmup controlla se la catena è stazionaria
# calcola i zeta score, maggiore di |2| indica non stazionarietà


# Intervallo di credibilità
ci(param_GLM3, method = "ETI")
ci(param_GLM3, method = "HDI")

# Valutazione del modello
log_lik_m3 <- extract(fit_modello3, pars = c("log_lik"))$log_lik
out_dens_m3 <- exp(log_lik_m3)

# Calcolo WAIC, LPML e LOO
WAIC_m3 <- my_waic(out_dens_m3)
LPML_m3 <- my_lpml(out_dens_m3)
LOO_m3 <- loo(log_lik_m3)$estimates["looic", "Estimate"]

# Calcolo del Bayes Factor (Necessario girare file STAN in questa sessione)
bs_m3 <- bridge_sampler(fit_modello3, silent = TRUE)
Log_Marginale_m3 <- bs_m3$logml

# Salvataggio risultati
metriche_m3 <- list(
  Modello = "Poisson_Modello_3",
  WAIC = WAIC_m3,
  LPML = LPML_m3,
  LOO_IC = LOO_m3,
  Log_Marginale = Log_Marginale_m3
)

print(metriche_m3)

# Confronto modello 2 vs modello 3
metriche_m2$WAIC; metriche_m3$WAIC  # scelgo il più piccolo
metriche_m2$LPML; metriche_m3$LPML  # scelgo il più grande
metriche_m2$LOO_IC; metriche_m3$LOO_IC # scelgo il più piccolo
exp(metriche_m2$Log_Marginale - metriche_m3$Log_Marginale) #bf a favore del modello 2

# Modello 2 e 3 simili, la scelta di m2 non giutifica la complessità di avere una variabile in più
# tengo età auto continua e età conducente categorica, quindi modello 3 è il modello di riferimento

#####################################
#                                   #
#     Analisi sovradispersione      #
#                                   #
#####################################

# Generazione posterior predictive modello 3
post_m3 <- extract(fit_modello3)
beta_m3 <- post_m3$beta
n_sim <- 1000
N <- nrow(df_stand)
y_rep <- matrix(NA, nrow = n_sim, ncol = N)
idx <- sample(nrow(beta_m3), n_sim)

for (s in seq_len(n_sim)) {
  lambda_stimato <- df_stand$esposizione * exp(matrice_disegno3 %*% beta_m3[idx[s], ])
  y_rep[s, ] <- rpois(N, lambda = lambda_stimato)
}

# Check dispersione
(disp_obs <- var(df_stand$num_sinistri) / mean(df_stand$num_sinistri))
disp_rep <- apply(y_rep, 1, function(x) var(x) / mean(x))

# Con quale frequenza il modello genera dati dispersi quanto o più degli osservati?
p_disp <- mean(disp_rep >= disp_obs)

cat("Dispersione osservata:", round(disp_obs, 3), "\n")
cat("p-value bayesiano (dispersione):", round(p_disp, 3), "\n")

# ce sovradispersione, modello di Poisson non va bene

# Check eccesso di zeri
prop_zero_obs <- mean(df_stand$num_sinistri == 0)
prop_zero_rep <- apply(y_rep, 1, function(x) mean(x == 0))
p_zero <- mean(prop_zero_rep >= prop_zero_obs)

# Il modello riesce a simulare i zeri del campione osservato?
cat("Prop. zeri osservata:", round(prop_zero_obs, 4), "\n")
cat("Prop. zeri media simulata:", round(mean(prop_zero_rep), 4), "\n")
cat("p-value bayesiano (zeri):", round(p_zero, 3), "\n")

# niente eccesso di zeri

# check poisson simula bene il campione?
tab <- rbind(
  reale = as.integer(table(factor(df_stand$num_sinistri, levels=0:3))),
  sim_Poisson = round(rowMeans(sapply(1:n_sim, function(s)
    table(factor(y_rep[s,], levels=0:3)))), 1)
)
colnames(tab) <- 0:3; print(tab)

# la poisson non riesce a catturare la coda a causa della sovradispersione
