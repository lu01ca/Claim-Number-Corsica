

####################################
#                                  #
#     Negative Binomial SSVS       #
#                                  #
####################################

# Predittori SEMPRE inclusi: intercetta + fascia_eta (categoria base = Adulti)
X_fix <- model.matrix(~ fascia_eta, data = df_stand)

# Predittori SOGGETTI a selezione (tutti standardizzati/dummy con base)
X_sel <- model.matrix(~ eta_auto + tipo_alimentazione +
                        log_densita + potenza_auto +
                        bonus_malus, data = df_stand)[, -1]
colnames(X_sel)

stan_data_ssvs <- list(
  N = nrow(df_stand),
  p_fix = ncol(X_fix),
  p_sel = ncol(X_sel),
  y = df_stand$num_sinistri,
  X_fix = X_fix,
  X_sel = X_sel,
  exposure = df_stand$esposizione,
  ngr = max(df_stand$zona_id),
  group = df_stand$zona_id,
  c = 0.001, tau2 = 10^3, theta = 0.5
)

fit_nb_ssvs <- stan(file = "nb_ssvs_corsica.stan",
                    data = stan_data_ssvs,
                    chains = 4,
                    warmup = 1000,
                    iter = 4000,
                    seed = 1234,
                    control = list(adapt_delta = 0.99))

print(fit_nb_ssvs, pars = c("beta_fix","beta_sel","gamma","sigma","phi"))
param_NB_ssvs <- As.mcmc.list(fit_nb_ssvs, pars = c("beta_fix","beta_sel","gamma","sigma","phi"))
summary(param_NB_ssvs)
geweke.diag(param_NB_ssvs)

ci(param_NB_ssvs, method = "ETI")
ci(param_NB_ssvs, method = "HDI")

# selezione delle variabili
eps <- sqrt(2 * (log(0.001) * 0.001^2) / (0.001^2 - 1))
(kappa <- sqrt(10^3) * eps)

beta_sel_dist <- extract(fit_nb_ssvs, pars = "beta_sel")$beta_sel
gamma_mat <- ifelse(abs(beta_sel_dist) > kappa, 1, 0)
pip <- colMeans(gamma_mat)
names(pip) <- colnames(X_sel)
print(round(pip, 4))

df_pip_plot <- data.frame(
  value = pip,
  var = factor(names(pip), levels = rev(names(pip)))
)

ggplot(df_pip_plot) + 
  geom_bar(aes(y = value, x = var), stat="identity", fill = "darkred", alpha = 0.7, col = 1) + 
  geom_hline(yintercept = 0.5, col = "black", lty = 2, linewidth = 1) +
  coord_flip() + 
  theme_minimal() + 
  labs(title = "Posterior Inclusion Probabilities (PIP) via SSVS",
       y = "PIP", x = "")

# HDP estimate
# seleziono i modelli unici
(unique_model <- unique(gamma_mat, MARGIN  = 1))
# conteggio le frequenze dei modelli unici
freq <- apply(unique_model, 1, function(b) sum(apply(gamma_mat, MARGIN = 1, function(a) all(a == b))))
# il modello finale è il modello più frequente
(HPD_model <- unique_model[which.max(freq),])

# MPM estimate
(MPM <- as.numeric(pip > 0.5)); names(MPM) <- colnames(X_sel); print(MPM)

# HS estimate
(HS_model <- as.numeric(colMeans(gamma_mat) == 1))

# ssvs seleziona le variabili potenza (buona rilevanza), log_densita (moderata influenza)
# e alimentazione (moderata rilevanza), oltre alla fascia d'età

# Per non aggiungere al modello un predittore in più poco rievante, rimuovo soltanto l'alimentazione
