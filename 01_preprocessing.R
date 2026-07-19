library(dplyr)
library(rstan)
library(bayesplot)
library(ggplot2)
library(loo)
library(bridgesampling)
library(CASdatasets)
library(coda)
library(bayestestR)


####################################
#                                  #
#      FUNZIONI PER WAIC E LPML    #
#                                  #
####################################

my_waic <- function(out_pred_dens){
  LPPD <- sum(log(colMeans(out_pred_dens))) 
  p_waic <- sum(apply(out_pred_dens, 2, function(x) var(log(x))))
  return(-2 * LPPD + 2 * p_waic)
}

my_lpml <- function(out_pred_dens){
  out <- sum(log(1 / colMeans(1 / out_pred_dens)))
  return(out)
}

# Caricamento dei modelli salvati
# fit_base <- readRDS("modelli_salvati/fit_base.rds")
# fit_modello2 <- readRDS("modelli_salvati/fit_modello2.rds")
# fit_modello3 <- readRDS("modelli_salvati/fit_modello3.rds")
# fit_nb <- readRDS("modelli_salvati/fit_nb.rds")
# fit_glmm <- readRDS("modelli_salvati/fit_glmm.rds")
# fit_nb_ssvs <- readRDS("modelli_salvati/fit_nb_ssvs.rds")
# fit_finale <- readRDS("modelli_salvati/fit_finale.rds")
# metriche <- readRDS("modelli_salvati/metriche.rds")


data("freMTPL2freq")
df <- freMTPL2freq
head(df)

#############################
#                           #
#     Pre-processing        #
#                           #
#############################

names(df)
names(df) = c('ID_pol','num_sinistri','esposizione',
              'potenza_auto','eta_auto', 'eta_conducente',
              'bonus_malus', 'marca_veicolo', 'tipo_alimentazione',
              'tipo_zona', 'densita_popolazione', 'regione')

names(df)

df <- df[, c(2:7, 11, 8:10, 12)]

# escludo la marca veicoli e considero solo la regione della corsica
df_corsica <- df[which(df$regione == "Corse"), -c(8, 11)]


##################
#                #
#     EDA        #
#                #
##################

# statistiche descrittive
summary(df_corsica)
table(df_corsica$num_sinistri)

df_corsica$tipo_zona <- factor(df_corsica$tipo_zona, levels = c("A", "B", "C", "D"))

# trasformata logaritmica alle densità
hist(df_corsica$densita_popolazione, freq = F)
df_corsica$log_densita = log(df_corsica$densita_popolazione)
hist(df_corsica$log_densita, freq = F)

# media incidenti per tipologia di carburante
aggregate(num_sinistri ~ tipo_alimentazione, data = df_corsica, FUN = mean)

# media incidenti per area geografica
aggregate(num_sinistri ~ tipo_zona, data = df_corsica, FUN = mean)

# Standardizzare i predittori
df_stand <- scale(df_corsica[, -c(1, 2, 7, 8, 9)])
df_stand <- as.data.frame(df_stand)
df_stand$tipo_alimentazione <- df_corsica$tipo_alimentazione
df_stand$tipo_alimentazione <- relevel(factor(df_stand$tipo_alimentazione), ref = "Diesel")
df_stand$tipo_zona <- df_corsica$tipo_zona
df_stand$esposizione <- df_corsica$esposizione
df_stand$num_sinistri <- df_corsica$num_sinistri

# creazione delle fasce d'età per conducenti e veicoli da usare come confronto
df_corsica <- df_corsica %>%
  mutate(fascia_eta = case_when(
    eta_conducente <= 25 ~ "Giovani",
    eta_conducente > 25 & eta_conducente <= 65 ~ "Adulti",
    eta_conducente > 65 ~ "Anziani"
  ))
df_stand$fascia_eta <- factor(df_corsica$fascia_eta, levels = c("Adulti", "Giovani", "Anziani"))

df_corsica <- df_corsica %>%
  mutate(fascia_eta_auto = case_when(
    eta_auto <= 1 ~ "meno di 1 anno",
    eta_auto > 1 & eta_auto <= 3 ~ "1-3 anni",
    eta_auto > 3 ~ "più di 3 anni"
  ))
df_stand$fascia_eta_auto <- factor(df_corsica$fascia_eta_auto,
                                   levels = c("meno di 1 anno", "1-3 anni", "più di 3 anni"))

df_stand <- df_stand[, c(1:5, 8, 6, 7, 10, 11, 9)]

# Analisi della correlazione
R <- cor(df_stand[, c(1:6, 11)])
ggcorrplot::ggcorrplot(R)
# moderata correlaazione negativa tra età conducente e bonus/malus, il resto ok
