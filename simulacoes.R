# install.packages("devtools")
# install.packages("tidyverse")
detach("package:mixRegEM", unload=TRUE)
devtools::install_github("https://github.com/PedroHCAlmeida/mixRegEM",
                         force = TRUE)

library(tidyverse)
library(mixsmsn)
library(mixRegEM)

matrizP <- function(alpha, R){
  P1 <- exp(R%*%alpha)/(1 + rowSums(exp(R%*%alpha)))
  P <- cbind(P1, 1 - rowSums(P1))
  return(P)
}

# Função auxiliar
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

precision = function(x){
  p = sum(diag(x))/sum(x)
  if(p<0.5){
    p = sum(diag(x[,c(2,1)]))/sum(x)
  }
  p
}

# Simulando dados

set.seed(123)
n = 3000
X = cbind(rep(1, n), rnorm(n, 5, 3), rnorm(n, 10, 5))
x = X[,-1]

R = cbind(rep(1, n), runif(n, -2, 1), runif(n, -1, 1))
r = R[,-1]

beta01 <- c(5, 6, -1)
alpha01 <- c(0.7, 1, 2)
sigma01 <- 100
mu01 <- X%*%beta01

beta02 <- c(-10, 1, 4)
sigma02 <- 30
mu02 <- X%*%beta02

alpha <- matrix(alpha01, byrow = T, nrow = 3)
P <- matrizP(alpha, R)
table(apply(P, 1, which.max))

obs = do.call(rbind, 
              lapply(1:n,
                     function(i){
                       arg1 <- list(mu = mu01[i,], sigma2 = sigma01, shape = 0)
                       arg2 <- list(mu = mu02[i,], sigma2 = sigma02, shape = 0)
                       
                       obs <- mixsmsn::rmix(1, pii = P[i,], family = 'Normal', 
                                   arg = list(arg1, arg2), cluster = T)
                       
                       return(obs)
                     }))

rm(beta01, alpha01, sigma01, mu01, beta02, sigma02, mu02, alpha)
y = unlist(obs[,1])
grupo = unlist(obs[,2]) 

table(grupo)
ng = 2

# Analisando Grupos

data.frame(grupo = as.factor(grupo), y = y) |>
  ggplot(aes(x = y, fill = grupo)) +
  geom_histogram(bins = 15, alpha=0.6, position = 'identity') +
  theme_minimal()

# Regressão usual

ml = lm(y ~ ., data.frame(cbind(x, y = y)))
summary(ml)

# Misturas

regMoEN = regEM(y, cbind(x, runif(n,0,0.5)), r = r, g = ng, tol = 1E-6,
                grupoReal = grupo,
                family = "MoENormal", initGrupo = "KMeans", min_iter = 1)

regMoEN$eps
regMoEN$tabela
precision(regMoEN$tabela)
regMoEN$l
regMoEN$metricas
regMoEN$Parametros

regMixN = regEM(y, x, g = ng, tol = 1E-6, grupoReal = grupo,
                family = "MixNormal")
regMixN$eps
regMixN$tabela
regMixN$metricas
precision(regMixN$tabela)
regMixN$Parametros

# Análise resíduos

dadosMetodos = data.frame(do.call(cbind.fill, regMoEN$residuos)) |>
  pivot_longer(everything()) |> 
  filter(!is.na(value)) |>
  mutate(Metodo = "MoE") |>
  bind_rows(
    data.frame(do.call(cbind.fill, regMixN$residuos)) |>
      pivot_longer(everything()) |> 
      filter(!is.na(value)) |>
      mutate(Metodo = "Mix")) |>
  rename(Grupo = name) |>
  mutate(Grupo = str_replace(Grupo, "X", ""))

dadosMetodos |>
  ggplot(aes(value, ..density.., fill = Grupo)) +
  facet_grid(vars(Grupo), vars(Metodo)) +
  geom_histogram() +
  theme_minimal()

dadosMetodos |>
  ggplot(aes(Grupo, value, fill = Grupo)) +
  facet_wrap(~Metodo, nrow = 2) +
  geom_boxplot() +
  theme_minimal()

#################################################

set.seed(seed)
train = sample(1:n, 
               floor(n*0.7), replace = F)
test = (1:n)[!1:n %in% train]

regMoENtrain = regEM(y[train], x[train,], 
                     r = r[train,], g = ng, tol = 1E-6,
                     family = "MoENormal", initGrupo = "KMeans")

predicao_teste_1 = predictMix(regMoENtrain, x[test,],
                              r[test,], type = 1)
plot(y[test], predicao_teste_1)

predicao_teste_2 = predictMix(regMoENtrain, x[test,],
                              r[test,], type = 2)
plot(y[test], predicao_teste_2)

erros_1 = predicao_teste_1 - y[test]
hist(erros_1)
mean(erros_1)
mean(abs(erros_1))
sqrt(mean(erros_1**2))

erros_2 = predicao_teste_2 - y[test]
hist(erros_2)
mean(erros_2)
mean(abs(erros_2))
sqrt(mean(erros_2**2))

ml = lm(y ~ ., data.frame(cbind(x[train,], 
                     y = y[train])))
erros_ml = y[test] - (ml$coefficients %in% X[test,])
hist(erros_ml)
mean(erros_ml)
mean(abs(erros_ml))
sqrt(mean(erros_ml**2))

#################################################

library(mixtools)

regmix = regmixEM(y, x)
reghhme = hmeEM(y, x)

################################################
