# Mistura dos experts ST com resposta censurada à esquerda
library(MomTrunc) # Sempre passar a sigma² como parâmetro nas funções
library(sn)
library(moments)
library(mixsmsn)
library(numDeriv)
#library(mixRegEM)

matrizP2 <- function(alpha, R){
  P1 <- exp(R%*%t(alpha))/(1 + rowSums(exp(R%*%t(alpha))))
  P <- cbind(P1, 1 - rowSums(P1))
  return(P)
}

n = c(1000)
nivelC = c(0)
g = 2
tol = 1E-5

beta01 <- c(0, -1, -2, -3)
alpha01 <- c(0.7, 1, 2)
sigma2_01 <- 1
lambda01 <- -1

beta02 <- c(-1, 1, 2, 3)
sigma2_02 <- 2
lambda02 <- 3

alpha <- matrix(alpha01, byrow = T, nrow = 1)

nu = 3

rMoeEM = function(ni, ci, tol = 1E-4, verbose = F){

  covergiu = F

    while(!covergiu){

      X <- cbind(rep(1, ni), runif(ni, 1, 5), runif(ni, -2, 2), runif(ni, 1, 4))
      R <- cbind(rep(1, ni), runif(ni, -2, 1), runif(ni, -1, 1))

      P <- matrizP2(alpha, R)
      mu01 <- as.numeric(X%*%beta01)
      mu02 <- as.numeric(X%*%beta02)

      grupo <- numeric(ni)
      y <- numeric(ni)
      for(ii in 1:ni){
        arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01, nu = nu)
        arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02, nu = nu)

        obs <- rmix(1, pii = P[ii,], family = 'Skew.t', arg = list(arg1, arg2), cluster = T)

        y[ii] <- obs$y
        grupo[ii] <- obs$cluster
      }

      ki <- as.numeric(quantile(y, probs = ci))
      y[y <= ki] <- ki
      phi <- as.numeric(y == ki)
      resultados = NULL
      try({
        resultados = regEM(
          y,
          X[,-1],
          r = R[,-1],
          family = "MoECenST",
          phi = phi,
          c1 = -Inf,
          c2 = ki,
          g = 2,
          showSE = F,
          verbose = F,
          tol = tol,
          max_iter = 10000,
          min_iter = 1
        )
        })
      if(is.null(resultados)) convergiu = F
      else convergiu = resultados$convergiu
      }
    return(resultados$Parametros)
    }
set.seed(123)
start = Sys.time()
resultadosMoECenSN = n |>
  sapply(
    function(ni){
      nivelC |>
        lapply(
          function(ci) replicate(100, rMoeEM(ni, ci, tol = 1E-5, verbose = F))
        ) |>
        setNames(paste("Cen =", nivelC))
    }, simplify = "array"
  )
end = Sys.time()
print(end-start)

parametros =cbind(c(
  beta01, sqrt(sigma2_01), lambda01, alpha01
),
c(
  beta02, sqrt(sigma2_02), lambda02, rep(NA, 3)
)
)

rownames(parametros) = c(paste0("beta", 1:4), "sigma", "lambda", paste0("alpha", 1:3))

inverte_col = function(df){
  df_inv = df
  df_inv[,1] = df[,2]
  df_inv[,2] = df[,1]

  df_inv[startsWith(rownames(df), "alpha"),1] = -df_inv[startsWith(rownames(df), "alpha"),2]
  df_inv[startsWith(rownames(df), "alpha"),2] = NA

  df_inv
}

erros = apply(
  resultadosMoECenSN[[1]],
  3,
  function(x){
    print(x)
    res1 = x[startsWith(rownames(x), "beta"),1] - parametros[startsWith(rownames(parametros), "beta"),1]
    res2 = x[startsWith(rownames(x), "beta"),1] - parametros[startsWith(rownames(parametros), "beta"),2]

    if(sum(res1**2) > sum(res2**2)) x = inverte_col(x)

    params = rownames(x)
    x = x[!params %in% c("delta", "gama"),]

    #x = x - parametros

    x = data.frame(x)
    x$parametros1 = parametros[,1]
    x$parametros2 = parametros[,2]
    x$n = n
    x$cen = nivelC
    x$params = params[!params %in% c("delta", "gama")]

    return(x)
  }
) |>
  dplyr::bind_rows() |>
  dplyr::as_tibble()

errosMoECenSN = erros |>
  tidyr::pivot_longer(paste0("X", 1:g), names_to = "grupo", values_to = "valor") |>
  dplyr::mutate(parametro_real = ifelse(grupo == "X1", parametros1, parametros2))

library(tidyverse)

errosMoECenSN |>
  filter(grupo == "X1") |>
  ggplot() +
  aes(
    x = as.factor(n),
    y = valor,
    fill = as.factor(cen)
  ) +
  geom_hline(aes(yintercept = parametro_real)) +
  geom_boxplot() +
  facet_wrap(~params, scales = "free") +
  labs(x = "Tamanho amostral") +
  guides(fill=guide_legend(title="Censura"))
