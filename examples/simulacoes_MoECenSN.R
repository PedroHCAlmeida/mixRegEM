# Mistura dos experts ST com resposta censurada à esquerda
library(MomTrunc) # Sempre passar a sigma² como parâmetro nas funções
library(sn)
library(moments)
library(mixsmsn)
library(numDeriv)
library(mixRegEM)

matrizP2 <- function(alpha, R){
  P1 <- exp(R%*%t(alpha))/(1 + rowSums(exp(R%*%t(alpha))))
  P <- cbind(P1, 1 - rowSums(P1))
  return(P)
}

set.seed(123)

n <- c(50, 100, 200, 500)
nivelC <- c(0.075, 0.15, 0.3)
g <- 2

beta01 <- c(0, -1, -2, -3)
alpha01 <- c(0.7, 1, 2)
sigma2_01 <- 1
lambda01 <- -1

beta02 <- c(-1, 1, 2, 3)
sigma2_02 <- 2
lambda02 <- 3

alpha <- matrix(alpha01, byrow = T, nrow = 1)

rMoeEM = function(ni, ci, tol = 1E-4, verbose = F){
  X <- cbind(rep(1, ni), runif(ni, 1, 5), runif(ni, -2, 2), runif(ni, 1, 4))
  R <- cbind(rep(1, ni), runif(ni, -2, 1), runif(ni, -1, 1))

  P <- matrizP2(alpha, R)

  mu01 <- as.numeric(X%*%beta01)
  mu02 <- as.numeric(X%*%beta02)

  grupo <- numeric(ni)
  y <- numeric(ni)
  for(ii in 1:ni){
    arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01)
    arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02)

    obs <- rmix(1, pii = P[ii,], family = 'Skew.normal', arg = list(arg1, arg2), cluster = T)

    y[ii] <- obs$y
    grupo[ii] <- obs$cluster
  }

  ki <- as.numeric(quantile(y, probs = ci))
  y[y <= ki] <- ki
  phi <- as.numeric(y == ki)
  resultados =
    regEM(
      y,
      X[,-1],
      r = R[,-1],
      family = "MoECenSN",
      phi = phi,
      c1 = -Inf,
      c2 = ki,
      g = 2,
      showSE = F,
      verbose = verbose,
      tol = tol
    )$Parametros
}

rMoeEM(500, 0.075, tol = 1E-6, verbose = T)

resultadosMoECenSN = n |>
  plyr::laply(
    function(ni){
      nivelC |>
        lapply(
          function(ci) replicate(500, rMoeEM(ni, ci))
        ) |>
        setNames(paste("Cen =", nivelC))
    }, .progress = "text"
  )

rownames(resultadosMoECenSN) = paste("n =", n)

parametros =cbind(c(
  beta01, sqrt(sigma2_01), lambda01, alpha01
),
c(
  beta02, sqrt(sigma2_02), lambda02, rep(NA, 3)
)
)

inverte_col = function(df){
  df_inv = df
  df_inv[,1] = df[,2]
  df_inv[,2] = df[,1]

  df_inv[startsWith(rownames(df), "alpha"),1] = -df_inv[startsWith(rownames(df), "alpha"),2]
  df_inv[startsWith(rownames(df), "alpha"),2] = NA

  df_inv
}

erros = n |>
  plyr::laply(
    function(ni){
      nivelC |>
        lapply(
          function(ci){
            lapply(
              1:dim(resultadosMoECenSN[paste("n =", ni), paste("Cen =", ci)][[1]])[3],
              function(i){
                x = resultadosMoECenSN[paste("n =", ni), paste("Cen =", ci)][[1]][,,i]

                res1 = x[startsWith(rownames(x), "beta"),1] - parametros[startsWith(rownames(parametros), "beta"),1]
                res2 = x[startsWith(rownames(x), "beta"),1] - parametros[startsWith(rownames(parametros), "beta"),2]

                if(sum(res1**2) > sum(res2**2)) x = inverte_col(x)

                x = x[!rownames(x) %in% c("delta", "gama"),]

                return(x - parametros)
              }
              )
          }
        ) |>
        setNames(paste("Cen =", nivelC))
    }, .progress = "text"
  )

rownames(erros) = paste("n =", n)

erros_MoECenSN = erros

erros = n |>
  plyr::laply(
    function(ni){
      nivelC |>
        lapply(
          function(ci){
            lapply(
              erros_MoECenSN[paste("n =", ni),paste("Cen =", ci)][[1]],
              function(x){
                params = rownames(x)
                x = data.frame(x)

                x$cen = ci
                x$n = ni
                x$params = params

                return(x)
              }
            )
          }
        )
    }, .progress = "text"
  ) |>
  dplyr::bind_rows() |>
  dplyr::as_tibble()

errosMoECenSN = erros |>
  tidyr::pivot_longer(paste0("X", 1:g), names_to = "grupo", values_to = "vies")


save(errosMoECenSN, file = "errosMoECenSN.RData")





















