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

n = c(100, 200, 500, 1000)
nivelC = c(0, 0.075, 0.15, 0.3)
g = 2
tol = 1E-4

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

  nc = floor(ni*ci)+1
  ind_censored = sort(sample(1:ni, nc, replace = F))
  u = runif(nc)

  c1 = mapply(function(ic, i) max(c(y[ic] - u[i], y[ic]+u[i]-1)), ind_censored, 1:length(ind_censored))
  c2 = mapply(function(ic, i) min(c(y[ic] + u[i], y[ic]-u[i]+1)), ind_censored, 1:length(ind_censored))

  for(i in 1:length(ind_censored)){
    y[ind_censored[i]] = (c1[i]+c2[i])/2
  }
  phi <- as.numeric(1:length(y) %in% ind_censored)

  resultados = NULL
  tent = 0
  while(is.null(resultados) & tent < 3){
    tent = tent+1
    resultados = regEM(
      y,
      X[,-1],
      r = R[,-1],
      family = "MoECenSN",
      phi = phi,
      c1 = c1,
      c2 = c2,
      g = 2,
      showSE = F,
      verbose = verbose,
      tol = tol
    )$Parametros
  }
  return(resultados)
}

start = Sys.time()
set.seed(123)
resultadosMoECenSN_int = n |>
  plyr::laply(
    function(ni){
      nivelC |>
        lapply(
          function(ci) replicate(500, rMoeEM(ni, ci, tol = tol))
        ) |>
        setNames(paste("Cen =", nivelC))
    }, .progress = "text"
  )
end = Sys.time()
print(end-start)

rownames(resultadosMoECenSN_int) = paste("n =", n)

parametros =cbind(c(
  beta01, sqrt(sigma2_01), lambda01, alpha01
),
c(
  beta02, sqrt(sigma2_02), lambda02, rep(NA, 3)
)
)

rownames(parametros) = c(paste0("beta", 1:4), "sigma", "lambda", paste0("alpha", 1:3))

rownames(resultadosMoECenSN) = paste("n =", n)


inverte_col = function(df){
  df_inv = df
  df_inv[,1] = df[,2]
  df_inv[,2] = df[,1]

  df_inv[startsWith(rownames(df), "alpha"),1] = -df_inv[startsWith(rownames(df), "alpha"),2]
  df_inv[startsWith(rownames(df), "alpha"),2] = NA

  df_inv
}

erros_int = n |>
  plyr::laply(
    function(ni){
      nivelC |>
        lapply(
          function(ci){
            apply(
              resultadosMoECenSN_int[paste("n =", ni), paste("Cen =", ci)][[1]],
              3,
              function(x){

                res1 = x[startsWith(rownames(x), "beta"),1] - parametros[startsWith(rownames(parametros), "beta"),1]
                res2 = x[startsWith(rownames(x), "beta"),1] - parametros[startsWith(rownames(parametros), "beta"),2]

                if(sum(res1**2) > sum(res2**2)) x = inverte_col(x)

                params = rownames(x)
                x = x[!params %in% c("delta", "gama"),]

                #x = x - parametros

                x = data.frame(x)
                x$parametros1 = parametros[,1]
                x$parametros2 = parametros[,2]
                x$n = ni
                x$cen = ci
                x$params = params[!params %in% c("delta", "gama")]

                return(x)
              }
            )
          }
        )
    }, .progress = "text"
  ) |>
  dplyr::bind_rows() |>
  dplyr::as_tibble()

erros_intMoECenSN = erros_int |>
  tidyr::pivot_longer(paste0("X", 1:2), names_to = "grupo", values_to = "valor") |>
  dplyr::mutate(parametro_real = ifelse(grupo == "X1", parametros1, parametros2))


save(erros_intMoECenSN, file = "erros_int_MoECenSN_intervalar.RData")













