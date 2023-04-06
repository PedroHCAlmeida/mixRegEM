# Mistura dos experts ST com resposta censurada à esquerda
library(MomTrunc) # Sempre passar a sigma² como parâmetro nas funções
library(sn)
library(moments)
library(mixsmsn)
library(numDeriv)
library(mixRegEM)
library(future.apply)

matrizP2 <- function(alpha, R){
  P1 <- exp(R%*%t(alpha))/(1 + rowSums(exp(R%*%t(alpha))))
  P <- cbind(P1, 1 - rowSums(P1))
  return(P)
}

set.seed(123)

n = c(1000, 500)
nivelC = c(0, 0.075, 0.15, 0.3)
g = 2
tol = 1E-4

beta01 <- c(0, -1, -2, -3)
alpha01 <- c(0.7, 1, 2)
sigma2_01 <- 1
lambda01 <- -5

beta02 <- c(-1, 1, 2, 3)
sigma2_02 <- 2
lambda02 <- 3

alpha <- matrix(alpha01, byrow = T, nrow = 1)

rMoeEM = function(ni, ci, tol, verbose = F){

  files = list.files("R/")
  lapply(paste0("R/", files), source)


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
  resultados = NULL
  tryCatch({
    resultados = regEM(
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
      tol = tol,
      max_iter = 10000
    )
  }, error = function(e) print(e))

  if(is.null(resultados)){
    convergiu = F
  }
  else{
    convergiu = resultados$Convergiu
  }

  if(convergiu == F){
    return(NULL)
  }

  return(resultados$Parametros)
}

resultadosMoECenSN = vector(mode = 'list', length = length(n)) |>
  setNames(paste("n =", n))

plan(multisession, workers = 4)

set.seed(123)
for(ni in n){
  start = Sys.time()
  resultadosNi = list()

  for(ci in nivelC){
    resultadosNi = future_replicate(100, rMoeEM(ni, ci, tol = tol, verbose = F), future.seed = T)
    cat("Pronto: Censura ", ci, " Tamanho ", ni, "\n")
    resultadosMoECenSN[[paste("n =", ni)]][[paste("cen =", ci)]] = resultadosNi
    save(resultadosMoECenSN, file = "resultadosMoECenSN_10_2.RData")
  }
  sprintf("Pronto: %s", n)
  end = Sys.time()
  print(end-start)
  save(resultadosMoECenSN, file = "resultadosMoECenSN_10_2.RData")
}
