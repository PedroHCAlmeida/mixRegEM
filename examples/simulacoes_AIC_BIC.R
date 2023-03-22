# Mistura dos experts ST com resposta censurada à esquerda

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

n = c(100, 200, 500, 1000)
nivelC = c(0, 0.075, 0.15, 0.3)
g = 3
tol = 1E-4
M = 100

beta01 <- c(0, -1, -2, -3)
alpha01 <- c(0.7, 1, 2)
sigma2_01 <- 1
lambda01 <- -1

beta02 <- c(-1, 1, 2, 3)
alpha02 <- c(1, 0, 3)
sigma2_02 <- 2
lambda02 <- 3

beta03 <- c(3, 5, -3, -1)
sigma2_03 <- 4
lambda03 <- 1

nu = c(2, 4, 6)

alpha <- matrix(c(alpha01, alpha02), byrow = T, nrow = 2)

rMoeEMSTCrit = function(ni, ci, tol = 1E-4, verbose = F){
  X <- cbind(rep(1, ni), runif(ni, 1, 5), runif(ni, -2, 2), runif(ni, 1, 4))
  R <- cbind(rep(1, ni), runif(ni, -2, 1), runif(ni, -1, 1))

  P <- matrizP2(alpha, R)

  mu01 = as.numeric(X%*%beta01)
  mu02 = as.numeric(X%*%beta02)
  mu03 = as.numeric(X%*%beta03)

  grupo <- numeric(ni)
  y <- numeric(ni)
  for(ii in 1:ni){
    arg1 <- list(mu = mu01[ii], sigma2 = sigma2_01, shape = lambda01, nu = nu[1])
    arg2 <- list(mu = mu02[ii], sigma2 = sigma2_02, shape = lambda02, nu = nu[2])
    arg3 <- list(mu = mu03[ii], sigma2 = sigma2_03, shape = lambda03, nu = nu[3])

    obs <- rmix(1, pii = P[ii,], family = 'Skew.t', arg = list(arg1, arg2, arg3), cluster = T)

    y[ii] <- obs$y
    grupo[ii] <- obs$cluster
  }

  c2 = as.numeric(quantile(y, probs = ci))
  y[y<=c2] = c2

  phi <- as.numeric(y == c2)

  critST = NULL
  tent = 0
  while(is.null(critST) & tent < 3){
    tent = tent+1
    try({
      critST = lapply(1:5, function(g){
        resultados = regEM(
          y,
          X[,-1],
          r = R[,-1],
          family = "MoECenST",
          phi = phi,
          c1 = -Inf,
          c2 = c2,
          g = g,
          showSE = F,
          verbose = verbose,
          tol = tol
      )
       c(resultados[["AIC"]], resultados[["BIC"]])
      })
    })
  }

  critSN = NULL
  tent = 0
  while(is.null(critSN) & tent < 3){
    tent = tent+1
    try({
      critSN = lapply(1:5, function(g){
        resultados = regEM(
          y,
          X[,-1],
          r = R[,-1],
          family = "MoECenSN",
          phi = phi,
          c1 = -Inf,
          c2 = c2,
          g = g,
          showSE = F,
          verbose = verbose,
          tol = tol
        )
        c(resultados[["AIC"]], resultados[["BIC"]])
      }
      )
    })
  }

  critT = NULL
  tent = 0
  while(is.null(critT) & tent < 3){
    tent = tent+1
    try({
      critT = lapply(1:5, function(g){
        resultados = regEM(
          y,
          X[,-1],
          r = R[,-1],
          family = "MoECenST",
          phi = phi,
          c1 = -Inf,
          c2 = c2,
          g = g,
          showSE = F,
          verbose = verbose,
          tol = tol,
          lambda = rep(0, g)
        )
        c(resultados[["AIC"]], resultados[["BIC"]])
      }
      )
    })
  }

  critN = NULL
  tent = 0
  while(is.null(critN) & tent < 3){
    tent = tent+1
    try({
      critN = lapply(1:5, function(g){
        resultados = regEM(
          y,
          X[,-1],
          r = R[,-1],
          family = "MoECenSN",
          phi = phi,
          c1 = -Inf,
          c2 = c2,
          g = g,
          showSE = F,
          verbose = verbose,
          tol = tol,
          lambda = rep(0, g)
        )
        c(resultados[["AIC"]], resultados[["BIC"]])
      }
      )
    })
  }

  if(!is.null(critST)){
    critST = as.data.frame(do.call(rbind, critST))
    critST$g = 1:5
    critST$model = "ST"
  }

  if(!is.null(critSN)){
    critSN = as.data.frame(do.call(rbind,critSN))
    critSN$g = 1:5
    critSN$model = "SN"
  }

  if(!is.null(critT)){
    critT = as.data.frame(do.call(rbind,critT))
    critT$g = 1:5
    critT$model = "T"
  }

  if(!is.null(critN)){
    critN = as.data.frame(do.call(rbind,critN))
    critN$g = 1:5
    critN$model = "N"
  }

  return(rbind(critST, critSN, critT, critN))
}

ntotal = length(n)*length(nivelC)
it = 0
start = Sys.time()
set.seed(123)
MoECenSTCrit = lapply(1:length(n), function(j) lapply(1:length(nivelC), function(i) matrix(rep(0, 4*5*5), ncol = 4)))
for(i in seq_along(n)){
  for(j in seq_along(nivelC)){
    MoECenSTCrit[[i]][[j]] = replicate(100, rMoeEMSTCrit(n[i], nivelC[j], tol = tol, verbose = F))
    cat("%s", it/ntotal)
    it = it+1
    save(MoECenSTCrit, file = "MoECenSTCrit.RData")
  }
}
end = Sys.time()
print(end-start)
