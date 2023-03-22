matrizP <- function(alpha, R){

  if(any(is.na(alpha))){
    P = as.matrix(rep(1, nrow(R)))
  } else{
    P1 <- exp(R%*%alpha)/(1 + rowSums(exp(R%*%alpha)))
    P <- cbind(P1, 1 - rowSums(P1))
  }
  return(P)
}

calculaMetricas = function(y, medias){
  mse = mean((y-medias)**2)
  rmse = sqrt(mse)
  mae = mean(abs(y-medias))
  return(list(n = length(y), MSE = mse, RMSE = rmse, MAE = mae))
}

dMahalanobis = function(y, mu, sigma){
  ((y-mu)/sigma)**2
}

