matrizP <- function(alpha, R){

  if(any(is.na(alpha))){
    P = as.matrix(rep(1, nrow(R)))
  } else{
    exp_R_alpha = as.matrix(exp(R%*%alpha))
    P1 <- exp_R_alpha/(1 + rowSums(exp_R_alpha))
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

