matrizP <- function(alpha, R){
  P1 <- exp(R%*%alpha)/(1 + rowSums(exp(R%*%alpha)))
  P <- cbind(P1, 1 - rowSums(P1))
  return(P)
}

calculaMetricas = function(y, medias){
  mse = mean((y-medias)**2)
  rmse = sqrt(mse)
  mae = mean(abs(y-medias))
  return(list(n = length(y), MSE = mse, RMSE = rmse, MAE = mae))
}
