dMix = function(Y, X, params, ...){
  UseMethod("dMix")
}

dMix.MixNormal = function(y, X, beta, sigma, P){
  total = lapply(1:nrow(beta),
                 function(j) P[j]*dnorm(y, X%*%beta[j,], sigma[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MixNormal", dMix.MixNormal)

dMix.MoENormal = function(y, X, beta, sigma, P){
  total = lapply(1:nrow(beta),
                 function(j) P[,j]*dnorm(y, X%*%beta[j,], sigma[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MoENormal", dMix.MoENormal)
