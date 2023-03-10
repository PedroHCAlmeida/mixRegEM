dMix = function(Y, medias, params, ...){
  UseMethod("dMix")
}

dMix.MixNormal = function(y, medias, beta, sigma, P){
  total = lapply(1:nrow(beta),
                 function(j) P[j]*dnorm(y, medias[,j], sigma[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MixNormal", dMix.MixNormal)

dMix.MoENormal = function(y, medias, beta, sigma, P){
  total = lapply(1:nrow(beta),
                 function(j) P[,j]*dnorm(y, medias[,j], sigma[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MoENormal", dMix.MoENormal)

dMix.MixT = function(y, medias, beta, sigma, nu, P){
  total = lapply(1:nrow(beta),
                 function(j) P[j]*sn::dst(y, medias[,j], sigma[j], nu = nu[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MixT", dMix.MixT)

dMix.MoET = function(y, X, beta, sigma, nu, P){
  total = lapply(1:nrow(beta),
                 function(j) P[,j]*sn::dst(y, X%*%beta[j,], sigma[j], nu = nu[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MoET", dMix.MoET)

dMix.MixSN = function(y, medias, beta, sigma, lambda, delta, P){

  b = -sqrt(2/pi)

  total = lapply(1:nrow(beta),
                 function(j) P[j]*sn::dsn(y, medias[,j]+b*delta[j],
                                          omega = sigma[j],
                                          alpha = lambda[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MixSN", dMix.MixSN)

dMix.MoECenSN = function(y, medias, beta, sigma, lambda, P, args){

  phi1 = (args$phi == 1)
  return(sapply(1:args$g,
         function(j){
           total = numeric(args$n)
           total[!phi1] = P[!phi1,j]*sn::dsn(y[!phi1], medias[!phi1,j], sigma[j], lambda[j])
           total[phi1] = P[phi1,j]*(sn::psn(rep(args$c2, args$m),  medias[phi1,j], sigma[j], lambda[j])-sn::psn(rep(args$c1, args$m), medias[phi1,j], sigma[j],  lambda[j]))
           total
           }) |>
    rowSums())
}
.S3method("dMix", "MoECenSN", dMix.MoECenSN)

dMix.MoECenST = function(y, medias, beta, sigma, lambda, nu, P, args){

  phi1 = (args$phi == 1)
  return(sapply(1:args$g,
                function(j){
                  total = numeric(args$n)
                  total[!phi1] = P[!phi1, j]*sn::dst(x = y[!phi1], xi = medias[!phi1, j], omega = sigma[j], alpha = lambda[j], nu = nu[j])
                  total[phi1] = P[phi1, j]*(sn::pst(rep(args$c2, sum(phi1)),  medias[phi1,j], sigma[j], lambda[j], nu[j])-sn::pst(rep(args$c1, sum(phi1)), medias[phi1,j], sigma[j],  lambda[j], nu[j]))
                  total
                }) |>
           rowSums())
}
.S3method("dMix", "MoECenST", dMix.MoECenST)

