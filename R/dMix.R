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

dMix.MoET = function(y, medias, beta, sigma, nu, P){
  total = lapply(1:nrow(beta),
                 function(j) P[,j]*sn::dst(y, medias[j,], sigma[j], nu = nu[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MoET", dMix.MoET)

dMix.MixSN = function(y, medias, beta, sigma, lambda, delta, P){

  b = -sqrt(2/pi)

  total = lapply(1:nrow(beta),
                 function(j) P[j]*sn::dsn(y, medias[,j],
                                          omega = sigma[j],
                                          alpha = lambda[j]))
  return(rowSums(do.call(cbind, total)))
}
.S3method("dMix", "MixSN", dMix.MixSN)

dMix.MoESN = function(y, medias, beta, sigma, lambda, P, args){
  return(sapply(1:args$g,
                function(j){
                  P[,j]*sn::dsn(y, medias[,j], omega = sigma[j], alpha = lambda[j])

                }) |>
           rowSums())
}
.S3method("dMix", "MoESN", dMix.MoESN)

dMix.MoEST = function(y, medias, beta, sigma, lambda, nu, P, args){
  return(sapply(1:args$g,
                function(j){
                  sn::dst(x = y, xi = medias[, j], omega = sigma[j], alpha = lambda[j], nu = nu[j])

                }) |>
           rowSums())
}
.S3method("dMix", "MoEST", dMix.MoEST)

dMix.MixCenSN = function(y, medias, beta, sigma, lambda, P, args){

  phi1 = (args$phi[,1] == 1)
  return(sapply(1:args$g,
                function(j){
                  total = numeric(args$n)
                  total[!phi1] = P[j]*sn::dsn(y[!phi1], medias[!phi1,j], omega = sigma[j], alpha = lambda[j])
                  total[phi1] = P[j]*(sn::psn(args$c2,  medias[phi1,j], sigma[j], lambda[j])-sn::psn(args$c1, medias[phi1,j], sigma[j],  lambda[j]))
                  total
                }) |>
           rowSums())
}
.S3method("dMix", "MixCenSN", dMix.MixCenSN)

dMix.MoECenSN = function(y, medias, beta, sigma, lambda, P, args){

  phi1 = (args$phi[,1] == 1)
  return(sapply(1:args$g,
         function(j){
           total = numeric(args$n)
           total[!phi1] = P[!phi1,j]*sn::dsn(y[!phi1], medias[!phi1,j], omega = sigma[j], alpha = lambda[j])
           total[phi1] = P[phi1,j]*(sn::psn(args$c2,  medias[phi1,j], sigma[j], lambda[j])-sn::psn(args$c1, medias[phi1,j], sigma[j],  lambda[j]))
           total
           }) |>
    rowSums())
}
.S3method("dMix", "MoECenSN", dMix.MoECenSN)

dMix.MixCenST = function(y, medias, beta, sigma, lambda, nu, P, args){

  phi1 = (args$phi[,1] == 1)
  return(sapply(1:args$g,
                function(j){
                  total = numeric(args$n)
                  total[!phi1] = P[j]*sn::dst(x = y[!phi1], xi = medias[!phi1, j], omega = sigma[j], alpha = lambda[j], nu = nu[j])
                  total[phi1] = P[j]*(sn::pst(args$c2,  medias[phi1,j], sigma[j], lambda[j], nu[j])-sn::pst(args$c1, medias[phi1,j], sigma[j],  lambda[j], nu[j]))
                  total
                }) |>
           rowSums())
}
.S3method("dMix", "MixCenST", dMix.MixCenST)

dMix.MoECenST = function(y, medias, beta, sigma, lambda, nu, P, args){

  phi1 = (args$phi[,1] == 1)
  return(sapply(1:args$g,
                function(j){
                  total = numeric(args$n)
                  total[!phi1] = P[!phi1, j]*sn::dst(x = y[!phi1], xi = medias[!phi1, j], omega = sigma[j], alpha = lambda[j], nu = nu[j])
                  total[phi1] = P[phi1, j]*(sn::pst(args$c2,  medias[phi1,j], sigma[j], lambda[j], nu[j])-sn::pst(args$c1, medias[phi1,j], sigma[j],  lambda[j], nu[j]))
                  total
                }) |>
           rowSums())
}
.S3method("dMix", "MoECenST", dMix.MoECenST)

