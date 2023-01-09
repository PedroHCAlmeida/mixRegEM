etapaE = function(y, params, medias, args, ...){
  UseMethod("etapaE")
}

etapaE.MixNormal = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[j]*dnorm(y, medias[,j], params$params[j,"sigma"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))

  return(list(Z = Z))
}
.S3method("etapaE", "MixNormal", etapaE.MixNormal)

etapaE.MoENormal = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[,j]*dnorm(y, medias[,j], params$params[j,"sigma"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))
  return(list(Z = Z))
}
.S3method("etapaE", "MoENormal", etapaE.MoENormal)

etapaE.MixT = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[j]*sn::dst(y, X%*%params$params[j, startsWith(colnames(params$params), "beta")],
                                             params$params[j,"sigma"], nu = params$params[j,"nu"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))

  K = sapply(1:args$g,
             function(j) (params$params[j, "nu"]+1)/(params$params[j, "nu"]+dMahalanobis(y, medias[,j], params$params[j, "sigma"])))

  return(list(Z = Z, K = K))
}
.S3method("etapaE", "MixT", etapaE.MixT)

# microbenchmark(
#   "for" = {
#     Z = matrix(rep(0, args$n*args$g), ncol = args$g)
#     U = matrix(rep(0, args$n*args$g), ncol = args$g)
#     for(j in 1:args$g){
#       Z[, j] = params$P[j]*dt((y - medias[,j])/params$params[j,"sigma"], df = params$params[j,"nu"])
#       U[, j] = params$params[j, "nu"]+1/(params$params[j, "nu"]+dMahalanobis(y, medias[,j], params$params[j, "sigma"]))
#     }
#   },
#   "apply" = {
#     Z = sapply(1:args$g,
#                function(j) params$P[j]*dt((y - medias[,j])/params$params[j,"sigma"], df = params$params[j,"nu"]))
#     U = sapply(1:args$g,
#                function(j) params$params[j, "nu"]+1/(params$params[j, "nu"]+dMahalanobis(y, medias[,j], params$params[j, "sigma"])))
#   },
#    times = 100
#   )













