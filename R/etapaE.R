etapaE = function(Y, params, medias, args, ...){
  UseMethod("etapaE")
}

etapaE.MixNormal = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[j]*dnorm(y, medias[,j], params$params[j,"sigma"]))

  Z = t(sapply(1:n, function(i) Z[i,]/sum(Z[i,])))

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
