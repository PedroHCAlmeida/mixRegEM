chuteInicial = function(y, X, args, ...){
  UseMethod("chuteInicial")
}

chuteInicial.MixNormal = function(y, X, args, initGrupo = "KMeans"){

  dados = cbind(y, X[, -1])

  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T),
                  NULL = kmeans(dados, centers = args$g)$cluster
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F,
                                 MoreArgs = list(args = args)))

  P = prop.table(table(grupos))

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixNormal", chuteInicial.MixNormal)

chuteInicial.MoENormal = function(y, X, args, initGrupo = "KMeans"){

  k = ncol(args$R)
  dados = cbind(y, X[, -1])


  grupos = switch(initGrupo,
    "KMeans" = kmeans(dados, centers = args$g)$cluster,
    "Aleatório" = sample(1:args$g, args$n, replace = T),
    NULL = kmeans(dados, centers = args$g)$cluster
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F,
                                 MoreArgs = list(args = args)))

  P = matrix(rep(c(prop.table(table(grupos))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*k), rep(NA, k)), nrow = args$g, ncol = k, byrow = T)
  #alpha = matrix(c(nnet::multinom(grupos ~ R[,-1])$wts[-1], rep(NA, p)), nrow = g, ncol = k, byrow = T)
  colnames(alpha) = paste0("alpha", 1:k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoENormal", chuteInicial.MoENormal)

