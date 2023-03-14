chuteInicial = function(y, X, args, ...){
  UseMethod("chuteInicial")
}

chuteInicial.MixNormal = function(y, X, args){

  dados = cbind(y, X[, -1])

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  P = prop.table(table(grupos))

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixNormal", chuteInicial.MixNormal)

chuteInicial.MoENormal = function(y, X, args, initGrupo = "KMeans"){

  k = ncol(args$R)
  dados = cbind(y, X[, -1])


  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  P = matrix(rep(c(prop.table(table(grupos))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*k), rep(NA, k)), nrow = args$g, ncol = k, byrow = T)
  #alpha = matrix(c(nnet::multinom(grupos ~ R[,-1])$wts[-1], rep(NA, p)), nrow = g, ncol = k, byrow = T)
  colnames(alpha) = paste0("alpha", 1:k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoENormal", chuteInicial.MoENormal)

chuteInicial.MixT = function(y, X, args, initGrupo = "KMeans"){

  dados = cbind(y, X[, -1])

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos),
                                                  matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  params = as.matrix(cbind(params, "nu" = rep(5, args$g)))

  P = prop.table(table(grupos))

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixT", chuteInicial.MixT)

chuteInicial.MoET = function(y, X, args, initGrupo = "KMeans"){

  k = ncol(args$R)
  #dados = cbind(y, X[, -1])

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(X[,-1], centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos),
                                                  matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  params = as.matrix(cbind(params, "nu" = rep(5, args$g)))

  P = matrix(rep(c(prop.table(table(grupos))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*k), rep(NA, k)),
                 nrow = args$g, ncol = k, byrow = T)

  colnames(alpha) = paste0("alpha", 1:k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoET", chuteInicial.MoET)

chuteInicial.MixSN = function(y, X, args){

  dados = cbind(y, X[, -1])

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  lambda = sapply(
    1:args$g,
    function(j)
      moments::skewness(y[grupos == j]-medias[grupos == j, j])
    )

  params = cbind(
    params,
    "lambda" = lambda,
    "delta" = params[, "sigma"]*(lambda/sqrt(1 + lambda**2)),
    "gama" = (params[, "sigma"]**2)*(1 - (lambda/sqrt(1 + lambda**2))**2)
  )
  P = prop.table(table(grupos))

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixSN", chuteInicial.MixSN)

chuteInicial.MoECenSN = function(y, X, args){

  dados = cbind(y, X[, -1])

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  lambda = sapply(
    1:args$g,
    function(j)
      moments::skewness(y[grupos == j]-medias[grupos == j, j])
  )

  params = cbind(
    params,
    "lambda" = lambda,
    "delta" = params[, "sigma"]*(lambda/sqrt(1 + lambda**2)),
    "gama" = (params[, "sigma"]**2)*(1 - (lambda/sqrt(1 + lambda**2))**2)
  )

  P = matrix(rep(c(prop.table(table(grupos))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*args$k), rep(NA, args$k)),
                 nrow = args$g, ncol = args$k, byrow = T)

  colnames(alpha) = paste0("alpha", 1:args$k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoECenSN", chuteInicial.MoECenSN)

chuteInicial.MoECenST = function(y, X, args){

  # if(is.null(args$nuFixo)){
  #
  #   if(is.null(args$nu)){
  #     nuFixo = rep(5, args$g)
  #   } else{
  #     nuFixo = args$nu
  #   }

  regSN = regEM(
    y,
    X[,-1],
    family = "MoECenSN",
    phi = args$phi,
    c1 = args$c1,
    c2 = args$c1,
    g = args$g,
    showSE = F,
    verbose = F,
    tol = 1E-4,
    max_iter = 50
  )

    params = t(regSN$Parametros)
    P = regSN$P

    medias = X %*% t(params[, startsWith(colnames(params), "beta")])

    Q = function(NU){
      sum(log(dMix.MoECenST(
        y = y,
        medias = medias,
        sigma = params[,"sigma"],
        lambda = params[,"lambda"],
        nu = NU,
        P = P,
        args = args
      )))
    }

    nu = optim(rep(30, args$g),
               fn = Q,
               method = "L-BFGS-B",
               lower = 1,
               upper = 30,
               control = list(fnscale = -1)
    )$par

    paramsNovo = as.matrix(cbind(params, "nu" = as.numeric(nu)))

    return(list(params = params, P = P))
  # } else{
  #   chuteInicial.MoECenSTNuFixo(y, X, args)
  # }
}
.S3method("chuteInicial", "MoECenST", chuteInicial.MoECenST)

chuteInicial.MoECenSTNuFixo = function(y, X, args){

  dados = cbind(y)

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(dados, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  lambda = sapply(
    1:args$g,
    function(j)
      moments::skewness(y[grupos == j]-medias[grupos == j, j])
  )

  params = cbind(
    params,
    "lambda" = lambda,
    "delta" = params[, "sigma"]*(lambda/sqrt(1 + lambda**2)),
    "gama" = (params[, "sigma"]**2)*(1 - (lambda/sqrt(1 + lambda**2))**2),
    "nu" = ifelse(is.null(args$nu), rep(5, args$g), args$nu)
  )

  P = matrix(rep(c(prop.table(table(grupos))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*args$k), rep(NA, args$k)),
                 nrow = args$g, ncol = args$k, byrow = T)

  colnames(alpha) = paste0("alpha", 1:args$k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoECenST", chuteInicial.MoECenST)

