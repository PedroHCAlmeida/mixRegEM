chuteInicial = function(y, X, args, ...){
  UseMethod("chuteInicial")
}

chuteInicial.MixNormal = function(y, X, args){

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = do.call(function(...) kmeans(y, centers = args$g, ...), args$cluster_args)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.MixNormal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 Z = 1,
                                 lasso = ifelse(is.null(args$lasso), F, args$lasso),
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  if(args$chuteMahalanobis){
    try({
      if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
      if(is.null(args$tol_chute)) args$tol_chute = 0.01
      diff_grupos = 1
      i = 1
      while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
        i = i+1
        dma = do.call(
          cbind,
          lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
        )

        grupos_new = apply(dma, 1, which.min)
        diff_grupos = 1-mean(grupos_new == grupos_novo)

        if(min(table(grupos_new))/args$n <= 0.1){
          i = args$maxChuteIter
          break
        }

        dadosGrupos = lapply(list("X" = X, "y" = y),
                             function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                             grupos = grupos_new)

        params = do.call(rbind, mapply(estimaTeta.MixNormal,
                                       dadosGrupos$y,
                                       dadosGrupos$X,
                                       Z = 1,
                                       lasso = ifelse(is.null(args$lasso), F, args$lasso),
                                       SIMPLIFY = F))

        medias = estimaMedia(X, params, args)
        grupos_novo = grupos_new
      }
    })
  }

  P = prop.table(table(grupos))

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixNormal", chuteInicial.MixNormal)

chuteInicial.MoENormal = function(y, X, args, initGrupo = "KMeans"){

  k = ncol(args$R)
  dados = cbind(y)


  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = do.call(function(...) kmeans(dados, centers = args$g, ...), args$cluster_args)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]


  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.MixNormal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 Z = 1,
                                 lasso = ifelse(is.null(args$lasso), F, args$lasso),
                                 SIMPLIFY = F))


  medias = estimaMedia(X, params, args)

  if(args$chuteMahalanobis){

    if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
    if(is.null(args$tol_chute)) args$tol_chute = 0.01
    diff_grupos = 1
    i = 1
    while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
      i = i+1
      dma = do.call(
        cbind,
        lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
      )

      grupos_new = apply(dma, 1, which.min)
      diff_grupos = 1-mean(grupos_new == grupos_novo)

      if(min(table(grupos_new))/args$n <= 0.1){
        i = args$maxChuteIter
        break
      }

      dadosGrupos = lapply(list("X" = X, "y" = y),
                           function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                           grupos = grupos_new)

      params = do.call(rbind, mapply(estimaTeta.MixNormal,
                                     dadosGrupos$y,
                                     dadosGrupos$X,
                                     Z = 1,
                                     lasso = ifelse(is.null(args$lasso), F, args$lasso),
                                     SIMPLIFY = F))

      medias = estimaMedia(X, params, args)
      grupos_novo = grupos_new
    }
  }


  P = matrix(rep(c(prop.table(table(grupos))), args$n), byrow = T, ncol = args$g)
  #alpha = matrix(c(nnet::multinom(grupos ~ R[,-1])$wts[-1], rep(NA, p)), nrow = g, ncol = k, byrow = T)
  alpha = matrix(c(rep(0, (args$g-1)*k), rep(NA, k)), nrow = args$g, ncol = k, byrow = T)

  colnames(alpha) = paste0("alpha", 1:k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoENormal", chuteInicial.MoENormal)

chuteInicial.MixT = function(y, X, args, initGrupo = "KMeans"){


  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = do.call(function(...) kmeans(y, centers = args$g, ...), args$cluster_args)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos),
                                                  matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  if(args$chuteMahalanobis){

    if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
    if(is.null(args$tol_chute)) args$tol_chute = 0.01
    diff_grupos = 1
    i = 1
    while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
      i = i+1
      dma = do.call(
        cbind,
        lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
      )

      grupos_new = apply(dma, 1, which.min)
      diff_grupos = 1-mean(grupos_new == grupos_novo)

      if(min(table(grupos_new))/args$n <= 0.1){
        i = args$maxChuteIter
        break
      }

      dadosGrupos = lapply(list("X" = X, "y" = y),
                           function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                           grupos = grupos_new)

      params = do.call(rbind, mapply(estimaTeta.Normal,
                                     dadosGrupos$y,
                                     dadosGrupos$X,
                                     SIMPLIFY = F))


      medias = estimaMedia(X, params, args)
      grupos_novo = grupos_new
    }
  }

  params = as.matrix(cbind(params, "nu" = rep(5, args$g)))

  P = prop.table(table(grupos_novo))

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixT", chuteInicial.MixT)

chuteInicial.MoET = function(y, X, args, initGrupo = "KMeans"){

  k = ncol(args$R)

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = kmeans(y, centers = args$g)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos),
                                                  matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  # if(args$chuteMahalanobis){
  #
  #   if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
  #   if(is.null(args$tol_chute)) args$tol_chute = 0.01
  #   diff_grupos = 1
  #   i = 1
  #   while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
  #     i = i+1
  #     dma = do.call(
  #       cbind,
  #       lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
  #     )
  #
  #     grupos_new = apply(dma, 1, which.min)
  #     diff_grupos = 1-mean(grupos_new == grupos_novo)
  #
  #     if(min(table(grupos_new))/args$n <= 0.1){
  #       i = args$maxChuteIter
  #       break
  #     }
  #
  #     dadosGrupos = lapply(list("X" = X, "y" = y),
  #                          function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
  #                          grupos = grupos_new)
  #
  #     params = do.call(rbind, mapply(estimaTeta.Normal,
  #                                    dadosGrupos$y,
  #                                    dadosGrupos$X,
  #                                    SIMPLIFY = F))
  #
  #
  #     medias = estimaMedia(X, params, args)
  #     grupos_novo = grupos_new
  #   }
  # }

  params = as.matrix(cbind(params, "nu" = rep(1, args$g)))

  P = matrix(rep(c(prop.table(table(grupos_novo))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*k), rep(NA, k)),
                 nrow = args$g, ncol = k, byrow = T)

  colnames(alpha) = paste0("alpha", 1:k)
  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoET", chuteInicial.MoET)

chuteInicial.MixSN = function(y, X, args){

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = do.call(function(...) kmeans(y, centers = args$g, ...), args$cluster_args)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  if(args$chuteMahalanobis){

    if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
    if(is.null(args$tol_chute)) args$tol_chute = 0.01
    diff_grupos = 1
    i = 1
    while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
      i = i+1
      dma = do.call(
        cbind,
        lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
      )

      grupos_new = apply(dma, 1, which.min)
      diff_grupos = 1-mean(grupos_new == grupos_novo)

      if(min(table(grupos_new))/args$n <= 0.1){
        i = args$maxChuteIter
        break
      }

      dadosGrupos = lapply(list("X" = X, "y" = y),
                           function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                           grupos = grupos_new)

      params = do.call(rbind, mapply(estimaTeta.Normal,
                                     dadosGrupos$y,
                                     dadosGrupos$X,
                                     SIMPLIFY = F))


      medias = estimaMedia(X, params, args)
      grupos_novo = grupos_new
    }
  }

  lambda = sapply(
    1:args$g,
    function(j)
      moments::skewness(y[grupos_novo == j]-medias[grupos_novo == j, j])
    )

  params = cbind(
    params,
    "lambda" = lambda,
    "delta" = params[, "sigma"]*(lambda/sqrt(1 + lambda**2)),
    "gama" = (params[, "sigma"]**2)*(1 - (lambda/sqrt(1 + lambda**2))**2)
  )
  P = prop.table(table(grupos_novo))
  if(args$Pequal){
    P = rep(1/args$g, args$g)
  }

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixSN", chuteInicial.MixSN)

chuteInicial.MixCenSN = function(y, X, args){

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = do.call(function(...) kmeans(y, centers = args$g, ...), args$cluster_args)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )

  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  if(args$chuteMahalanobis){

    if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
    if(is.null(args$tol_chute)) args$tol_chute = 0.01
    diff_grupos = 1
    i = 1
    while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
      i = i+1
      dma = do.call(
        cbind,
        lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
      )

      grupos_new = apply(dma, 1, which.min)
      diff_grupos = 1-mean(grupos_new == grupos_novo)

      if(min(table(grupos_new))/args$n <= 0.1){
        i = args$maxChuteIter
        break
      }

      dadosGrupos = lapply(list("X" = X, "y" = y),
                           function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                           grupos = grupos_new)

      params = do.call(rbind, mapply(estimaTeta.Normal,
                                     dadosGrupos$y,
                                     dadosGrupos$X,
                                     SIMPLIFY = F))

      medias = estimaMedia(X, params, args)
      grupos_novo = grupos_new
    }
  }

  P = prop.table(table(grupos_novo))

  if(is.null(args$lambda)){
    lambda = sapply(
      1:args$g,
      function(j)
        moments::skewness(y[grupos_novo == j]-medias[grupos_novo == j, j])
    )
  }
  else{
    lambda = args$lambda
  }


  params = cbind(
    params,
    "lambda" = lambda,
    "delta" = params[, "sigma"]*(lambda/sqrt(1 + lambda**2)),
    "gama" = (params[, "sigma"]**2)*(1 - (lambda/sqrt(1 + lambda**2))**2)
  )

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixCenSN", chuteInicial.MixCenSN)

chuteInicial.MoECenSN = function(y, X, args){

  if(is.null(args$initGrupo)) args$initGrupo = "KMeans"
  grupos = switch(args$initGrupo,
                  "KMeans" = do.call(function(...) kmeans(y, centers = args$g, ...), args$cluster_args)$cluster,
                  "Aleatório" = sample(1:args$g, args$n, replace = T)
  )
  props = prop.table(table(grupos))
  ordem = order(-props)
  grupos_novo = ordem[grupos]

  dadosGrupos = lapply(list("X" = X, "y" = y),
                       function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                       grupos = grupos_novo)

  params = do.call(rbind, mapply(estimaTeta.Normal,
                                 dadosGrupos$y,
                                 dadosGrupos$X,
                                 SIMPLIFY = F))

  medias = estimaMedia(X, params, args)

  if(args$chuteMahalanobis){

    if(is.null(args$maxChuteIter)) args$maxChuteIter = 10
    if(is.null(args$tol_chute)) args$tol_chute = 0.01
    diff_grupos = 1
    i = 1
    while(diff_grupos > args$tol_chute & i < args$maxChuteIter){
      i = i+1
      dma = do.call(
        cbind,
        lapply(1:args$g, function(i) dMahalanobis(y, medias[,i], params[i, "sigma"]))
      )

      grupos_new = apply(dma, 1, which.min)
      diff_grupos = 1-mean(grupos_new == grupos_novo)

      if(min(table(grupos_new))/args$n <= 0.1){
        i = args$maxChuteIter
        break
      }

      dadosGrupos = lapply(list("X" = X, "y" = y),
                           function(x, grupos) lapply(split(x, grupos), matrix, ncol=dim(as.matrix(x))[2]),
                           grupos = grupos_new)

      params = do.call(rbind, mapply(estimaTeta.Normal,
                                     dadosGrupos$y,
                                     dadosGrupos$X,
                                     SIMPLIFY = F))

      medias = estimaMedia(X, params, args)
      grupos_novo = grupos_new
    }
  }



  P = matrix(rep(c(prop.table(table(grupos_novo))), args$n), byrow = T, ncol = args$g)
  alpha = matrix(c(rep(0, (args$g-1)*args$k), rep(NA, args$k)),
                 nrow = args$g, ncol = args$k, byrow = T)

  colnames(alpha) = paste0("alpha", 1:args$k)

  if(is.null(args$lambda)){
    lambda = sapply(
      1:args$g,
      function(j)
        moments::skewness(y[grupos_novo == j]-medias[grupos_novo == j, j])
    )
  }
  else{
    lambda = args$lambda
  }

  params = cbind(
    params,
    "lambda" = lambda,
    "delta" = params[, "sigma"]*(lambda/sqrt(1 + lambda**2)),
    "gama" = (params[, "sigma"]**2)*(1 - (lambda/sqrt(1 + lambda**2))**2)
  )

  params = cbind(params, alpha = alpha)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoECenSN", chuteInicial.MoECenSN)

chuteInicial.MoEST = function(y, X, args){

  if(is.null(args$initial_iter)) args$initial_iter = 2

  class(y) = ""
  modeloSN = regEM(
    y,
    X[, -1],
    r = args$R[,-1],
    g = args$g,
    min_iter = args$initial_iter,
    max_iter = args$initial_iter,
    family = "MoESN",
    lambda = args$lambda,
    verbose = args$verbose,
    tol = args$tol,
    mcFirst = F,
    chuteMahalanobis = args$chuteMahalanobis,
    initGrupo = args$initGrupo,
    varEqual = args$varEqual
  )

  params = modeloSN$Parametros
  if(args$g>1) params = t(params)
  medias = estimaMedia(X, params, args)
  P = modeloSN$P

  Q = function(NU){

    if(length(NU) == 1) NU = rep(NU, args$g)

    ll = sum(log(dMix.MoEST(
      y = y,
      medias = medias,
      sigma = params[,"sigma"],
      lambda = params[,"lambda"],
      nu = NU,
      P = P,
      args = args
    )))
    if(!is.finite(ll)) -.Machine$double.xmax
    else ll
  }

  nu = tryCatch({
    if(is.null(args$nuFixo)){
      if(is.null(args$nuIgual) || (args$nuIgual != T)){
        nu = optim(c(3, 3),
                   fn = Q,
                   method = "L-BFGS-B",
                   lower = 2.01,
                   upper = 150,
                   abstol = 1E-6,
                   control = list(fnscale = -1)
        )$par
      }else{
        nu = optimize(
          Q,
          c(1, 30),
          tol = 1E-6,
          maximum = T
        )$maximum
        nu = rep(nu, args$g)
      }
    }
    else nu = args$nuFixo
  },
  error = function(e) {return(rep(3, args$g))}
  )

  if(!is.null(args$nuFixo)) nu = args$nuFixo
  params = cbind(params, nu = nu)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoEST", chuteInicial.MoEST)

chuteInicial.MixCenST = function(y, X, args){

  if(is.null(args$initial_iter)) args$initial_iter = 2

  modeloSN = regEM(
    y,
    X[, -1],
    phi = args$phi,
    c1 = args$c1,
    c2 = args$c2,
    g = args$g,
    min_iter = args$initial_iter,
    max_iter = args$initial_iter,
    family = "MixCenSN",
    lambda = args$lambda,
    verbose = args$verbose,
    tol = args$tol,
    initGrupo = args$initGrupo,
    mcFirst = F
  )

  params = modeloSN$Parametros
  if(args$g>1) params = t(params)
  medias = estimaMedia(X, params, args)
  P = modeloSN$P

  Q = function(NU){

    if(length(NU) == 1) NU = rep(NU, args$g)

    ll = sum(log(dMix.MixCenST(
      y = y,
      medias = medias,
      sigma = params[,"sigma"],
      lambda = params[,"lambda"],
      nu = NU,
      P = P,
      args = args
    )))
    if(!is.finite(ll)) -.Machine$double.xmax
    else ll
  }

  nu = tryCatch({
    if(is.null(args$nuFixo)){
      if(is.null(args$nuIgual) || (args$nuIgual != T)){
        nu = optim(rep(30, args$g),
                   fn = Q,
                   method = "L-BFGS-B",
                   lower = 2.01,
                   upper = 150,
                   abstol = 1E-6,
                   control = list(fnscale = -1)
        )$par
      }else{
        nu = optimize(
          Q,
          tol = 1E-6,
          c(2.01, 150),
          maximum = T
        )$maximum
        nu = rep(nu, args$g)
      }
    }
    else nu = args$nuFixo
  },
  error = function(e) {return(rep(3, args$g))}
  )

  if(!is.null(args$nuFixo)) nu = args$nuFixo
  params = cbind(params, nu = nu)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MixCenST", chuteInicial.MixCenST)

chuteInicial.MoECenST = function(y, X, args){

  if(is.null(args$initial_iter)) args$initial_iter = 2

  modeloSN = regEM(
    y,
    X[, -1],
    r = args$R[,-1],
    phi = args$phi,
    c1 = args$c1,
    c2 = args$c2,
    g = args$g,
    min_iter = args$initial_iter,
    max_iter = args$initial_iter,
    family = "MoECenSN",
    lambda = args$lambda,
    verbose = args$verbose,
    initGrupo = args$initGrupo,
    mcFirst = F
  )

  params = modeloSN$Parametros
  if(args$g>1) params = t(params)
  medias = estimaMedia(X, params, args)

  P = modeloSN$P

  Q = function(NU){

    if(length(NU) == 1) NU = rep(NU, args$g)

    ll = sum(log(dMix.MoECenST(
      y = y,
      medias = medias,
      sigma = params[,"sigma"],
      lambda = params[,"lambda"],
      nu = NU,
      P = P,
      args = args
    )))
    if(!is.finite(ll)) -.Machine$double.xmax
    else ll
  }

  nu = tryCatch({
    if(is.null(args$nuFixo)){
      if(is.null(args$nuIgual) || (args$nuIgual != T)){
        nu = optim(c(3, 3),
                   fn = Q,
                   method = "L-BFGS-B",
                   lower = 2.01,
                   upper = 150,
                   abstol = 1E-6,
                   control = list(fnscale = -1)
        )$par
      }else{
        nu = optimize(
          Q,
          c(2.01, 150),
          tol = 1E-6,
          maximum = T
        )$maximum
        nu = rep(nu, args$g)
      }
    }
    else nu = args$nuFixo
  },
  error = function(e) return(rep(3, args$g))
  )

  if(!is.null(args$nuFixo)) nu = args$nuFixo
  params = cbind(params, nu = nu)

  return(list(params = params, P = P))
}
.S3method("chuteInicial", "MoECenST", chuteInicial.MoECenST)

