#' @include classes.R
#' @include estimaMedia.R
#' @include auxFuncs.R

#' @name predictMix
#' @title predictMix
#' @param reg objeto de regressão de misturas
#' @param x variável explicativas
#' @param r variáveis explicativas
#' @param class classificar

#' @export
#' @return y
predictMix = function(reg, x = NULL, r = NULL, grupos = NULL, class = T, real = NULL, estimator = "mean"){
  UseMethod("predictMix")
}
#' @export
predictMix.MoENormal = function(reg, x, r, grupos = NULL, class = T, real = NULL){

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), as.matrix(x))
  R = cbind(rep(1, args$n), as.matrix(r))

  X = MoECenST(X)

  alpha = matrix(reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")],
                 nrow = args$g
  )

  P = matrizP(matrix(alpha[-args$g,], nrow=ncol(R), byrow = T), R)
  mu = estimaMedia(X, reg$Parametros, args)

  if(is.null(grupos)) grupos = apply(P, 1, which.max)

  if(!class) y = apply(P*mu, 1, sum)
  else{
    y = sapply(1:args$n,
               function(i) mu[i, grupos[i]])
  }
  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MoENormal", predictMix.MoENormal)

#' @export
predictMix.MoET = function(reg, x, r, class = T, real = NULL){

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), x)
  R = cbind(rep(1, args$n), r)

  X = MoENormal(X)

  alpha = reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")]

  P = matrizP(matrix(t(alpha[-nrow(alpha),]), nrow = ncol(R), byrow = T), R)
  medias = estimaMedia(X, reg$Parametros, args)

  if(!class) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MoET", predictMix.MoET)

predictMix.MixCenSN = function(reg, x, r, class = T, grupos = NULL, real = NULL){

  b = -sqrt(2/pi)

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), as.matrix(x))
  X = MixCenST(X)

  P = reg$P
  mu = estimaMedia(X, reg$Parametros, args)
  medias = sapply(
    1:args$g,
    function(i) mu[,i]-b*reg$Parametros[i,"delta"]
  )

  y = sapply(1:nrow(medias), function(i) sum(P*medias[i,]))

  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MixCenSN", predictMix.MixCenSN)

predictMix.MoECenSN = function(reg, x, r, class = T, grupos = NULL, real = NULL){

  b = -sqrt(2/pi)

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), as.matrix(x))
  R = cbind(rep(1, args$n), as.matrix(r))

  X = MoECenST(X)

  alpha = matrix(reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")],
                 nrow = args$g
  )

  P = matrizP(matrix(alpha[-args$g,], nrow=ncol(R), byrow = T), R)
  mu = estimaMedia(X, reg$Parametros, args)
  medias = sapply(
    1:args$g,
    function(i) mu[,i] -b*reg$Parametros[i,"delta"]
  )

  if(!class) y = apply(P*medias, 1, sum)
  else{
    if(is.null(grupos)) grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MoECenSN", predictMix.MoECenSN)

predictMix.MoECenST = function(reg, x, r, class = T, grupos = NULL, real = NULL){

  b = -sqrt(2/pi)

  reg$Parametros = t(reg$Parametros)
  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), as.matrix(x))
  R = cbind(rep(1, args$n), as.matrix(r))

  X = MoECenST(X)

  alpha = matrix(reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")],
                 nrow = args$g
  )

  P = matrizP(matrix(alpha[-args$g,], nrow=ncol(R), byrow = T), R)
  mu = estimaMedia(X, reg$Parametros, args)
  medias = sapply(
    1:args$g,
    function(i) mu[,i]-b*reg$Parametros[i,"delta"]
  )

  if(!class) y = apply(P*medias, 1, sum)
  else{
    if(is.null(grupos)) grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MoECenST", predictMix.MoECenST)


predictMix.MoECenSN = function(reg, x, r, class = T, grupos = NULL, real = NULL, estimator = "mean"){

  b = -sqrt(2/pi)

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), as.matrix(x))
  R = cbind(rep(1, args$n), as.matrix(r))

  X = MoECenST(X)

  alpha = matrix(reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")],
                 nrow = args$g
  )

  P = matrizP(matrix(alpha[-args$g,], nrow=ncol(R), byrow = T), R)
  mu = estimaMedia(X, reg$Parametros, args)
  medias = sapply(
    1:args$g,
    function(i){
      if(estimator == "mean") mu[,i]-b*reg$Parametros[i,"delta"]
      else mu[,i]
    }
  )

  if(is.null(grupos)) grupos = apply(P, 1, which.max)

  if(!class) y = apply(P*medias, 1, sum)
  else{
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MoECenSN", predictMix.MoECenSN)

predictMix.MoECenST = function(reg, x, r, class = T, grupos = NULL, real = NULL){

  b = -sqrt(2/pi)

  reg$Parametros = t(reg$Parametros)
  args = list()
  args$n = nrow(as.matrix(x))
  args$g = reg$g
  X = cbind(rep(1, args$n), as.matrix(x))
  R = cbind(rep(1, args$n), as.matrix(r))

  X = MoECenST(X)

  alpha = matrix(reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")],
                 nrow = args$g
                 )

  P = matrizP(matrix(alpha[-args$g,], nrow=ncol(R), byrow = T), R)
  mu = estimaMedia(X, reg$Parametros, args)
  medias = sapply(
    1:args$g,
    function(i) mu[,i]-b*reg$Parametros[i,"delta"]
  )

  if(!class) y = apply(P*medias, 1, sum)
  else{
    if(is.null(grupos)) grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  if(!is.null(real)){
    res = y-real
    rmse = sqrt(mean(res**2))
  } else{
    res = NULL
    rmse = NULL
  }

  return(list(y_hat = y, grupos = grupos, res = res, rmse = rmse))
}
.S3method("predictMix", "MoECenST", predictMix.MoECenST)

