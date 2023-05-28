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
predictMix = function(reg, x = NULL, r = NULL, class = T){
  UseMethod("predictMix")
}
#' @export
predictMix.MoENormal = function(reg, x, r, class = T){

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(x)
  args$g = reg$g
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(r)), r)

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
  return(y)
}
.S3method("predictMix", "MoENormal", predictMix.MoENormal)

#' @export
predictMix.MoET = function(reg, x, r, class = T){

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(x)
  args$g = reg$g
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(r)), r)

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
  return(y)
}
.S3method("predictMix", "MoET", predictMix.MoET)

predictMix.MoECenST = function(reg, x, r, class = T){

  reg$Parametros = t(reg$Parametros)

  args = list()
  args$n = nrow(x)
  args$g = reg$g
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(r)), r)

  X = MoECenST(X)

  alpha = reg$Parametros[,startsWith(colnames(reg$Parametros), "alpha")]

  P = matrizP(matrix(alpha[-nrow(alpha),], nrow=ncol(R), byrow = T), R)
  medias = estimaMedia(X, reg$Parametros, args)

  if(!class) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  return(y)
}
.S3method("predictMix", "MoECenST", predictMix.MoECenST)
