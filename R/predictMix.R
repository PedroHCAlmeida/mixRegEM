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

  args = list()
  args$n = nrow(x)
  args$g = reg$g
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(r)), r)

  X = MoENormal(X)

  alpha = reg$params$params[,startsWith(colnames(reg$params$params), "alpha")]

  P = matrizP(matrix(alpha[-nrow(alpha),], nrow=3, byrow = T), R)
  medias = estimaMedia(X, reg$params$params, args)

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

  args = list()
  args$n = nrow(x)
  args$g = reg$g
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(r)), r)

  X = MoENormal(X)

  alpha = reg$params$params[,startsWith(colnames(reg$params$params), "alpha")]

  P = matrizP(matrix(alpha[-nrow(alpha),], nrow=3, byrow = T), R)
  medias = estimaMedia(X, reg$params$params, args)

  if(!class) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  return(y)
}
.S3method("predictMix", "MoET", predictMix.MoET)
