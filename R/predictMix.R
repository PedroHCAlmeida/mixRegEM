#' @param reg objeto de regressão de misturas
#' @export
predictMix = function(reg, x = NULL, r = NULL, class = T){
  UseMethod("predictMix")
}
#' @export
predictMix.MoENormal = function(reg, x, r, class = T){

  args = list()
  args$n = nrow(x)
  args$g = reg$g
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(x)), r)

  alpha = reg$Parametros[startsWith(rownames(reg$Parametros), "alpha"),]

  P = matrizP(alpha[,-ncol(alpha)], R)
  medias = estimaMedia.MoENormal(X, reg$params$params, args)
  if(!class) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:args$n,
               function(i) medias[i, grupos[i]])
  }
  return(y)
}
.S3method("predictMix", "MoENormal", predictMix.MoENormal)
