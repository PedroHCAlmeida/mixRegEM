#' @export predictMix
predictMix = function(reg, ...){
  UseMethod("predictMix")
}

predictMix.MoENormal = function(reg, x, r, type){

  args = list()
  args$n = nrow(x)

  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(x)), r)

  alpha = reg$Parametros[startsWith(rownames(reg$Parametros), "alpha"),]

  P = matrizP(alpha[,-ncol(alpha)], R)
  medias = estimaMedia.MoENormal(X, reg$params$params, args)
  if(type == 1) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:n,
               function(i) medias[i, grupos[i]])
  }
  return(y)
}
.S3method("predictMix", "MoENormal", predictMix.MoENormal)
