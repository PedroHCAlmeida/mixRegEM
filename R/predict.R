#' @export predictMix
predictMix = function(reg, ...){
  UseMethod("predictMix")
}

predictMix.MoENormal = function(reg, x, r, type){

  n = nrow(x)
  X = cbind(rep(1, nrow(x)), x)
  R = cbind(rep(1, nrow(x)), r)

  alpha = reg$Parâmetros %>%
    t() %>%
    as_tibble() %>%
    dplyr::select(starts_with("alpha")) %>%
    t()

  P = matrizP(alpha[,-ncol(alpha)], R)
  medias = estimaMedia.MoENormal(X, reg$params$params)
  if(type == 1) y = apply(P*medias, 1, sum)
  else{
    grupos = apply(P, 1, which.max)
    y = sapply(1:n,
               function(i) medias[i, grupos[i]])
  }
  return(y)
}
.S3method("predictMix", "MoENormal", predictMix.MoENormal)
