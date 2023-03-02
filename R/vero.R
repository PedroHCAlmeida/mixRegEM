vero = function(y, medias, params, args, ...){
  UseMethod("vero")
}

vero.MixNormal = function(y, medias, params, args){
  sum(log(dMix.MixNormal(y, medias,
                         params$params[, startsWith(colnames(params$params), "beta")],
                         params$params[,"sigma"], params$P)))
}
.S3method("vero", "MixNormal", vero.MixNormal)

vero.MoENormal = function(y, medias, params, args){
  sum(log(dMix(y, medias,
               params$params[, startsWith(colnames(params$params), "beta")],
               params$params[,"sigma"], params$P)))
}
.S3method("vero", "MoENormal", vero.MoENormal)

vero.MixT = function(y, medias, params, args){
  sum(log(dMix(
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    nu = params$params[,"nu"],
    P = params$P
    )))
}
.S3method("vero", "MixT", vero.MixT)

vero.MoET = function(y, medias, params, args){
  sum(log(dMix(
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    nu = params$params[, "nu"],
    P = params$P
  )))
}
.S3method("vero", "MoET", vero.MoET)

vero.MixSN = function(y, X, params, args){
  sum(log(dMix.MixSN(
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    lambda = params$params[, "lambda"],
    delta = params$params[, "delta"],
    P = params$P)))
}
.S3method("vero", "MixSN", vero.MixSN)

vero.MoECenSN = function(y, medias, params, args){
  sum(log(dMix.MoECenSN(
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    lambda = params$params[, "lambda"],
    P = params$P,
    args = args
    )))
}
.S3method("vero", "MoECenSN", vero.MoECenSN)

vero.MoECenST = function(y, medias, params, args){
  sum(log(dMix.MoECenST(
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    lambda = params$params[, "lambda"],
    nu = params$params[, "nu"],
    P = params$P,
    args = args
  )))
}
.S3method("vero", "MoECenST", vero.MoECenST)
