vero = function(f, y, medias, params, args, ...){
  UseMethod("vero")
}

vero.MixNormal = function(f, y, medias, params, args){
  sum(log(dMix.MixNormal(f, y, medias,
                         params$params[, startsWith(colnames(params$params), "beta")],
                         params$params[,"sigma"], params$P)))
}
.S3method("vero", "MixNormal", vero.MixNormal)

vero.MoENormal = function(f, y, medias, params, args){
  sum(log(dMix(f, y, medias,
               params$params[, startsWith(colnames(params$params), "beta")],
               params$params[,"sigma"], params$P)))
}
.S3method("vero", "MoENormal", vero.MoENormal)

vero.MixT = function(f, y, medias, params, args){
  sum(log(dMix(
    f = f,
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    nu = params$params[,"nu"],
    P = params$P
    )))
}
.S3method("vero", "MixT", vero.MixT)

vero.MoET = function(f, y, medias, params, args){
  sum(log(dMix(
    f = f,
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    nu = params$params[, "nu"],
    P = params$P
  )))
}
.S3method("vero", "MoET", vero.MoET)

vero.MixSN = function(f, y, X, params, args){
  sum(log(dMix.MixSN(
    f = f,
    y = y,
    medias = medias,
    beta = params$params[, startsWith(colnames(params$params), "beta")],
    sigma = params$params[,"sigma"],
    lambda = params$params[, "lambda"],
    delta = params$params[, "delta"],
    P = params$P)))
}
.S3method("vero", "MixSN", vero.MixSN)

vero.MoECenSN = function(f, y, medias, params, args){
  sum(log(dMix.MoECenSN(
    f = f,
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

vero.MoECenST = function(f, y, medias, params, args){
  sum(log(dMix.MoECenST(
    f = f,
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
