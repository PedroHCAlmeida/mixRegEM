vero = function(Y, X, params, ...){
  UseMethod("vero")
}

vero.MixNormal = function(Y, X, params){
  sum(log(dMix.MixNormal(Y, X,
                         params$params[, startsWith(colnames(params$params), "beta")],
                         params$params[,"sigma"], params$P)))
}
.S3method("vero", "MixNormal", vero.MixNormal)

vero.MoENormal = function(Y, X, params){
  sum(log(dMix(Y, X,
               params$params[, startsWith(colnames(params$params), "beta")],
               params$params[,"sigma"], params$P)))
}
.S3method("vero", "MoENormal", vero.MoENormal)
