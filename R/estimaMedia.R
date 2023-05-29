estimaMedia = function(f, X, params, args, ...){
  UseMethod("estimaMedia")
}

estimaMedia.Normal = function(f, X, params, args){
  X%*%params[startsWith(names(params), "beta")]
}
.S3method("estimaMedia", "Normal", estimaMedia.Normal)

estimaMedia.MixNormal = function(f, X, params, args){
  do.call(
    cbind,
    lapply(1:args$g,
           function(j) estimaMedia.Normal(f, X, params[j,], args))
  )
}
.S3method("estimaMedia", "MixNormal", estimaMedia.MixNormal)

estimaMedia.MoENormal = function(f, X, params, args){
  do.call(
    cbind,
    lapply(1:args$g,
         function(j) estimaMedia.Normal(f, X, params[j,], args = args)))
}
.S3method("estimaMedia", "MoENormal", estimaMedia.MoENormal)

estimaMedia.MixT = function(f, X, params, args){
  do.call(
    cbind,
    lapply(1:args$g,
         function(j) estimaMedia.Normal(f, X, params[j,], args)))
}
.S3method("estimaMedia", "MixT", estimaMedia.MixT)

estimaMedia.MoET = function(f, X, params, args){
  do.call(
    cbind,
    lapply(1:args$g,
         function(j) estimaMedia.Normal(f, X, params[j,], args = args)))
}
.S3method("estimaMedia", "MoET", estimaMedia.MoET)

estimaMedia.MixSN = function(f, X, params, args){
  do.call(
    cbind,
    lapply(1:args$g,
         function(j) estimaMedia.Normal(f, X, params[j,], args = args)))
}
.S3method("estimaMedia", "MixSN", estimaMedia.MixSN)



