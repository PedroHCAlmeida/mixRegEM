estimaMedia = function(X, params, args, ...){
  UseMethod("estimaMedia")
}

estimaMedia.Normal = function(X, params, args){
  X%*%params[startsWith(names(params), "beta")]
}
.S3method("estimaMedia", "Normal", estimaMedia.Normal)

estimaMedia.MixNormal = function(X, params, args){
  sapply(1:args$g,
         function(j) estimaMedia.Normal(X, params[j,], args))}
.S3method("estimaMedia", "MixNormal", estimaMedia.MixNormal)

estimaMedia.MoENormal = function(X, params, args){
  sapply(1:args$g,
         function(j) estimaMedia.Normal(X, params[j,], args = args))
}
.S3method("estimaMedia", "MoENormal", estimaMedia.MoENormal)
