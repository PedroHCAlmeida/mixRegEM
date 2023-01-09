estimaSe = function(X, params, args, weights = 1, ...){
  UseMethod("estimaSe")
}

estimaSe.Normal = function(X, params, args, weights = 1, ...){

  beta = params[startsWith(names(params), "beta")]
  sigma = params["sigma"]

  C = solve(t(weights*X)%*%(weights*X))

  se = sigma*sqrt(diag(C))

  t_0 = beta/se
  p_values = 2*pt(-abs(t_0), df = args$n-args$p)

  return(data.frame(list(se = se,
                         t_0 = t_0,
                         p_values = p_values))
  )
}
.S3method("estimaSe", "Normal", estimaSe.Normal)

estimaSe.MixNormal = function(X, params, args, U){
  lapply(1:args$g,
         function(j) estimaSe.Normal(X, params$params[j,], args = args, weights = U$Z[,j])
  )
}
.S3method("estimaSe.MixNormal", "MixNormal", estimaSe.MixNormal)

estimaSe.MoENormal = function(X, params, args, U){
  lapply(1:args$g,
         function(j) estimaSe.Normal(X, params$params[j,],
                                     args = args, weights = U$Z[,j])
  )
}
.S3method("estimaSe.MoENormal", "MoENormal", estimaSe.MoENormal)

estimaSe.MixT = function(X, params, args, U){
  lapply(1:args$g,
         function(j) estimaSe.Normal(X, params$params[j,], args = args, weights = U$Z[,j]*U$K[,j])
  )
}
.S3method("estimaSe.MixT", "MixT", estimaSe.MixT)
