estimaEp = function(X, params, args, weights = 1, ...){
  UseMethod("estimaEp")
}

estimaEp.Normal = function(X, params, args, weights = 1, ...){

  beta = params[startsWith(names(params), "beta")]
  sigma = params["sigma"]

  C = solve(t(weights*X)%*%(weights*X))

  Eps = sigma*sqrt(diag(C))

  t_0 = beta/Eps
  p_values = 2*pt(-abs(t_0), df = args$n-args$p)

  return(list(Eps = Eps,
              t_0 = t_0,
              p_values = p_values
  ))
}
.S3method("estimaEp", "Normal", estimaEp.Normal)

estimaEp.MixNormal = function(X, params, args, U){
  lapply(1:args$g,
         function(j) estimaEp.Normal(X, params$params[j,], args = args, weights = U$Z[,j])
  )
}
.S3method("estimaEp.MixNormal", "MixNormal", estimaEp.MixNormal)

estimaEp.MoENormal = function(X, params, args, U){
  lapply(1:args$g,
         function(j) estimaEp.Normal(X, params$params[j,],
                                     args = args, weights = U$Z[,j])
  )
}
.S3method("estimaEp.MoENormal", "MoENormal", estimaEp.MoENormal)
