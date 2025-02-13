etapaM = function(y, X, U, params, args, ...){
  UseMethod("etapaM")
}

etapaM.MixNormal = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixNormal(y = y, X = X, Z = U$Z[,j], lasso = args$lasso)))


  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"sigma"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"sigma"]))/args$n
    }
  }

  P = colMeans(U$Z)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixNormal", etapaM.MixNormal)

etapaM.MoENormal = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    sapply(1:args$g,
           function(j) estimaTeta.MoENormal(y = y, X = X, Z = U$Z[,j], R = args$R,
                                            alpha = params$params[j,startsWith(colnames(params$params), "alpha")],
                                            P = params$P[,j], class = apply(U$Z, 1, which.max)==j, lasso = args$lasso)))

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"sigma"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"sigma"]))/args$n
    }
  }

  P = matrizP(t(paramsNovo)[startsWith(colnames(paramsNovo), "alpha"), 1:(args$g-1)], args$R)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MoENormal", etapaM.MoENormal)

etapaM.MoEKernelNormal = function(y, X, U, params, ...){

  paramsNovo = do.call(
    rbind,
    lapply(1:g,
           function(j) estimaTeta.MixNormal(y = y, X = X, Z = U$Z[,j])))

  kjList = lapply(1:g,
                  function(j){
                    sum_j = sum(U$Z[,j])
                    ks::kde(r, eval.points = r,
                            w = U$Z[,j]*n/sum_j)$estimate*sum_j/n
                  }
  )
  kj = do.call(cbind, kjList)
  kj_std = t(sapply(1:n,
                    function(i) kj[i,]/sum(kj[i,])))

  return(list(params = paramsNovo, P = kj_std))
}
.S3method("etapaM", "MoEKernelNormal", etapaM.MoEKernelNormal)

etapaM.MixT = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixT(y = y, X = X,
                                       Z = U$Z[,j], K = U$K[,j])))


  medias = estimaMedia(X, paramsNovo, args)
  Q = function(nu){
    sum(log(dMix.MixT(y, medias, paramsNovo[, startsWith(colnames(paramsNovo), "beta")],
              paramsNovo[,"sigma"], nu = nu, params$P)))
  }

  nu = optim(params$params[,"nu"],
             fn = Q,
             method = "L-BFGS-B",
             lower = 1,
             upper = 100,
             control = list(fnscale = -1)
             )$par

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }

  P = colMeans(U$Z)

  paramsNovo = as.matrix(cbind(paramsNovo, "nu" = as.numeric(nu)))

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixT", etapaM.MixT)

etapaM.MoET = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MoET(
             y = y,
             X = X,
             Z = U$Z[,j],
             K = U$K[,j],
             R = args$R,
             alpha = params$params[j,startsWith(colnames(params$params), "alpha")],
             P = params$P[,j]
             )
           )
    )


  Q = function(nu){
    sum(log(dMix.MixT(y, X, paramsNovo[, startsWith(colnames(paramsNovo), "beta")],
                      paramsNovo[,"sigma"], nu = nu, params$P)))
  }

  nu = optim(params$params[,"nu"],
             fn = Q,
             method = "L-BFGS-B",
             lower = 1,
             upper = 100,
             control = list(fnscale = -1)
  )$par

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }

  P = matrizP(t(paramsNovo)[startsWith(colnames(paramsNovo), "alpha"), 1:(args$g-1)], args$R)

  paramsNovo = as.matrix(cbind(paramsNovo, "nu" = as.numeric(nu)))

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MoET", etapaM.MoET)

etapaM.MoEST = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MoECenST(
             y = y,
             X = X,
             R = args$R,
             Z = U$Z[,j],
             e00 = U$e00[,j],
             e01 = U$e01[,j],
             e02 = U$e02[,j],
             e10 = U$e10[,j],
             e20 = U$e20[,j],
             e11 = U$e11[,j],
             delta = params$params[j,startsWith(colnames(params$params), "delta")],
             P = params$P[,j],
             alpha = params$params[j,startsWith(colnames(params$params), "alpha")],
             lambda = args$lambda[j],
             sigma = unname(params$params[j, "sigma"])
           )
    )
  )

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }


  if(any(is.nan(c(paramsNovo)))){
    return(params)
  } else{
    paramsAtual = paramsNovo
  }

  medias = estimaMedia(X, paramsNovo, args)
  Q = function(NU){
    if(length(NU) == 1) NU = rep(NU, args$g)

    ll = sum(log(dMix.MoEST(
      y = y,
      medias = medias,
      sigma = paramsNovo[,"sigma"],
      lambda = paramsNovo[,"lambda"],
      nu = NU,
      P = params$P,
      args = args
    )))
    if(!is.finite(ll)) -.Machine$double.xmax
    else ll
  }

  if(is.null(args$nuFixo)){
    if(is.null(args$nuIgual) || (args$nuIgual != T)){
      nu = optim(params$params[,"nu"],
                 fn = Q,
                 method = "L-BFGS-B",
                 lower = 1,
                 upper = 100,
                 control = list(fnscale = -1)
      )$par
    }else{
      nu = optimize(
        Q,
        c(1, 100),
        maximum = T
      )$maximum
      nu = rep(nu, args$g)
    }
  }
  else nu = params$params[,"nu"]

  paramsNovo = as.matrix(cbind(paramsNovo, "nu" = as.numeric(nu)))
  P = matrizP(t(paramsNovo)[startsWith(colnames(paramsNovo), "alpha"), 1:(args$g-1)], args$R)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MoEST", etapaM.MoEST)

etapaM.MixSN = function(y, X, U, params, args){

  medias = estimaMedia.MixSN(X, params$params, args)

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixSN(
             y = y,
             X = X,
             medias = medias[,j],
             Z = U$Z[,j],
             t1 = U$t1[,j],
             t2 = U$t2[,j],
             delta = params$params[j,"delta"],
             lasso = args$lasso
             )
           )
    )

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }


  P = colMeans(U$Z)
  if(args$Pequal){
    P = rep(1/args$g, args$g)
  }

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixSN", etapaM.MixSN)

etapaM.MixCenSN = function(y, X, U, params, args){
  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixCenSN(
             y = y,
             X = X,
             Z = U$Z[,j],
             e01 = U$e01[,j],
             e02 = U$e02[,j],
             e10 = U$e10[,j],
             e20 = U$e20[,j],
             e11 = U$e11[,j],
             delta = params$params[j,startsWith(colnames(params$params), "delta")],
             P = params$P[j],
             lambda = args$lambda[j],
             sigma = unname(params$params[j, "sigma"])
           )
    )
  )

  if(any(is.nan(c(paramsNovo)))){
    return(params)
  } else{
    paramsAtual = paramsNovo
  }

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }

  P = colMeans(U$Z)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixCenSN", etapaM.MixCenSN)

etapaM.MoECenSN = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MoECenSN(
             y = y,
             X = X,
             R = args$R,
             Z = U$Z[,j],
             e01 = U$e01[,j],
             e02 = U$e02[,j],
             e10 = U$e10[,j],
             e20 = U$e20[,j],
             e11 = U$e11[,j],
             delta = params$params[j,startsWith(colnames(params$params), "delta")],
             P = params$P[,j],
             alpha = params$params[j,startsWith(colnames(params$params), "alpha")],
             lambda = args$lambda[j],
             sigma = unname(params$params[j, "sigma"])
           )
    )
  )

  if(any(is.nan(c(paramsNovo)))){
    return(params)
  } else{
    paramsAtual = paramsNovo
  }

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }

  P = matrizP(t(paramsNovo)[startsWith(colnames(paramsNovo), "alpha"), 1:(args$g-1)], args$R)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MoECenSN", etapaM.MoECenSN)

etapaM.MixCenST = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixCenST(
             y = y,
             X = X,
             Z = U$Z[,j],
             e00 = U$e00[,j],
             e01 = U$e01[,j],
             e02 = U$e02[,j],
             e10 = U$e10[,j],
             e20 = U$e20[,j],
             e11 = U$e11[,j],
             delta = params$params[j,startsWith(colnames(params$params), "delta")],
             P = params$P[j],
             lambda = args$lambda[j],
             sigma = unname(params$params[j, "sigma"])
           )
    )
  )


  if(any(is.nan(c(paramsNovo)))){
    return(params)
  } else{
    paramsAtual = paramsNovo
  }

  medias = estimaMedia(X, paramsNovo, args)
  Q = function(NU){
    if(length(NU) == 1) NU = rep(NU, args$g)

    ll = sum(log(dMix.MixCenST(
      y = y,
      medias = medias,
      sigma = paramsNovo[,"sigma"],
      lambda = paramsNovo[,"lambda"],
      nu = NU,
      P = params$P,
      args = args
    )))
    if(!is.finite(ll)) -.Machine$double.xmax
    else ll
  }

  if(is.null(args$nuFixo)){
    if(is.null(args$nuIgual) || (args$nuIgual != T)){
      nu = optim(params$params[,"nu"],
                 fn = Q,
                 method = "L-BFGS-B",
                 lower = 1,
                 upper = 100,
                 control = list(fnscale = -1)
      )$par
    }else{
      nu = optimize(
        Q,
        c(1, 100),

        maximum = T
      )$maximum
      nu = rep(nu, args$g)
    }
  }
  else nu = params$params[,"nu"]

  paramsNovo = as.matrix(cbind(paramsNovo, "nu" = as.numeric(nu)))

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }

  P = colMeans(U$Z)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixCenST", etapaM.MixCenST)

etapaM.MoECenST = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MoECenST(
             y = y,
             X = X,
             R = args$R,
             Z = U$Z[,j],
             e00 = U$e00[,j],
             e01 = U$e01[,j],
             e02 = U$e02[,j],
             e10 = U$e10[,j],
             e20 = U$e20[,j],
             e11 = U$e11[,j],
             delta = params$params[j,startsWith(colnames(params$params), "delta")],
             P = params$P[,j],
             alpha = params$params[j,startsWith(colnames(params$params), "alpha")],
             lambda = args$lambda[j],
             sigma = unname(params$params[j, "sigma"])
           )
    )
  )


  if(any(is.nan(c(paramsNovo)))){
    return(params)
  } else{
    paramsAtual = paramsNovo
  }

  medias = estimaMedia(X, paramsNovo, args)
  Q = function(NU){
    if(length(NU) == 1) NU = rep(NU, args$g)

    ll = sum(log(dMix.MoECenST(
      y = y,
      medias = medias,
      sigma = paramsNovo[,"sigma"],
      lambda = paramsNovo[,"lambda"],
      nu = NU,
      P = params$P,
      args = args
    )))
    if(!is.finite(ll)) -.Machine$double.xmax
    else ll
  }

  if(is.null(args$nuFixo)){
    if(is.null(args$nuIgual) || (args$nuIgual != T)){
      nu = optim(params$params[,"nu"],
                 fn = Q,
                 method = "L-BFGS-B",
                 lower = 1,
                 upper = 100,
                 control = list(fnscale = -1)
      )$par
    }else{
      nu = optimize(
        Q,
        c(1, 100),

        maximum = T
      )$maximum
      nu = rep(nu, args$g)
    }
  }
  else nu = params$params[,"nu"]

  paramsNovo = as.matrix(cbind(paramsNovo, "nu" = as.numeric(nu)))

  if(!is.null(args$varEqual)){
    if(args$varEqual){
      paramsNovo[,"gama"] = sum(apply(U$Z, 1, function(x) x*paramsNovo[,"gama"]))/args$n
      paramsNovo[,"lambda"] = paramsNovo[,"delta"]/sqrt(paramsNovo[,"gama"])
      paramsNovo[,"sigma"] = sqrt(paramsNovo[,"delta"]**2 + paramsNovo[,"gama"])
    }
  }

  P = matrizP(t(paramsNovo)[startsWith(colnames(paramsNovo), "alpha"), 1:(args$g-1)], args$R)

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MoECenST", etapaM.MoECenST)


