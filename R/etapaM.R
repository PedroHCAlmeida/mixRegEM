etapaM = function(y, X, U, params, args, ...){
  UseMethod("etapaM")
}

etapaM.MixNormal = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixNormal(y = y, X = X, Z = U$Z[,j])))


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
                                            P = params$P[,j])))

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
  print(table(apply(kj_std, 1, which.max)))

  return(list(params = paramsNovo, P = kj_std))
}
.S3method("etapaM", "MoEKernelNormal", etapaM.MoEKernelNormal)

etapaM.MixT = function(y, X, U, params, args){

  paramsNovo = do.call(
    rbind,
    lapply(1:args$g,
           function(j) estimaTeta.MixT(y = y, X = X,
                                       Z = U$Z[,j], K = U$K[,j])))


  Q = function(nu){
    sum(log(dMix.MixT(y, X, paramsNovo[, startsWith(colnames(paramsNovo), "beta")],
              paramsNovo[,"sigma"], nu = nu, params$P)))
  }

  nu = optim(params$params[,"nu"],
             fn = Q,
             method = "L-BFGS-B",
             lower = 0.001,
             upper = 30,
             control = list(fnscale = -1)
             )$par

  P = colMeans(U$Z)

  paramsNovo = as.matrix(cbind(paramsNovo, "nu" = as.numeric(nu)))

  return(list(params = paramsNovo, P = P))
}
.S3method("etapaM", "MixT", etapaM.MixT)






