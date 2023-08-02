etapaE = function(y, params, medias, args, ...){
  UseMethod("etapaE")
}

etapaE.MixNormal = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[j]*dnorm(y, medias[,j], params$params[j,"sigma"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))

  return(list(Z = Z))
}
.S3method("etapaE", "MixNormal", etapaE.MixNormal)

etapaE.MoENormal = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[,j]*dnorm(y, medias[,j], params$params[j,"sigma"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))
  return(list(Z = Z))
}
.S3method("etapaE", "MoENormal", etapaE.MoENormal)

etapaE.MixT = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[j]*sn::dst(y, X%*%params$params[j, startsWith(colnames(params$params), "beta")],
                                             params$params[j,"sigma"], nu = params$params[j,"nu"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))

  K = sapply(1:args$g,
             function(j) (params$params[j, "nu"]+1)/(params$params[j, "nu"]+dMahalanobis(y, medias[,j], params$params[j, "sigma"])))

  return(list(Z = Z, K = K))
}
.S3method("etapaE", "MixT", etapaE.MixT)

etapaE.MoET = function(y, X, params, medias, args, ...){

  Z = sapply(1:args$g,
             function(j) params$P[,j]*sn::dst(y, X%*%params$params[j, startsWith(colnames(params$params), "beta")],
                                              params$params[j,"sigma"], nu = params$params[j,"nu"]))

  Z = t(sapply(1:args$n, function(i) Z[i,]/sum(Z[i,])))

  K = sapply(1:args$g,
             function(j) (params$params[j, "nu"]+1)/(params$params[j, "nu"]+dMahalanobis(y, medias[,j], params$params[j, "sigma"])))

  return(list(Z = Z, K = K))
}
.S3method("etapaE", "MoET", etapaE.MoET)

etapaE.MixSN = function(y, X, params, medias, args, ...){

  b = -sqrt(2/pi)

  U_list = sapply(1:args$g,
             function(j){
               Z = params$P[j]*sn::dsn(y, medias[,j], params$params[j,"sigma"], params$params[j,"lambda"])
               Mu = (params$params[j,"delta"]/(params$params[j,"delta"]**2+params$params[j,"gama"]))*(y-medias[,j]-b*params$params[j,"delta"])
               M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
               MT = sqrt(M2T)
               p = pnorm(Mu/MT)
               p = ifelse(p == 0, .Machine$double.xmin, p)
               tau = (dnorm(Mu/MT))/p
               t1 = (Mu + b + MT*tau)
               t2 = ((Mu + b)^2 + M2T + (Mu + 2*b)*MT*tau)
               return(list(Z, t1, t2))
               })


  U = lapply(1:3, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "t1", "t2"))

  U$Z = t(sapply(1:args$n, function(i) U$Z[i,]/sum(U$Z[i,])))
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)
  U$t1 = U$Z*U$t1
  U$t2 = U$Z*U$t2

  return(U)
}
.S3method("etapaE", "MixSN", etapaE.MixSN)

etapaE.MoESN = function(y, X, params, medias, args, ...){

  medias = estimaMedia(X, params$params, args)

  U_list = sapply(
    1:args$g,
    function(j){
      Z = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)

      Z = params$P[, j]*sn::dsn(y, medias[, j], params$params[j,"sigma"], params$params[j,"lambda"])

      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)
      Mu = M2T*params$params[j,"delta"]*(y-medias[,j])/params$params[j,"gama"]

      e01 = y
      e02 = y**2

      aij = Mu/MT
      p = pnorm(aij)
      p = ifelse(p == 0, .Machine$double.xmin, p)
      tau_gama = (dnorm(aij))/p

      e10 = (Mu + MT*tau_gama)
      e20 = Mu**2 + M2T + Mu*MT*tau_gama
      e11 = e01*e10

      return(list(Z, e01, e02, e10, e20, e11))
    })
  U = lapply(1:6, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e01", "e02", "e10", "e20", "e11"))

  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)
  return(U)
}
.S3method("etapaE", "MoESN", etapaE.MoESN)

etapaE.MixST = function(y, X, params, medias, args, ...){

  U_list = sapply(
    1:args$g,
    function(j){

      Z = e00 = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)
      d2ij = dMahalanobis(y, medias[, j], params$params[j,"sigma"])

      Z = params$P[j]*sn::dst(y, medias[, j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])

      MuAux = (params$params[j,"delta"]/(params$params[j,"delta"]**2+params$params[j,"gama"]))
      Mu = MuAux*(y-medias[,j])
      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)

      aij = params$params[j,"lambda"]*(y - medias[,j])/params$params[j,"sigma"]

      e00Aux = 4*(params$params[j,"nu"]**(params$params[j,"nu"]/2))*gamma((params$params[j,"nu"]+3)/2)*
        ((params$params[j,"nu"]+d2ij)**(-(params$params[j,"nu"]+3)/2))/
        (sqrt(pi)*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2))*
        pt(aij*sqrt((params$params[j,"nu"] + 3)/(params$params[j,"nu"] + d2ij)), params$params[j,"nu"] + 3)

      pU = sn::dst(y, medias[,j], params$params[j,"sigma"], params$params[j,"lambda"], params$params[j,"nu"])

      e00 = e00Aux/ifelse(pU == 0, .Machine$double.xmin, pU)
      e01 = e00*y

      e02 = e00*(y**2)

      aij = params$params[j,"lambda"]*(y - medias[,j])/params$params[j,"sigma"]

      auxU = (2*MT*params$params[j,"nu"]**(params$params[j,"nu"]/2)*gamma((params$params[j,"nu"]+2)/2)*
                (params$params[j,"nu"]+d2ij+(aij**2))**(-(params$params[j,"nu"]+2)/2))/
        (pi*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2)*pU)

      cv = 2*gamma((params$params[j,"nu"]+1)/2)/(gamma((params$params[j,"nu"])/2)*sqrt(params$params[j,"nu"]*(1+(params$params[j,"lambda"]**2))*pi))

      e10 = (Mu*e00) + auxU
      e20 = ((Mu**2)*e00) + M2T + Mu*auxU
      e11 = e10*y

      return(list(Z, e00, e01, e02, e10, e20, e11))
    })

  U = lapply(1:7, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e00", "e01", "e02", "e10", "e20", "e11"))
  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)

  return(U)
}
.S3method("etapaE", "MixST", etapaE.MixST)

etapaE.MoEST = function(y, X, params, medias, args, ...){

  U_list = sapply(
    1:args$g,
    function(j){

      Z = e00 = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)
      d2ij = dMahalanobis(y, medias[, j], params$params[j,"sigma"])

      Z = params$P[, j]*sn::dst(y, medias[, j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])

      MuAux = (params$params[j,"delta"]/(params$params[j,"delta"]**2+params$params[j,"gama"]))
      Mu = MuAux*(y-medias[,j])
      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)

      aij = params$params[j,"lambda"]*(y - medias[,j])/params$params[j,"sigma"]

      e00Aux = 4*(params$params[j,"nu"]**(params$params[j,"nu"]/2))*gamma((params$params[j,"nu"]+3)/2)*
        ((params$params[j,"nu"]+d2ij)**(-(params$params[j,"nu"]+3)/2))/
        (sqrt(pi)*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2))*
        pt(aij*sqrt((params$params[j,"nu"] + 3)/(params$params[j,"nu"] + d2ij)), params$params[j,"nu"] + 3)

      pU = sn::dst(y, medias[,j], params$params[j,"sigma"], params$params[j,"lambda"], params$params[j,"nu"])

      e00 = e00Aux/ifelse(pU == 0, .Machine$double.xmin, pU)
      e01 = e00*y

      e02 = e00*(y**2)

      aij = params$params[j,"lambda"]*(y - medias[,j])/params$params[j,"sigma"]

      auxU = (2*MT*params$params[j,"nu"]**(params$params[j,"nu"]/2)*gamma((params$params[j,"nu"]+2)/2)*
                (params$params[j,"nu"]+d2ij+(aij**2))**(-(params$params[j,"nu"]+2)/2))/
        (pi*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2)*pU)

      cv = 2*gamma((params$params[j,"nu"]+1)/2)/(gamma((params$params[j,"nu"])/2)*sqrt(params$params[j,"nu"]*(1+(params$params[j,"lambda"]**2))*pi))

      e10 = (Mu*e00) + auxU
      e20 = ((Mu**2)*e00) + M2T + Mu*auxU
      e11 = e10*y

      return(list(Z, e00, e01, e02, e10, e20, e11))
    })

  U = lapply(1:7, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e00", "e01", "e02", "e10", "e20", "e11"))
  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)

  return(U)
}
.S3method("etapaE", "MoEST", etapaE.MoEST)

etapaE.MixCenSN = function(y, X, params, medias, args, ...){

  phi1 = (args$phi[,1] == 1)

  medias = estimaMedia(X, params$params, args)

  U_list = sapply(
    1:args$g,
    function(j){
      Z = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)
      R0 = sn::psn(args$c2, medias[phi1, j], params$params[j,"sigma"], params$params[j,"lambda"])-sn::psn(args$c1, medias[phi1,j], params$params[j,"sigma"], params$params[j,"lambda"])
      R0 = ifelse(R0 == 0, .Machine$double.xmin, R0)

      Z[!phi1] = params$P[j]*sn::dsn(y[!phi1], medias[!phi1, j], params$params[j,"sigma"], params$params[j,"lambda"])
      Z[phi1] = params$P[j]*R0


      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)
      Mu = (y-medias[,j])*params$params[j,"delta"]/(params$params[j,"delta"]**2+params$params[j,"gama"])

      moments = mapply(
        function(ic, i){
          if(params$params[j,"lambda"] == 0){
            return(MomTrunc::meanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                        Sigma = params$params[j,"sigma"]^2,
                                        dist = "normal")[c(1,2)])
          }else{
            return(MomTrunc::meanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                        Sigma = params$params[j,"sigma"]^2,
                                        lambda = params$params[j,"lambda"], dist = "SN")[c(1,2)])
          }}
        , which(phi1), 1:args$m
      )

      if(length(moments) == 0) moments = rbind(0, 0)


      e01[!phi1] = y[!phi1]
      try({
        e01[phi1] = unlist(moments[1,])
        e02[phi1] = unlist(moments[2,])
      }, silent = T)
      e02[!phi1] = y[!phi1]**2

      w0[phi1] = mapply(
        function(ic, i) MomTrunc::meanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                             Sigma = params$params[j,"gama"], dist = 'normal')$mean
        , which(phi1), 1:args$m
      )

      aij = params$params[j,"lambda"]*(y[!phi1] - medias[!phi1,j])/params$params[j,"sigma"]
      p = pnorm(aij)
      p = ifelse(p == 0, .Machine$double.xmin, p)
      tau_gama[!phi1] = (dnorm(aij))/p

      try({
        P0 = pnorm(args$c2, medias[phi1,j], sqrt(params$params[j,"gama"]))-pnorm(args$c1, medias[phi1,j], sqrt(params$params[j,"gama"]))
        tau_gama[phi1] = P0/(sqrt(pi*(1 + params$params[j,"lambda"]^2)/2)*R0)
      }, silent = T)

      e10[!phi1] = (Mu[!phi1] + MT*tau_gama[!phi1])
      try({
        e10[phi1] = ((M2T*params$params[j,"delta"]/params$params[j,"gama"])*(e01[phi1] - medias[phi1, j]))+MT*tau_gama[phi1]
      }, silent = T)
      e20[!phi1] = Mu[!phi1]**2 + M2T + Mu[!phi1]*MT*tau_gama[!phi1]
      try({
        e20[phi1] = (M2T*params$params[j,"delta"]/params$params[j,"gama"])**2*(e02[phi1] - 2*medias[phi1, j]*e01[phi1] + medias[phi1, j]**2) +
          (MT**3*params$params[j,"delta"]/params$params[j,"gama"]*(w0[phi1]-medias[phi1, j])*tau_gama[phi1]) + M2T
      }, silent = T)
      e11[!phi1] = e01[!phi1]*e10[!phi1]
      try({
        e11[phi1] = (M2T*params$params[j,"delta"]/params$params[j,"gama"])*(e02[phi1] -medias[phi1, j]*e01[phi1]) +
          (MT*w0[phi1]*tau_gama[phi1])
      }, silent = T)

      return(list(Z, e01, e02, e10, e20, e11))
    })
  U = lapply(1:6, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e01", "e02", "e10", "e20", "e11"))

  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)

  return(U)
}
.S3method("etapaE", "MixCenSN", etapaE.MixCenSN)

etapaE.MoECenSN = function(y, X, params, medias, args, ...){

  phi1 = (args$phi[,1] == 1)

  medias = estimaMedia(X, params$params, args)

  U_list = sapply(
    1:args$g,
    function(j){
      Z = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)
      R0 = sn::psn(args$c2, medias[phi1, j], params$params[j,"sigma"], params$params[j,"lambda"])-sn::psn(args$c1, medias[phi1,j], params$params[j,"sigma"], params$params[j,"lambda"])
      R0 = ifelse(R0 == 0, .Machine$double.xmin, R0)

      Z[!phi1] = params$P[!phi1, j]*sn::dsn(y[!phi1], medias[!phi1, j], params$params[j,"sigma"], params$params[j,"lambda"])
      Z[phi1] = params$P[phi1, j]*R0


      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)
      Mu = M2T*params$params[j,"delta"]*(y-medias[,j])/params$params[j,"gama"]

      moments = mapply(
        function(ic, i){
          if(params$params[j,"lambda"] == 0){
            return(MomTrunc::meanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                 Sigma = params$params[j,"sigma"]^2,
                                 dist = "normal")[c(1,2)])
          }else{
            return(MomTrunc::meanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                               Sigma = params$params[j,"sigma"]^2,
                               lambda = params$params[j,"lambda"], dist = "SN")[c(1,2)])
          }}
        , which(phi1), 1:args$m
      )

      if(length(moments) == 0) moments = rbind(0, 0)


      e01[!phi1] = y[!phi1]
      try({
        e01[phi1] = unlist(moments[1,])
        e02[phi1] = unlist(moments[2,])
      }, silent = T)
      e02[!phi1] = y[!phi1]**2

      w0[phi1] = mapply(
        function(ic, i) MomTrunc::meanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                             Sigma = params$params[j,"gama"], dist = 'normal')$mean
        , which(phi1), 1:args$m
      )

      aij = params$params[j,"lambda"]*(y[!phi1] - medias[!phi1,j])/params$params[j,"sigma"]
      p = pnorm(aij)
      p = ifelse(p == 0, .Machine$double.xmin, p)
      tau_gama[!phi1] = (dnorm(aij))/p

      try({
        P0 = pnorm(args$c2, medias[phi1,j], sqrt(params$params[j,"gama"]))-pnorm(args$c1, medias[phi1,j], sqrt(params$params[j,"gama"]))
        tau_gama[phi1] = P0/(sqrt(pi*(1 + params$params[j,"lambda"]^2)/2)*R0)
      }, silent = T)

      e10[!phi1] = (Mu[!phi1] + MT*tau_gama[!phi1])
      try({
        e10[phi1] = ((M2T*params$params[j,"delta"]/params$params[j,"gama"])*(e01[phi1] - medias[phi1, j]))+MT*tau_gama[phi1]
      }, silent = T)
      e20[!phi1] = Mu[!phi1]**2 + M2T + Mu[!phi1]*MT*tau_gama[!phi1]
      try({
        e20[phi1] = (M2T*params$params[j,"delta"]/params$params[j,"gama"])**2*(e02[phi1] - 2*medias[phi1]*e01[phi1] + medias[phi1, j]**2) +
              (MT**3*params$params[j,"delta"]/params$params[j,"gama"]*(w0[phi1]-medias[phi1, j])*tau_gama[phi1]) + M2T
      }, silent = T)
      e11[!phi1] = e01[!phi1]*e10[!phi1]
      try({
        e11[phi1] = (M2T*params$params[j,"delta"]/params$params[j,"gama"])*(e02[phi1] -medias[phi1, j]*e01[phi1]) +
                    (MT*w0[phi1]*tau_gama[phi1])
      }, silent = T)

      return(list(Z, e01, e02, e10, e20, e11))
      })
  U = lapply(1:6, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e01", "e02", "e10", "e20", "e11"))

  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)
  return(U)
}
.S3method("etapaE", "MoECenSN", etapaE.MoECenSN)

etapaE.MixCenST = function(y, X, params, medias, args, ...){

  phi1 = (args$phi[,1] == 1)

  U_list = sapply(
    1:args$g,
    function(j){

      Z = e00 = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)

      F0 = sn::pst(args$c2, medias[phi1, j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])-
        sn::pst(args$c1, medias[phi1,j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])

      d2ij = dMahalanobis(y, medias[, j], params$params[j,"sigma"])

      Z[!phi1] = params$P[j]*sn::dst(y[!phi1], medias[!phi1, j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])
      Z[phi1] = params$P[j]*F0

      MuAux = (params$params[j,"delta"]/(params$params[j,"delta"]**2+params$params[j,"gama"]))
      Mu = MuAux*(y-medias[,j])
      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)

      aij = params$params[j,"lambda"]*(y[!phi1] - medias[!phi1,j])/params$params[j,"sigma"]

      e00Aux = 4*(params$params[j,"nu"]**(params$params[j,"nu"]/2))*gamma((params$params[j,"nu"]+3)/2)*
        ((params$params[j,"nu"]+d2ij[!phi1])**(-(params$params[j,"nu"]+3)/2))/
        (sqrt(pi)*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2))*
        pt(aij*sqrt((params$params[j,"nu"] + 3)/(params$params[j,"nu"] + d2ij[!phi1])), params$params[j,"nu"] + 3)

      pU = sn::dst(y[!phi1], medias[!phi1,j], params$params[j,"sigma"], params$params[j,"lambda"], params$params[j,"nu"])

      sigma_ =  (params$params[j,"sigma"]**2)*(params$params[j,"nu"]/(params$params[j,"nu"]+2))
      P0 = sn::pst(args$c2, medias[phi1, j], sqrt(sigma_), params$params[j,"lambda"], nu = params$params[j,"nu"]+2)-
        sn::pst(args$c1, medias[phi1,j], sqrt(sigma_), params$params[j,"lambda"], nu = params$params[j,"nu"]+2)

      sigma__ = (params$params[j,"sigma"]**2)*sqrt(params$params[j,"nu"]/((params$params[j,"nu"]+1)*(1+params$params[j,"lambda"]**2)))

      R0 = sn::pst(args$c2, medias[phi1, j], sqrt(sigma__),  0, params$params[j,"nu"]+1)-
        sn::pst(args$c1, medias[phi1,j], sqrt(sigma__), 0, params$params[j,"nu"]+1)

      R0_F0 = R0/ifelse(F0 == 0, .Machine$double.xmin, F0)

      e00[!phi1] = e00Aux/ifelse(pU == 0, .Machine$double.xmin, pU)
      e00[phi1] = P0/ifelse(F0 == 0, .Machine$double.xmin, F0)

      moments = mapply(
        function(ic, i){
          if(params$params[j,"lambda"] == 0){
            return(MomTrunc::MCmeanvarTMD(
              lower = args$c1[i],
              upper = args$c2[i],
              mu = medias[ic, j],
              Sigma = as.matrix(sigma_),
              nu = as.matrix(round(params$params[j,"nu"])+2),
              dist = 't')[c(1,2)])
          }else{
            return(MomTrunc::MCmeanvarTMD(
              lower = args$c1[i],
              upper = args$c2[i],
              mu = medias[ic, j],
              Sigma = as.matrix(sigma_),
              lambda = as.matrix(params$params[j,"lambda"]),
              nu = as.matrix(round(params$params[j,"nu"])+2),
              dist = 'ST')[c(1,2)])
          }}
        , which(phi1), 1:args$m
      )

      if(length(moments) == 0) moments = rbind(0, 0)

      e01[!phi1] = e00[!phi1]*y[!phi1]
      e01[phi1] = e00[phi1]*unlist(moments[1,])

      e02[!phi1] = e00[!phi1]*(y[!phi1]**2)
      e02[phi1] = e00[phi1]*unlist(moments[2,])

      w0[phi1] = mapply(
        function(ic, i) MomTrunc::MCmeanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                               Sigma = as.matrix(sigma__), nu = as.matrix(round(params$params[j,"nu"])+1), dist = 't')$mean
        , which(phi1), 1:args$m
      )

      if(sum(phi1) == 0) w0 = 0

      aij = params$params[j,"lambda"]*(y[!phi1] - medias[!phi1,j])/params$params[j,"sigma"]

      auxU = (2*MT*params$params[j,"nu"]**(params$params[j,"nu"]/2)*gamma((params$params[j,"nu"]+2)/2)*
                (params$params[j,"nu"]+d2ij[!phi1]+(aij**2))**(-(params$params[j,"nu"]+2)/2))/
        (pi*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2)*pU)

      cv = 2*gamma((params$params[j,"nu"]+1)/2)/(gamma((params$params[j,"nu"])/2)*sqrt(params$params[j,"nu"]*(1+(params$params[j,"lambda"]**2))*pi))

      e10[!phi1] = (Mu[!phi1]*e00[!phi1]) + auxU
      e10[phi1] = MuAux*(e01[phi1]-e00[phi1]*medias[phi1, j])+MT*cv*(R0_F0)

      e20[!phi1] = ((Mu[!phi1]**2)*e00[!phi1]) + M2T + Mu[!phi1]*auxU
      e20[phi1] = (MuAux**2)*(e02[phi1] - 2*e01[phi1]*medias[phi1, j]+ e00[phi1]*(medias[phi1, j]**2)) + MuAux*(w0[phi1]-medias[phi1,j])*MT*cv*(R0_F0)+M2T

      e11[!phi1] = e10[!phi1]*y[!phi1]
      e11[phi1] = MuAux*(e02[phi1]-e01[phi1]*medias[phi1, j])+
        MT*cv*(R0_F0)*w0[phi1]
      return(list(Z, e00, e01, e02, e10, e20, e11))
    })

  U = lapply(1:7, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e00", "e01", "e02", "e10", "e20", "e11"))
  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)

  return(U)
}
.S3method("etapaE", "MixCenST", etapaE.MixCenST)

etapaE.MoECenST = function(y, X, params, medias, args, ...){

  phi1 = (args$phi[,1] == 1)

  U_list = sapply(
    1:args$g,
    function(j){

      Z = e00 = e01 = e02 = w0 = tau_gama = e10 = e11 = e20 = numeric(args$n)

      F0 = sn::pst(args$c2, medias[phi1, j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])-
        sn::pst(args$c1, medias[phi1,j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])

      d2ij = dMahalanobis(y, medias[, j], params$params[j,"sigma"])

      Z[!phi1] = params$P[!phi1, j]*sn::dst(y[!phi1], medias[!phi1, j], params$params[j,"sigma"], params$params[j,"lambda"], nu = params$params[j,"nu"])
      Z[phi1] = params$P[phi1, j]*F0

      MuAux = (params$params[j,"delta"]/(params$params[j,"delta"]**2+params$params[j,"gama"]))
      Mu = MuAux*(y-medias[,j])
      M2T = params$params[j,"gama"]/(params$params[j,"delta"]**2+params$params[j,"gama"])
      MT = sqrt(M2T)

      aij = params$params[j,"lambda"]*(y[!phi1] - medias[!phi1,j])/params$params[j,"sigma"]

      e00Aux = 4*(params$params[j,"nu"]**(params$params[j,"nu"]/2))*gamma((params$params[j,"nu"]+3)/2)*
            ((params$params[j,"nu"]+d2ij[!phi1])**(-(params$params[j,"nu"]+3)/2))/
            (sqrt(pi)*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2))*
            pt(aij*sqrt((params$params[j,"nu"] + 3)/(params$params[j,"nu"] + d2ij[!phi1])), params$params[j,"nu"] + 3)

      pU = sn::dst(y[!phi1], medias[!phi1,j], params$params[j,"sigma"], params$params[j,"lambda"], params$params[j,"nu"])

      sigma_ = (params$params[j,"sigma"]**2)*(params$params[j,"nu"]/(params$params[j,"nu"]+2))
      P0 = sn::pst(args$c2, medias[phi1, j], sqrt(sigma_), params$params[j,"lambda"], nu = params$params[j,"nu"]+2)-
        sn::pst(args$c1, medias[phi1,j], sqrt(sigma_), params$params[j,"lambda"], nu = params$params[j,"nu"]+2)

      sigma__ = (params$params[j,"sigma"]**2)*(params$params[j,"nu"]/((params$params[j,"nu"]+1)*(1+params$params[j,"lambda"]**2)))

      R0 = sn::pst(args$c2, medias[phi1, j], sqrt(sigma__),  0, params$params[j,"nu"]+1)-
        sn::pst(args$c1, medias[phi1,j], sqrt(sigma__), 0, params$params[j,"nu"]+1)

      R0_F0 = R0/ifelse(F0 == 0, .Machine$double.xmin, F0)

      e00[!phi1] = e00Aux/ifelse(pU == 0, .Machine$double.xmin, pU)
      e00[phi1] = P0/ifelse(F0 == 0, .Machine$double.xmin, F0)

      moments = mapply(
        function(ic, i){
          if(params$params[j,"lambda"] == 0){
            return(
              MomTrunc::MCmeanvarTMD(
                lower = args$c1[i],
                upper = args$c2[i],
                mu = medias[ic, j],
                Sigma = as.matrix(sigma_),
                nu = as.matrix(round(params$params[j,"nu"])+2),
                dist = 't')[c(1,2)])
          }else{
            return(
              MomTrunc::MCmeanvarTMD(
                lower = args$c1[i],
                upper = args$c2[i],
                mu = medias[ic, j],
                Sigma = as.matrix(sigma_),
                lambda = as.matrix(params$params[j,"lambda"]),
                nu = as.matrix(round(params$params[j,"nu"])+2),
                dist = 'ST')[c(1,2)])
          }}
        , which(phi1), 1:args$m
      )

      if(length(moments) == 0) moments = rbind(0, 0)

      e01[!phi1] = e00[!phi1]*y[!phi1]
      e01[phi1] = e00[phi1]*unlist(moments[1,])

      e02[!phi1] = e00[!phi1]*(y[!phi1]**2)
      e02[phi1] = e00[phi1]*unlist(moments[2,])

      w0[phi1] = mapply(
            function(ic, i) MomTrunc::MCmeanvarTMD(lower = args$c1[i], upper = args$c2[i], mu = medias[ic, j],
                                             Sigma = as.matrix(sigma__), nu = as.matrix(round(params$params[j,"nu"])+1), dist = 't')$mean
            , which(phi1), 1:args$m
          )

      if(sum(phi1) == 0) w0 = 0

      aij = params$params[j,"lambda"]*(y[!phi1] - medias[!phi1,j])/params$params[j,"sigma"]

      auxU = (2*MT*params$params[j,"nu"]**(params$params[j,"nu"]/2)*gamma((params$params[j,"nu"]+2)/2)*
                   (params$params[j,"nu"]+d2ij[!phi1]+(aij**2))**(-(params$params[j,"nu"]+2)/2))/
            (pi*params$params[j,"sigma"]*gamma(params$params[j,"nu"]/2)*pU)

      cv = 2*gamma((params$params[j,"nu"]+1)/2)/(gamma((params$params[j,"nu"])/2)*sqrt(params$params[j,"nu"]*(1+(params$params[j,"lambda"]**2))*pi))

      e10[!phi1] = (Mu[!phi1]*e00[!phi1]) + auxU
      e10[phi1] = MuAux*(e01[phi1]-e00[phi1]*medias[phi1, j])+MT*cv*(R0_F0)

      e20[!phi1] = ((Mu[!phi1]**2)*e00[!phi1]) + M2T + Mu[!phi1]*auxU
      e20[phi1] = (MuAux**2)*(e02[phi1] - 2*e01[phi1]*medias[phi1, j]+ e00[phi1]*(medias[phi1, j]**2)) + MuAux*(w0[phi1]-medias[phi1,j])*MT*cv*(R0_F0)+M2T

      e11[!phi1] = e10[!phi1]*y[!phi1]
      e11[phi1] = MuAux*(e02[phi1]-e01[phi1]*medias[phi1, j])+
        MT*cv*(R0_F0)*w0[phi1]
      return(list(Z, e00, e01, e02, e10, e20, e11))
      })

  U = lapply(1:7, function(i) do.call(cbind, U_list[i,])) |>
    setNames(c("Z", "e00", "e01", "e02", "e10", "e20", "e11"))
  soma_z = apply(U$Z, 1, sum)
  U$Z = sapply(1:args$g, function(j) U$Z[,j]/soma_z)
  if(dim(U$Z)[1] != args$n) U$Z = t(U$Z)

  return(U)
}
.S3method("etapaE", "MoECenST", etapaE.MoECenST)
