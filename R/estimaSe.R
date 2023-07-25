estimaSe = function(reg, y, X, params = NULL, args = NULL, weights = 1, ...){
  UseMethod("estimaSe")
}

estimaSe.Normal = function(reg, y, X, params = NULL, args = NULL, weights = 1, ...){

 if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)

    args = reg$args
  }

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

estimaSe.MixNormal = function(reg, y, X, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
  }

  lapply(1:args$g,
         function(j) estimaSe.Normal(X, params[j,], args = args, weights = U$Z[,j])
  )
}
.S3method("estimaSe", "MixNormal", estimaSe.MixNormal)

estimaSe.MoENormal = function(reg, y, X, params = NULL, args = NULL, U = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    U = reg$U

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
  }

  lapply(1:args$g,
         function(j) estimaSe.Normal(y, X, params[j,],
                                     args = args, weights = U$Z[,j])
  )
}
.S3method("estimaSe", "MoENormal", estimaSe.MoENormal)

estimaSe.MixT = function(reg, y, X, params = NULL, args = NULL, U = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    U = reg$U

    X = cbind(rep(1, nrow(X)), x)
    R = cbind(rep(1, args$n), as.matrix(x))
    args = reg$args
  }

  lapply(1:args$g,
         function(j) estimaSe.Normal(y, X, params[j,], args = args, weights = U$Z[,j]*U$K[,j])
  )
}
.S3method("estimaSe", "MixT", estimaSe.MixT)

estimaSe.MoET = function(reg, y, X, params = NULL, args = NULL, U = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    U = reg$U

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
  }

  lapply(1:args$g,
         function(j) estimaSe.Normal(y, X, params[j,], args = args, weights = U$Z[,j]*U$K[,j])
  )
}
.S3method("estimaSe", "MoET", estimaSe.MoET)

estimaSe.MoESN = function(reg, y, X, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    U = reg$U

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args

    R = cbind(rep(1, args$n), as.matrix(R))
    P = reg$P
  }

  dMoeSE = function(yi, Xi, beta, sigma, lambda, nu, P, args){
    sum(log(sapply(1:args$g, function(j) P[j]*sn::dsn(yi, Xi%*%beta[j,], sigma[j], lambda[j]))))
  }

  veroSE = function(theta, yi, Xi, Ri, args){

    beta = matrix(theta[1:(args$p*args$g)], ncol = args$p, byrow = F)
    sigma = theta[(args$p*args$g+1):((args$p*args$g)+args$g)]
    lambda = theta[(((args$p*args$g)+args$g)+1):((args$p*args$g)+2*args$g)]
    alpha = theta[((args$p*args$g)+2*args$g+1):length(theta)]

    P = matrizP(matrix(alpha, ncol = args$g-1, byrow = T), Ri)
    return(dMoeSE(yi, Xi, beta, sigma, lambda, nu, P, args))
  }

  theta = matrix(params[, -c(args$p+1, args$p+2)], ncol = args$g)
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]
  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))

  InfObs = lapply(
    1:args$n,
    function(i){
      score = numDeriv::grad(veroSE, theta_total[!is.na(theta_total)], yi = y[i], Xi = X[i,], Ri = R[i,], args = args)
      score %*% t(score)
    }
  )

  SE = sqrt(diag(solve(Reduce("+", InfObs))))
  names(SE) = unlist(sapply(colnames(theta), function(x) paste0(x, 1:ifelse(startsWith(x, "alpha"), args$g-1, args$g))))

  estat = theta_total[!is.na(theta_total)]/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoESN", estimaSe.MoESN)

estimaSe.MixCenSN = function(reg, y, X, params = NULL, args = NULL, U = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    U = reg$U
  }

  P = colMeans(U$Z)

  dMixCen = function(yi, Xi, beta, sigma, lambda, P, args){

    phi1 = (args$phi[,1] == 1)
    return(sapply(1:args$g,
                  function(j){
                    medias = (Xi) %*% beta[j,]
                    total = numeric(args$n)
                    total[!phi1] = P[j]*sn::dsn(x = y[!phi1], xi = medias[!phi1], omega = sigma[j], alpha = lambda[j])
                    total[phi1] = P[j]*(sn::psn(args$c2,  medias[phi1], sigma[j], lambda[j])-sn::psn(args$c1, medias[phi1], sigma[j],  lambda[j]))
                    total
                  }) |>
             rowSums())
    }

  veroSE = function(theta, yi, Xi, phi, args, P){

    theta = matrix(theta, byrow = T, nrow = args$g)

    beta = theta[,1:(args$p)]
    sigma = theta[, args$p+1]
    if(is.null(args$lambda)) lambda = theta[,args$p+2]
    else lambda = args$lambda

    return(
      sum(log(dMixCen(yi, Xi, beta, sigma, lambda, P, args)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))

  theta_final = theta_total[!is.na(theta_total)]

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, P = P),
    as.vector(theta_total[!is.na(theta_total)])
    )

  SE = sqrt(-diag(solve(H)))
  names(SE) = paste0(rep(rownames(theta), args$g),
                     sapply(1:args$g, function(x) rep(x, length(rownames(theta)))))

  estat = (theta_total[!is.na(theta_total)]-theta_total[!is.na(theta_total)])/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MixCenSN", estimaSe.MixCenSN)

estimaSe.MoECenSN = function(reg, y, X, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    U = reg$U

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args

    R = cbind(rep(1, args$n), as.matrix(R))
    P = reg$P
  }

  dMoECen = function(yi, Xi, beta, sigma, lambda, P, args){

    phi1 = (args$phi[,1] == 1)
    return(sapply(1:args$g,
                  function(j){
                    medias = (Xi) %*% beta[j,]
                    total = numeric(args$n)
                    total[!phi1] = P[!phi,j]*sn::dsn(x = y[!phi1], xi = medias[!phi1], omega = sigma[j], alpha = lambda[j])
                    total[phi1] = P[phi,j]*(sn::psn(args$c2,  medias[phi1], sigma[j], lambda[j])-sn::psn(args$c1, medias[phi1], sigma[j],  lambda[j]))
                    total
                  }) |>
             rowSums())
  }

  veroSE = function(theta, yi, Xi, phi, args, P){

    theta = matrix(theta, byrow = T, nrow = args$g)

    beta = theta[,1:(args$p)]
    sigma = theta[, args$p+1]
    alpha = theta[,(args$p+2):ncol(theta)]
    if(is.null(args$lambda)) lambda = theta[,args$p+2]
    else lambda = args$lambda

    P = matrizP(matrix(alpha, ncol = args$g-1, byrow = T), Ri)

    return(
      sum(log(dMoECen(yi, Xi, beta, sigma, lambda, P, args)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))
  theta_final = theta_total[!is.na(theta_total)]

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, P = P),
    as.vector(theta_total[!is.na(theta_total)])
  )

  SE = sqrt(-diag(solve(H)))
  names(SE) = paste0(rep(rownames(theta), args$g),
                     sapply(1:args$g, function(x) rep(x, length(rownames(theta)))))

  estat = theta_total[!is.na(theta_total)]/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoECenSN", estimaSe.MoECenSN)

estimaSe.MixCenST = function(reg, y, X, params = NULL, args = NULL, U = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    U = reg$U
  }
  P = colMeans(U$Z)

  dMixCen = function(yi, Xi, beta, sigma, lambda, P, args, nu){

    phi1 = (args$phi[,1] == 1)
    return(sapply(1:args$g,
                  function(j){
                    medias = (Xi) %*% beta[j,]
                    total = numeric(args$n)
                    total[!phi1] = P[j]*sn::dst(x = y[!phi1], xi = medias[!phi1], omega = sigma[j], alpha = lambda[j], nu = nu[j])
                    total[phi1] = P[j]*(sn::pst(args$c2,  medias[phi1], sigma[j], lambda[j], nu = nu[j])-sn::pst(args$c1, medias[phi1], sigma[j],  lambda[j], nu = nu[j]))
                    total
                  }) |>
             rowSums())
  }

  veroSE = function(theta, yi, Xi, phi, args, P){

    theta = matrix(theta, byrow = T, nrow = args$g)

    beta = theta[,1:(args$p)]
    sigma = theta[, args$p+1]
    if(is.null(args$lambda)) lambda = theta[,args$p+2]
    else lambda = args$lambda

    nu = theta[,ncol(theta)]

    return(
      sum(log(dMixCen(yi, Xi, beta, sigma, lambda, P, args, nu = nu)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))

  theta_final = theta_total[!is.na(theta_total)]

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, P = P),
    as.vector(theta_total[!is.na(theta_total)])
  )

  SE = sqrt(-diag(solve(H)))
  names(SE) = paste0(rep(rownames(theta), args$g),
                     sapply(1:args$g, function(x) rep(x, length(rownames(theta)))))

  estat = theta_total[!is.na(theta_total)]/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MixCenST", estimaSe.MixCenST)

estimaSe.MoECenST = function(reg, y, X, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    R = cbind(rep(1, args$n), as.matrix(R))
    P = reg$P
  }

  dMixCen = function(yi, Xi, beta, sigma, lambda, P, args, nu){

    phi1 = (args$phi[,1] == 1)
    return(sapply(1:args$g,
                  function(j){
                    medias = (Xi) %*% beta[j,]
                    total = numeric(args$n)
                    total[!phi1] = P[!phi1,j]*sn::dst(x = y[!phi1], xi = medias[!phi1], omega = sigma[j], alpha = lambda[j], nu = nu[j])
                    total[phi1] = P[phi1,j]*(sn::pst(args$c2,  medias[phi1], sigma[j], lambda[j], nu = nu[j])-sn::pst(args$c1, medias[phi1], sigma[j],  lambda[j], nu = nu[j]))
                    total
                  }) |>
             rowSums())
  }

  veroSE = function(theta, yi, Xi, phi, args, P){

    theta = matrix(theta, byrow = T, nrow = args$g)

    beta = theta[,1:(args$p)]
    sigma = theta[, args$p+1]
    alpha = theta[,(args$p+2):(ncol(theta)-1)]
    if(is.null(args$lambda)) lambda = theta[,args$p+2]
    else lambda = args$lambda

    P = matrizP(matrix(alpha, ncol = args$g-1, byrow = T), Ri)
    nu = theta[,ncol(theta)]

    return(
      sum(log(dMixCen(yi, Xi, beta, sigma, lambda, P, args, nu = nu)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))

  theta_final = theta_total[!is.na(theta_total)]

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, P = P),
    as.vector(theta_total[!is.na(theta_total)])
  )

  SE = sqrt(-diag(solve(H)))
  names(SE) = paste0(rep(rownames(theta), args$g),
                     sapply(1:args$g, function(x) rep(x, length(rownames(theta)))))

  estat = theta_total[!is.na(theta_total)]/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoECenST", estimaSe.MoECenST)
