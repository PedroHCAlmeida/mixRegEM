#' @title estimaSe
#' @param reg modelo estimado da classe regEM
#' @param y variável respostas
#' @param x variáveis explicativas
#' @return se
#' @export
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

estimaSe.MoENormal = function(reg, y, X, r, params = NULL, args = NULL, U = NULL, weights = 1, ...){

  # if(is.null(params)){
  #   params = reg$Parametros |> t()
  #
  #   U = reg$U
  #
  #   X = cbind(rep(1, nrow(X)), X)
  #   args = reg$args
  # }
  #
  # lapply(1:args$g,
  #        function(j) estimaSe.Normal(y, X, params[j,],
  #                                    args = args, weights = U$Z[,j])
  # )
  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    R = cbind(rep(1, args$n), as.matrix(r))
    P = reg$P
  }

  dMoE = function(y, X, beta, sigma, lambda, P, args, nu){
    sum(log(rowSums(sapply(1:args$g, function(j) P[,j]*dnorm(y, X%*%beta[j,], sigma[j])))))
  }

  veroSE = function(theta, yi, Xi, phi, args, R, naGrupo){

    if(args$varEqual == T && !is.null(args$varEqual)){
      theta = append(
        theta,
        rep(theta[((args$g*args$p)+1)], args$g-1),
        ((args$g*args$p)+1)
      )
    }

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    alpha = matrix(
      theta[((args$g*args$p)+1+args$g):((args$g*args$p)+args$g+args$k*(args$g-1))],
      nrow = args$g-1
    )

    P_aux = matrizP(t(alpha), R)
    P = P_aux
    P[,!1:args$g %in% naGrupo] = P_aux[,1:(args$g-1)]
    P[,naGrupo] = P_aux[,args$g]

    return(
      dMoE(yi, Xi, beta, sigma, lambda, args = args, P = P)
    )
  }

  theta = as.matrix(t(params))
  theta_total = c(t(theta))
  theta_final = theta_total[!is.na(theta_total)]
  if(args$varEqual == T && !is.null(args$varEqual)){
    theta_final = unique(theta_final)
  }

  naGrupo = which(params |> apply(1, function(x) sum(is.na(x))) > 0)

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, R = R, naGrupo = naGrupo),
    theta_final
  )

  #SE = theta_final
  SE = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE))) try({SE = HelpersMG::SEfromHessian(-H)})


  names(SE) = rownames(theta) |>
    sapply(function(x){
      if(startsWith(x, "beta")){
        return(paste0(rep(x, args$g), 1:args$g))
      }
      if(startsWith(x, "sigma")){
        if(args$varEqual == T && !is.null(args$varEqual)) return(x)
        return(paste0(rep(x, args$g), 1:args$g))
      }
      if(startsWith(x, "alpha")){
        return(paste0(rep(x, args$g-1), 1:(args$g-1)))
      }
      }) |>
    unlist()

  estat = theta_final/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
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

estimaSe.MoESN = function(reg, y, X, r, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    R = cbind(rep(1, args$n), as.matrix(r))
    P = reg$P
  }

  dMoE = function(y, X, beta, sigma, lambda, P, args, nu){
    sum(log(rowSums(sapply(1:args$g, function(j) P[,j]*sn::dsn(y, X%*%beta[j,], sigma[j], lambda[j])))))
  }

  veroSE = function(theta, yi, Xi, phi, args, R, naGrupo){

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    if(is.null(args$lambda)){
      lambda = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+2*args$g)]
      alpha = matrix(
        theta[((args$g*args$p)+1+2*args$g):((args$g*args$p)+2*args$g+args$k*(args$g-1))],
        nrow = args$g-1
      )
    } else{
      lambda = args$lambda
      alpha = matrix(
        theta[((args$g*args$p)+1+args$g):((args$g*args$p)+args$g+args$k*(args$g-1))],
        nrow = args$g-1
      )
    }
    P_aux = matrizP(t(alpha), R)
    P = P_aux
    P[,!1:args$g %in% naGrupo] = P_aux[,1:(args$g-1)]
    P[,naGrupo] = P_aux[,args$g]

    return(
      dMoE(yi, Xi, beta, sigma, lambda, args = args, P = P)
    )
  }


  if(any(colnames(params) %in% c("delta", "gama", "nu"))){
    theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama", "nu"))]))
  }else{
    theta = as.matrix(t(params))
  }
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = c(t(theta))
  theta_final = theta_total[!is.na(theta_total)]

  naGrupo = which(params |> apply(1, function(x) sum(is.na(x))) > 0)

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, R = R, naGrupo = naGrupo),
    theta_final
  )

  SE = theta_total
  SE[!is.na(theta_total)] = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE[!is.na(theta_total)]))) try({SE[!is.na(theta_total)] = HelpersMG::SEfromHessian(-H)})

  names(SE) = rownames(theta) |>
    sapply(function(x) paste0(rep(x, args$g), 1:args$g))

  estat = theta_total/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoESN", estimaSe.MoESN)

estimaSe.MoET = function(reg, y, X, params = NULL, args = NULL, U = NULL, weights = 1, ...){
  estimaSe.MoEST(reg, y, X, params, args, U, weights, ...)
}
.S3method("estimaSe", "MoET", estimaSe.MoET)

estimaSe.MoEST = function(reg, y, X, r, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    R = cbind(rep(1, args$n), as.matrix(r))
    P = reg$P
  }

  dMoE = function(y, X, beta, sigma, lambda, P, args, nu){
    sum(log(rowSums(sapply(1:args$g, function(j) P[,j]*sn::dst(y, X%*%beta[j,], sigma[j], lambda[j], nu = nu[j])))))
  }

  veroSE = function(theta, yi, Xi, phi, args, R, naGrupo, nu){

    if(args$varEqual == T && !is.null(args$varEqual)){
      theta = append(
        theta,
        rep(theta[((args$g*args$p)+1)], args$g-1),
        ((args$g*args$p)+1)
      )
    }

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    if(is.null(args$lambda)){
      lambda = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+2*args$g)]
      alpha = matrix(
        theta[((args$g*args$p)+1+2*args$g):((args$g*args$p)+2*args$g+args$k*(args$g-1))],
        nrow = args$g-1
      )
    } else{
      lambda = args$lambda
      alpha = matrix(
        theta[((args$g*args$p)+1+args$g):((args$g*args$p)+args$g+args$k*(args$g-1))],
        nrow = args$g-1
        )
    }
    P_aux = matrizP(t(alpha), R)
    P = P_aux
    P[,!1:args$g %in% naGrupo] = P_aux[,1:(args$g-1)]
    P[,naGrupo] = P_aux[,args$g]

    return(
      dMoE(yi, Xi, beta, sigma, lambda, args = args, P = P, nu = nu)
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama", "nu"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = c(t(theta))
  theta_final = theta_total[!is.na(theta_total)]

  if(args$varEqual == T && !is.null(args$varEqual)){
    theta_final = unique(theta_final)
  }

  naGrupo = which(params |> apply(1, function(x) sum(is.na(x))) > 0)

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, R = R, naGrupo = naGrupo, nu = params[,"nu"]),
    theta_final
  )

  #SE = theta_final
  SE = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE))) try({SE = HelpersMG::SEfromHessian(-H)})


  names(SE) = rownames(theta) |>
    sapply(function(x){
      if(startsWith(x, "beta")){
        return(paste0(rep(x, args$g), 1:args$g))
      }
      if(startsWith(x, "sigma")){
        if(args$varEqual == T && !is.null(args$varEqual)) return(x)
        return(paste0(rep(x, args$g), 1:args$g))
      }
      if(startsWith(x, "alpha")){
        return(paste0(rep(x, args$g-1), 1:(args$g-1)))
      }
    }) |>
    unlist()

  estat = theta_final/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoEST", estimaSe.MoEST)

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

  veroSE = function(theta, yi, Xi, phi, args){

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    if(is.null(args$lambda)){
      lambda = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+2*args$g)]
    } else{
      lambda = args$lambda
    }

    p = theta[(length(theta)-args$g+2):length(theta)]
    P = c(p, 1-sum(p))

    return(
      sum(log(dMixCen(yi, Xi, beta, sigma, lambda, P, args)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = c(t(theta))
  theta_final = c(theta_total[!is.na(theta_total)], P[1:(args$g-1)])

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args),
    as.vector(theta_final)
    )

  SE = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE))) try({SE = HelpersMG::SEfromHessian(-H)})

  names(SE) = c(rownames(theta) |>
      sapply(function(x) paste0(rep(x, args$g), 1:args$g)),
    paste0("p", 1:(args$g-1))
  )

  estat = (theta_total[!is.na(theta_total)]-theta_total[!is.na(theta_total)])/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MixCenSN", estimaSe.MixCenSN)

estimaSe.MoECenSN = function(reg, y, X, r, params = NULL, args = NULL, weights = 1, ...){

  if(is.null(params)){
    params = reg$Parametros |> t()

    U = reg$U

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args

    R = cbind(rep(1, args$n), as.matrix(r))
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

  veroSE = function(theta, yi, Xi, phi, args, R, naGrupo){

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    if(is.null(args$lambda)){
      lambda = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+2*args$g)]
      alpha = theta[((args$g*args$p)+1+2*args$g):((args$g*args$p)+2*args$g+args$k*(args$g-1))]
    } else{
      lambda = args$lambda
      alpha = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+args$g+args$k*(args$g-1))]
    }

    P_aux = matrizP(alpha, R)
    P = P_aux
    P[,!1:args$g %in% naGrupo] = P_aux[,1:(args$g-1)]
    P[,naGrupo] = P_aux[,args$g]

    return(
      sum(log(dMoECen(yi, Xi, beta, sigma, lambda, P, args)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = c(t(theta))
  theta_final = theta_total[!is.na(theta_total)]

  naGrupo = which(params |> apply(1, function(x) sum(is.na(x))) > 0)
  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, R = R, naGrupo = naGrupo),
    theta_final
  )

  SE = theta_total
  SE[!is.na(theta_total)] = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE))) try({SE[!is.na(theta_total)] = HelpersMG::SEfromHessian(-H)})


  SE = theta_total
  SE[!is.na(theta_total)] = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE[!is.na(theta_total)]))) try({SE[!is.na(theta_total)] = HelpersMG::SEfromHessian(-H)})

  names(SE) = rownames(theta) |>
    sapply(function(x) paste0(rep(x, args$g), 1:args$g))

  estat = theta_total/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

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

    medias = apply(beta, 1, function(j) Xi%*%j)
    phi1 = (args$phi[,1] == 1)
    return(sapply(1:args$g,
                  function(j){
                    total = numeric(args$n)
                    total[!phi1] = P[j]*sn::dst(x = y[!phi1], xi = medias[!phi1, j], omega = sigma[j], alpha = lambda[j], nu = nu[j])
                    total[phi1] = P[j]*(sn::pst(args$c2,  medias[phi1,j], sigma[j], lambda[j], nu[j])-sn::pst(args$c1, medias[phi1,j], sigma[j],  lambda[j], nu[j]))
                    total
                  }) |>
             rowSums())
  }

  veroSE = function(theta, yi, Xi, phi, args, P, nu){

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    if(is.null(args$lambda)){
      lambda = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+2*args$g)]
    } else{
      lambda = args$lambda
    }
    p = theta[(length(theta)-args$g+2):length(theta)]
    P = c(p, 1-sum(p))

    #nu = theta[(length(theta)-args$g+1):length(theta)]
    l = sum(log(dMixCen(yi = yi, Xi = Xi, beta = beta, sigma = sigma, lambda = lambda, P = P, args = args, nu = nu)))
    return(l)
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama", "nu"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = c(t(theta))
  theta_final = c(theta_total[!is.na(theta_total)], P[1:(args$g-1)])

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, P = P, nu = params[,"nu"]),
    as.vector(theta_final)
  )

  SE = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE))) try({SE = HelpersMG::SEfromHessian(-H)})

  # if(!is.null(args$nuIgual)){
  #   if(args$nuIgual){
  #     SE = c(SE, rep(SE[length(SE)], args$g-1))
  #     theta_final = c(theta_final, rep(theta_final[length(theta_final)], args$g-1))
  #   }
  # }

  names(SE) = c(rownames(theta) |>
    sapply(function(x) paste0(rep(x, args$g), 1:args$g)),
    paste0("p", 1:(args$g-1))
    )

  estat = theta_final/SE
  valorP = round(1 - pnorm(abs(estat)), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MixCenST", estimaSe.MixCenST)

estimaSe.MoECenST = function(reg, y, X, r, params = NULL, args = NULL, weights = 1, ...){


  if(is.null(params)){
    params = reg$Parametros |> t()

    X = cbind(rep(1, nrow(X)), X)
    args = reg$args
    R = cbind(rep(1, args$n), as.matrix(r))
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

  veroSE = function(theta, yi, Xi, phi, args, R, naGrupo, nu){

    beta =  matrix(theta[1:(args$g*args$p)], nrow = args$g)
    sigma =  theta[((args$g*args$p)+1):((args$g*args$p)+args$g)]

    if(is.null(args$lambda)){
      lambda = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+2*args$g)]
      alpha = theta[((args$g*args$p)+1+2*args$g):((args$g*args$p)+2*args$g+args$k*(args$g-1))]
    } else{
      lambda = args$lambda
      alpha = theta[((args$g*args$p)+1+args$g):((args$g*args$p)+args$g+args$k*(args$g-1))]
    }

    P_aux = matrizP(alpha |> matrix(nrow = args$k, byrow = T), R)
    P = P_aux
    P[,!1:args$g %in% naGrupo] = P_aux[,1:(args$g-1)]
    P[,naGrupo] = P_aux[,args$g]

    return(
      sum(log(dMixCen(yi, Xi, beta, sigma, lambda, args = args, P = P, nu = nu)))
    )
  }

  theta = as.matrix(t(params[,-which(colnames(params) %in% c("delta", "gama", "nu"))]))
  if(!is.null(args$lambda)) theta = theta[-which(rownames(theta) == "lambda"), ]

  theta_total = c(t(theta))
  theta_final = theta_total[!is.na(theta_total)]

  naGrupo = which(params |> apply(1, function(x) sum(is.na(x))) > 0)

  H = numDeriv::hessian(
    function(x) veroSE(x, yi = y, Xi = X, phi = args$phi, args = args, R = R, naGrupo = naGrupo, nu = params[,"nu"]),
    theta_final
  )

  SE = theta_total
  SE[!is.na(theta_total)] = sqrt(diag(solve(-H)))
  if(!all(is.finite(SE[!is.na(theta_total)]))) try({SE[!is.na(theta_total)] = HelpersMG::SEfromHessian(-H)})

  names(SE) = rownames(theta) |>
    sapply(function(x) paste0(rep(x, args$g), 1:args$g))

  valorP = round(1 - pnorm(abs(theta_final), mean = 0, sd = SE), 4)
  names(valorP) = names(SE)

  valorP[!grepl("beta|alpha|lambda", names(SE))] = NA

  se = cbind("est" = theta_total, SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoECenST", estimaSe.MoECenST)
