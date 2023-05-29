estimaSe = function(f, y, X, params, args, weights = 1, ...){
  UseMethod("estimaSe")
}

estimaSe.Normal = function(f, y, X, params, args, weights = 1, ...){

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

estimaSe.MixNormal = function(f, y, X, params, args, weights = 1, ...){
  lapply(1:args$g,
         function(j) estimaSe.Normal(X, params$params[j,], args = args, weights = U$Z[,j])
  )
}
.S3method("estimaSe", "MixNormal", estimaSe.MixNormal)

estimaSe.MoENormal = function(f, y, X, params, args, U, weights = 1, ...){
  lapply(1:args$g,
         function(j) estimaSe.Normal(y, X, params$params[j,],
                                     args = args, weights = U$Z[,j])
  )
}
.S3method("estimaSe", "MoENormal", estimaSe.MoENormal)

estimaSe.MixT = function(f, y, X, params, args, U, weights = 1, ...){
  lapply(1:args$g,
         function(j) estimaSe.Normal(y, X, params$params[j,], args = args, weights = U$Z[,j]*U$K[,j])
  )
}
.S3method("estimaSe", "MixT", estimaSe.MixT)

estimaSe.MoET = function(f, y, X, params, args, U, weights = 1, ...){
  lapply(1:args$g,
         function(j) estimaSe.Normal(y, X, params$params[j,], args = args, weights = U$Z[,j]*U$K[,j])
  )
}
.S3method("estimaSe", "MoET", estimaSe.MoET)

estimaSe.MoECenSN = function(f, y, X, params, args, weights = 1, ...){

  dMoeSEUnc = function(f, yi, Xi, beta, sigma, lambda, nu, P, args){
    sum(sapply(1:args$g, function(j) P[j]*sn::dsn(yi, Xi%*%beta[j,], sigma[j], lambda[j])))
  }
  dMoeSECen = function(Xi, beta, sigma, lambda, nu, P, args){
    sum(sapply(1:args$g, function(j) P[j]*(sn::psn(args$c2,  Xi%*%beta[j,], sigma[j], lambda[j])-sn::pst(args$c1, Xi%*%beta[j,], sigma[j],  lambda[j]))))
  }

  veroSE = function(theta, yi, Xi, Ri, phi, args){

    beta = matrix(theta[1:(args$p*args$g)], ncol = args$p, byrow = F)
    sigma = theta[(args$p*args$g+1):((args$p*args$g)+args$g)]
    lambda = theta[(((args$p*args$g)+args$g)+1):((args$p*args$g)+2*args$g)]
    alpha = theta[((args$p*args$g)+2*args$g+1):length(theta)]

    P = matrizP(matrix(alpha, ncol = args$g-1, byrow = T), Ri)

    if(phi == 1){
      return(log(dMoeSECen(Xi, beta, sigma, lambda, nu, P, args)))
    }
    else{
      return(log(dMoeSEUnc(yi, Xi, beta, sigma, lambda, nu, P, args)))
    }
  }

  theta = params$params[, -c(args$p+1, args$p+2)]
  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))

  InfObs = lapply(
    1:args$n,
    function(i){
      score = grad(veroSE, theta_total[!is.na(theta_total)], yi = y[i], Xi = X[i,], Ri = args$R[i,], phi = args$phi[i], args = args)
      score %*% t(score)
    }
  )

  SE = sqrt(diag(solve(Reduce("+", InfObs))))
  names(SE) = unlist(sapply(colnames(theta), function(x) paste0(x, 1:ifelse(startsWith(x, "alpha"), args$g-1, args$g))))

  estat <- theta_total[!is.na(theta_total)]/SE
  valorP <- round(1 - pnorm(abs(estat)), 4)
  names(valorP) <- names(SE)

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoECenSN", estimaSe.MoECenSN)

estimaSe.MoECenST = function(f, y, X, params, args, weights = 1, ...){

  dMoeSEUnc = function(f, yi, Xi, beta, sigma, lambda, nu, P, args){
   sum(sapply(1:args$g, function(j) P[j]*sn::dst(yi, Xi%*%beta[j,], sigma[j], lambda[j], nu[j])))
  }
  dMoeSECen = function(Xi, beta, sigma, lambda, nu, P, args){
    sum(sapply(1:args$g, function(j) P[j]*(sn::pst(args$c2,  Xi%*%beta[j,], sigma[j], lambda[j], nu[j])-sn::pst(args$c1, Xi%*%beta[j,], sigma[j],  lambda[j], nu[j]))))
  }

  veroSE = function(theta, yi, Xi, Ri, phi, args){

    beta = matrix(theta[1:(args$p*args$g)], ncol = args$p, byrow = F)
    sigma = theta[(args$p*args$g+1):((args$p*args$g)+args$g)]
    lambda = theta[(((args$p*args$g)+args$g)+1):((args$p*args$g)+2*args$g)]
    alpha = theta[((args$p*args$g)+2*args$g+1):(length(theta)-(args$g))]
    nu = theta[(length(theta)-(args$g-1)):length(theta)]

    P = matrizP(matrix(alpha, ncol = args$g-1, byrow = T), Ri)

    if(phi == 1){
      return(log(dMoeSECen(Xi, beta, sigma, lambda, nu, P, args)))
    }
    else{
      return(log(dMoeSEUnc(yi, Xi, beta, sigma, lambda, nu, P, args)))
    }
  }

  theta = params$params[, -c(args$p+1, args$p+2)]
  theta_total = do.call(rbind, lapply(1:ncol(theta), function(i) cbind(theta[,i])))

  InfObs = lapply(
    1:args$n,
    function(i){
      score = grad(veroSE, theta_total[!is.na(theta_total)], yi = y[i], Xi = X[i,], Ri = args$R[i,], phi = args$phi[i], args = args)
      score %*% t(score)
    }
  )

  SE = sqrt(diag(solve(Reduce("+", InfObs))))
  names(SE) = unlist(sapply(colnames(theta), function(x) paste0(x, 1:ifelse(startsWith(x, "alpha"), args$g-1, args$g))))

  estat <- theta_total[!is.na(theta_total)]/SE
  valorP <- round(1 - pnorm(abs(estat)), 4)
  names(valorP) <- names(SE)

  se = cbind(SE, valorP)

  return(se)
}
.S3method("estimaSe", "MoECenST", estimaSe.MoECenST)



