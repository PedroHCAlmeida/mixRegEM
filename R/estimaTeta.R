estimaTeta = function(f, y, X, ...){
  UseMethod("estimaTeta")
}

estimaTeta.Normal = function(f, y, X){

  Xl = Matrix::t(X)
  beta = solve(Xl%*%X)%*%(Xl%*%y)

  sigma = sqrt(sum((y - (X%*%beta))^2)/(length(y) - ncol(X)))

  c(beta = as.matrix(beta), sigma = sigma)
}
.S3method("estimaTeta", "Normal", estimaTeta.Normal)

estimaTeta.MixNormal = function(f, y, X, Z){

  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%(Z*y))
  sigma = sqrt(sum(Z*(y - (X%*%beta))^2)/sum(Z))

  c(beta = as.matrix(beta), sigma = as.matrix(sigma))
}
.S3method("estimaTeta", "MixNormal", estimaTeta.MixNormal)

estimaTeta.MoENormal = function(f, y, X, Z, R, alpha, P){

  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%(Z*y))
  sigma = sqrt(sum(Z*(y - (X%*%beta))^2)/sum(Z))
  alphaNovo = alpha + 4*solve(Matrix::t(R)%*%R)%*%(Matrix::t(R)%*%(Z - P))

  list(params = c(beta = as.matrix(beta), sigma = sigma, alpha = as.matrix(alphaNovo)))
}
.S3method("estimaTeta", "MoENormal", estimaTeta.MoENormal)

estimaTeta.MixT = function(f, y, X, Z, K){

  beta = solve(Matrix::t(X)%*%Matrix::diag(Z*K)%*%X)%*%(Matrix::t(X)%*%(Z*K*y))
  sigma = sqrt(sum(Z*K*((y - (X%*%beta))^2))/sum(Z))

  c(beta = as.matrix(beta), sigma = sigma)
}
.S3method("estimaTeta", "MixT", estimaTeta.MixT)

estimaTeta.MoET = function(f, y, X, Z, K, R, alpha, P){

  beta = solve(Matrix::t(X)%*%Matrix::diag(Z*K)%*%X)%*%(Matrix::t(X)%*%(Z*K*y))
  sigma = sqrt(sum(Z*K*((y - (X%*%beta))^2))/sum(Z))
  alphaNovo = alpha + 4*solve(Matrix::t(R)%*%R)%*%(Matrix::t(R)%*%(Z - P))

  c(beta = as.matrix(beta), sigma = sigma, alpha = as.matrix(alphaNovo))
}
.S3method("estimaTeta", "MoET", estimaTeta.MoET)

estimaTeta.MixSN = function(f, y, X, medias, Z, t1, t2, deltaAtual){

  b = -sqrt(2/pi)

  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%(Z*(y-(deltaAtual*t1))))
  res = (y-X%*%beta)
  delta = sum(t1*res)/sum(t2)
  gama = sum(Z*(res**2)- (2*delta*res*t1) + (delta**2)*t2)/sum(Z)

  lambda = delta/sqrt(gama)
  sigma =sqrt(delta**2 + gama)

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda, sigma)
}
.S3method("estimaTeta", "MixSN", estimaTeta.MixSN)


estimaTeta.MoECenSN = function(f, y, X, R, Z, e01, e02, e10, e20, e11, delta, alpha, P, lambda, sigma){

  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%(Z*(e01-e10*delta)))
  medias = X%*%beta

  if(is.null(lambda)){
    delta = sum(Z*(e11-e10*medias))/sum(Z*e20)
  }
  else{
    if(lambda == 0){
      delta = 0
    }
  }

  delta = sum(Z*(e11-e10*medias))/sum(Z*e20)

  gama = sum(Z*(e02-2*e01*medias+medias**2+(delta**2)*e20-2*delta*e11+2*delta*e10*medias))/sum(Z)

  if(gama <= 0){
    gama = .Machine$double.xmin
  }

  lambda = delta/sqrt(gama)

  sigma = sqrt(delta**2 + gama)
  alphaNovo = alpha + 4*solve(Matrix::t(R)%*%R)%*%(Matrix::t(R)%*%(Z - P))

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda, alpha = as.matrix(alphaNovo))
}
.S3method("estimaTeta", "MoECenSN", estimaTeta.MoECenSN)

estimaTeta.MoECenST = function(f, y, X, R, Z, e00, e01, e02, e10, e20, e11, delta, alpha, P, lambda, sigma){

  beta = solve(Matrix::t(X)%*%Matrix::diag(Z*e00)%*%X)%*%(Matrix::t(X)%*%(Z*(e01-e10*delta)))
  medias = X%*%beta

  if(is.null(lambda)){
    delta = sum(Z*(e11-e10*medias))/sum(Z*e20)
  }
  else{
    if(lambda == 0){
      delta = 0
    }
  }

  gama = sum(Z*(e02-2*e01*medias+e00*(medias**2)+(delta**2)*e20-2*delta*e11+2*delta*e10*medias))/sum(Z)
  lambda = delta/sqrt(gama)
  sigma = sqrt(delta**2 + gama)

  alphaNovo = alpha + 4*solve(Matrix::t(R)%*%R)%*%(Matrix::t(R)%*%(Z - P))

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda, alpha = as.matrix(alphaNovo))
}
.S3method("estimaTeta", "MoECenST", estimaTeta.MoECenST)




