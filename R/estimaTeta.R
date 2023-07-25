estimaTeta = function(y, X, ...){
  UseMethod("estimaTeta")
}

estimaTeta.Normal = function(y, X){

  Xl = t(X)
  beta = solve(Xl%*%X)%*%(Xl%*%y)

  sigma = sqrt(sum((y - (X%*%beta))^2)/(length(y) - ncol(X)))

  c(beta = as.matrix(beta), sigma = sigma)
}
.S3method("estimaTeta", "Normal", estimaTeta.Normal)

estimaTeta.MixNormal = function(y, X, Z){

  X = Matrix::Matrix(X, sparse = T)
  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*y))
  sigma = sqrt(sum(Z*(y - (X%*%beta))^2)/sum(Z))

  c(beta = as.matrix(beta), sigma = sigma)
}
.S3method("estimaTeta", "MixNormal", estimaTeta.MixNormal)

estimaTeta.MoENormal = function(y, X, Z, R, alpha, P){

  X = Matrix::Matrix(X, sparse = T)
  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*y))
  sigma = sqrt(sum(Z*(y - (X%*%beta))^2)/sum(Z))
  alphaNovo = alpha + 4*solve(t(as.matrix(R))%*%as.matrix(R))%*%(t(as.matrix(R))%*%(Z - P))

  list(params = c(beta = as.matrix(beta), sigma = sigma, alpha = alphaNovo))
}
.S3method("estimaTeta", "MoENormal", estimaTeta.MoENormal)

estimaTeta.MixT = function(y, X, Z, K){

  X = Matrix::Matrix(X, sparse = T)
  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z*K)%*%X)%*%Matrix::Matrix(Matrix::t(X)%*%Matrix::Matrix(Z*K*y))
  sigma = sqrt(sum(Z*K*((y - (X%*%beta))^2))/sum(Z))

  c(beta = as.matrix(beta), sigma = sigma)
}
.S3method("estimaTeta", "MixT", estimaTeta.MixT)

estimaTeta.MoET = function(y, X, Z, K, R, alpha, P){

  X = Matrix::Matrix(X, sparse = T)
  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z*K)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*K*y))
  sigma = sqrt(sum(Z*K*((y - (X%*%beta))^2))/sum(Z))
  alphaNovo = alpha + 4*solve(t(as.matrix(R))%*%as.matrix(R))%*%(t(as.matrix(R))%*%(Z - P))

  c(beta = as.matrix(beta), sigma = sigma, alpha = alphaNovo)
}
.S3method("estimaTeta", "MoET", estimaTeta.MoET)

estimaTeta.MixSN = function(y, X, medias, Z, t1, t2, deltaAtual){

  b = -sqrt(2/pi)

  X = Matrix::Matrix(X, sparse = T)

  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*(y-(deltaAtual*t1))))
  res = (y-X%*%beta)
  delta = sum(t1*res)/sum(t2)
  gama = sum(Z*(res**2)- (2*delta*res*t1) + (delta**2)*t2)/sum(Z)

  lambda = delta/sqrt(gama)
  sigma =sqrt(delta**2 + gama)

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda)
}
.S3method("estimaTeta", "MixSN", estimaTeta.MixSN)

estimaTeta.MixCenSN = function(y, X, Z, e01, e02, e10, e20, e11, delta, P, lambda, sigma){

  X = Matrix::Matrix(X, sparse = T)
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

  if(is.null(lambda)) lambda = delta/sqrt(gama)

  sigma = sqrt(delta**2 + gama)

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda)
}
.S3method("estimaTeta", "MixCenSN", estimaTeta.MixCenSN)

estimaTeta.MoECenSN = function(y, X, R, Z, e01, e02, e10, e20, e11, delta, alpha, P, lambda, sigma){

  X = Matrix::Matrix(X, sparse = T)
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
  alphaNovo = alpha + 4*solve(t(as.matrix(R))%*%as.matrix(R))%*%(t(as.matrix(R))%*%(Z - P))

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda, alpha = alphaNovo)
}
.S3method("estimaTeta", "MoECenSN", estimaTeta.MoECenSN)


estimaTeta.MixCenST = function(y, X, R, Z, e00, e01, e02, e10, e20, e11, delta, P, lambda, sigma){

  X = Matrix::Matrix(X, sparse = T)
  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z*e00)%*%X)%*%(Matrix::t(X)%*%(Z*(e01-e10*delta)))
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
  if(is.null(lambda)) lambda = delta/sqrt(gama)
  sigma = sqrt(delta**2 + gama)

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda)
}
.S3method("estimaTeta", "MixCenST", estimaTeta.MixCenST)

estimaTeta.MoECenST = function(y, X, R, Z, e00, e01, e02, e10, e20, e11, delta, alpha, P, lambda, sigma){

  X = Matrix::Matrix(X, sparse = T)
  beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z*e00)%*%X)%*%(Matrix::t(X)%*%(Z*(e01-e10*delta)))
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

  alphaNovo = alpha + 4*solve(t(as.matrix(R))%*%as.matrix(R))%*%(t(as.matrix(R))%*%(Z - P))

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda, alpha = alphaNovo)
}
.S3method("estimaTeta", "MoECenST", estimaTeta.MoECenST)
