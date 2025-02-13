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

estimaTeta.MixNormal = function(y, X, Z, lasso){

  if(length(Z) == 1) Z = rep(1, nrow(X))
  X = Matrix::Matrix(X, sparse = T)

  if(is.null(lasso) || lasso == F){
    beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*y))
    penalty = 0
  }else{
    if(lasso){



      model_cv = glmnet::cv.glmnet(
        x = as.matrix(scale(X)[,-1]),
        y = as.vector(scale(y)),
        intercept = T,
        folds = 5,
        type.measure = "deviance",
        family = "gaussian",
        weights = Z,
        alpha = 1,
        tresh = 1E-10
        )

      penalty = model_cv$lambda.1se
      beta = coef(model_cv, s = "lambda.1se") |> as.vector()

      beta[beta!=0] = solve(Matrix::t(X[,beta!=0])%*%Matrix::Diagonal(x = Z)%*%X[,beta!=0])%*%(Matrix::t(X[,beta!=0])%*%Matrix::Matrix(Z*y))

    }else{
      beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*y))
      penalty = 0
    }
  }
  sigma = sqrt(sum(Z*(y - as.matrix(X%*%beta))^2)/sum(Z))

  c(beta = as.matrix(beta), sigma = sigma, penalty = penalty/(sigma**2))
}
.S3method("estimaTeta", "MixNormal", estimaTeta.MixNormal)

estimaTeta.MoENormal = function(y, X, Z, R, alpha, P, lasso, class){

  if(length(Z) == 1){
    Z = P = class = rep(1, nrow(X))
  }

  X = Matrix::Matrix(X, sparse = T)

  if(is.null(lasso) || lasso == F){
    beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*y))
    alphaNovo = alpha + 4*solve(t(as.matrix(R))%*%as.matrix(R))%*%(t(as.matrix(R))%*%(Z - P))
    penalty = 0
  }else{
    if(lasso){
      model_cv = glmnet::cv.glmnet(
        x = as.matrix(scale(X)[,-1]),
        y = as.vector(scale(y)),
        intercept = T,
        folds = 5,
        type.measure = "deviance",
        family = "gaussian",
        weights = Z,
        #lambda = seq((ncol(X)-1)/100, (ncol(X)-1)/10, length.out = 30),
        alpha = 1,
        tresh = 1E-10
      )

      penalty = model_cv$lambda.1se
      beta = coef(model_cv, s = "lambda.1se") |> as.vector()

      beta[beta!=0] = solve(Matrix::t(X[,beta!=0])%*%Matrix::Diagonal(x = Z)%*%X[,beta!=0])%*%(Matrix::t(X[,beta!=0])%*%Matrix::Matrix(Z*y))

      if(!any(is.na(alpha))){
        model_cv_alpha = glmnet::cv.glmnet(
          x = as.matrix(scale(R)[,-1]),
          y = as.numeric(class),
          intercept = T,
          folds = 5,
          type.measure = "deviance",
          family = "multinomial",
          alpha = 1,
          tresh = 1E-10
        )
        alphaNovo = coef(model_cv_alpha, s = "lambda.1se")
        alphaNovo = alphaNovo$`1` |> as.vector()
        alphaNovo[alphaNovo!=0] = alpha[alphaNovo!=0] + 4*solve(t(as.matrix(R[,alphaNovo!=0]))%*%as.matrix(R[,alphaNovo!=0]))%*%(t(as.matrix(R[,alphaNovo!=0]))%*%(Z - P))
      }else{
        alphaNovo = rep(NA, ncol(R))
      }

    }else{
      beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*y))
      penalty = 0
    }
  }

  sigma = sqrt(sum(Z*(y - as.vector(X%*%beta))^2)/sum(Z))

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

estimaTeta.MixSN = function(y, X, medias, Z, t1, t2, deltaAtual, lasso){

  X = Matrix::Matrix(X, sparse = T)

  if(is.null(lasso) || lasso == F){
    beta = solve(Matrix::t(X)%*%Matrix::Diagonal(x = Z)%*%X)%*%(Matrix::t(X)%*%Matrix::Matrix(Z*(y-(deltaAtual*t1))))
    penalty = 0
  }else{
    if(lasso){
      y_trans = Matrix::Matrix(y-(deltaAtual*t1))

      model = glmnet::cv.glmnet(x = as.matrix(scale(X)), y = as.vector(y_trans),
                                weights = Z, nfolds = 5, intercept=FALSE)
      beta = as.vector(coef(model))
      beta = beta[-1]
      penalty = model$lambda.min
    }
  }

  medias = (X%*%beta)
  res = y-as.numeric(medias)
  delta = sum(Z*t1*res)/sum(Z*t2)
  gama = sum(Z*((res**2)-(2*delta*res*t1) + (delta**2)*t2))/sum(Z)

  lambda = delta/sqrt(gama)
  sigma =sqrt(delta**2 + gama)

  c(beta = as.matrix(beta), delta = delta, gama = gama, sigma = sigma, lambda = lambda, penalty = penalty/gama)
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

  gama = sum(Z*(e02-2*e01*medias+medias**2+(delta**2)*e20-2*delta*e11+2*delta*e10*medias))/sum(Z)

  if(!is.finite(gama))

  if(gama <= 0 | !is.finite(gama)){
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

  gama = sum(Z*(e02-2*e01*medias+medias**2+(delta**2)*e20-2*delta*e11+2*delta*e10*medias))/sum(Z)

  if(gama <= 0 | !is.finite(gama)){
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
