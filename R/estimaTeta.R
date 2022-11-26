estimaTeta = function(Y, X, ...){
  UseMethod("estimaTeta")
}

estimaTeta.Normal = function(Y, X, args){

  Xl = t(X)
  beta = solve(Xl%*%X)%*%(Xl%*%Y)

  sigma = sqrt(sum((Y - (X%*%beta))^2)/(args$n - args$p))

  c(beta = beta, sigma = sigma)
}
.S3method("estimaTeta", "Normal", estimaTeta.Normal)

estimaTeta.MixNormal = function(Y, X, Z){

  beta = solve(t(X)%*%diag(Z)%*%X)%*%(t(X)%*%(Z*y))
  sigma = sqrt(sum(Z*(Y - (X%*%beta))^2)/sum(Z))

  c(beta = beta, sigma = sigma)
}
.S3method("estimaTeta", "MixNormal", estimaTeta.MixNormal)

estimaTeta.MoENormal = function(Y, X, Z, R, alpha, P){

  beta = solve(t(X)%*%diag(Z)%*%X)%*%(t(X)%*%(Z*y))
  sigma = sqrt(sum(Z*(Y - (X%*%beta))^2)/sum(Z))
  alphaNovo = alpha + 4*solve(t(R)%*%R)%*%(t(R)%*%(Z - P))

  list(params = c(beta = beta, sigma = sigma, alpha = alphaNovo))
}
.S3method("estimaTeta", "MoENormal", estimaTeta.MoENormal)

