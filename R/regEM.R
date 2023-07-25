#' @include classes.R
#' @include resultadosEM.R

#' @name regEM
#' @title regEM
#' @param y variável respostas
#' @param x variáveis explicativas
#' @param g número de grupos
#' @param tol tolerância
#' @param family família
#' @param grupoReal grupo real da classificação
#' @param max_iter número máximo de iterações
#' @param min_iter número mínimo de iterações
#' @param verbose mostrar a verossimilhança passo a passo
#' @return resultados finais
#' @export
regEM = function(y, x, g = 2, ..., tol = 1E-6, family = "MixNormal",
                 grupoReal = NULL, max_iter = 1000, min_iter = 1, verbose = F, showSE = F){

  if(length(g) > 1 | length(family) > 1)
    return(
      do.call(
        function(...) search(y, x, ..., g = g, tol = tol, family = family,
                             grupoReal = grupoReal, max_iter = max_iter,
                             min_iter = min_iter, verbose, showSE = showSE), list(...))
  )

  args = list(...)

  if(family %in% c("MixCenNormal", "MixCenT", "MoECenNormal", "MoECenT")){
    args$lambda = rep(0, g)
    family = switch(
      family,
      "MixCenNormal" = "MixCenSN",
      "MixCenT"= "MixCenST",
      "MoECenNormal" = "MoECenSN",
      "MoECenT" = "MoECenST"
    )
  }
  if(family %in% c("MixST")){
    args$phi = rep(0, nrow(x));args$c1 = -Inf;args$c2 = -Inf
    family = switch(
      family,
      "MixST" = "MixCenST"
    )
  }
  args$n = length(y)
  args$g = g
  args$verbose = verbose
  args$tol = tol
  if(!is.null(args$lambda) & length(args$lambda) == 1) args$lambda = rep(args$lambda, g)

  X = cbind(rep(1, args$n), x)
  gc(x, verbose = F)
  args$p = ncol(X)

  try({

    args$R = Matrix::Matrix(cbind(rep(1, args$n), args$r), sparse = T)
    if(identical(args$R, X) | is.null(args$r)){
      pointr::ptr("R", "X")
      args$R = R
    }
    args$k = ncol(args$R)
    gc(args$r, verbose = F)
  })

  try({
    args$m = sum(args$phi == 1)
    if(length(args$c1) != args$m){
      args$c1 = rep(args$c1, args$m)
    }
    if(length(args$c2) != args$m){
      args$c2 = rep(args$c2, args$m)
    }
    args$phi = Matrix::Matrix(args$phi, sparse = T)
  }, silent = T)

  y = eval(parse(text = family))(y)
  X = eval(parse(text = family))(X)

  paramsAtual = chuteInicial(y, X, args)
  medias = estimaMedia(X, paramsAtual$params, args)
  crit = 1
  it = 0
  ll = c(0,0,0)
  ll[3] = vero(y, medias, paramsAtual, args)
  ll[2] = ll[3]

  while(((crit > tol) & (it < max_iter)) | (it < min_iter)){

    ll = c(ll[-1], 0)

    # Etapa E
    U = etapaE(y, X, paramsAtual, medias, args)

    # Etapa M
    paramsNovo = etapaM(y, X, U, paramsAtual, args)

    # Estimando valores esperados
    medias = estimaMedia(X, paramsNovo$params, args)

    paramsAtual = paramsNovo

    # Calculando critério
    ll[3] = vero(y, medias, paramsAtual, args)
    #crit = abs((veroAtual-vero0)/(vero0))

    # Aitken Acceleration
    ck = (ll[3] - ll[2])/(ll[2] - ll[1])
    denom = max(1L - ck, .Machine$double.eps)
    llInf = ll[2]+(ll[3]-ll[2])/denom

    if((ll[1] == ll[2]) & (ll[2] == ll[3])){
      crit = 0
    } else{
      crit = abs(llInf - ll[3])
    }

    if(verbose){
      print(paramsNovo$params)
      cat('Loglikelihood =', ll[3], ' Critério:', crit, '\n')
    }
    it = it+1

    if(!is.finite(crit)) crit = .Machine$double.xmax
  }

  conv = T
  if((it == max_iter) && (it != min_iter)){
    if(verbose){
      print("Warning: O algortimo parou pelo máximo de itereções, e não convergiu")
    }
    conv = F
    }

  ll[3] = vero(y, medias, paramsAtual, args)

  Par = c(paramsAtual$params[,!colnames(paramsAtual$params) %in% c("delta", "gama")])
  nPar = length(Par[!is.na(Par)])-length(args$lambda)-length(args$nuFixo)

  if(is.null(args$nuFixo) & !is.null(args$nuIgual)){
    if(args$nuIgual == T) nPar = nPar-args$g+1
  }

  if(grepl("Mix", family)) nPar = nPar+args$g-1

  aic = -2*ll[3] + 2*nPar
  bic = -2*ll[3] + log(args$n)*nPar

  gruposEM = apply(U$Z, 1, which.max)

  props = prop.table(table(gruposEM))
  ordem = rank(-props)
  grupos_novo = ordem[gruposEM]

  if(all(gruposEM != grupos_novo)){
    paramsAtual$params = paramsAtual$params[ordem,]
    paramsAtual$P = matrix(paramsAtual$P, ncol = args$g)[,ordem]
    U = lapply(U, function(x) matrix(x, ncol = args$g)[,ordem])
  }

  rownames(paramsAtual$params) = 1:nrow(paramsAtual$params)
  if(showSE) se = estimaSe(y, X, paramsAtual, args = args, U = U) else se = NULL

  resultados = list(
    Iteracoes = it,
    Convergiu = conv,
    g = g,
    nPar = nPar,
    l = ll[3],
    AIC = aic,
    BIC = bic,
    Parametros = t(paramsAtual$params),
    U = U,
    se = se,
    P = paramsAtual$P,
    args = args
  )
  class(resultados) = c("resultadosEM", class(X))
  return(resultados)
}
#.S3method("regEM", "default", regEM.default)

# regEM.formula = function(form, data, g = 2, ..., tol = 1E-6, family = "MixNormal",
#                          grupoReal = NULL, max_iter = 1000, min_iter = 1, verbose = F, showSE = F){
#
#   v1 = rlang::f_lhs(form)
#   v2 = rlang::f_rhs(form)
#
#   y = data[,v1] |>
#     as.matrix()
#
#   x = data[,v2] |>
#     as.matrix()
#
#   print(x)
#
#   f = function(...)
#     regEM.default(y = y, x = x, g = g, ..., tol = tol, family = family,
#                   grupoReal = grupoReal, max_iter = max_iter, min_iter = min_ter,
#                   verbose = verbose, showSE = showSE)
#
#   do.call(f, list(...))
# }
# .S3method("regEM", "formula", regEM.formula)

