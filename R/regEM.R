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

  args = list(...)
  args$n = length(y)
  args$g = g
  args$verbose = verbose
  args$tol = tol

  X = cbind(rep(1, args$n), x)
  args$p = ncol(X)

  try({
    args$R = cbind(rep(1, args$n), args$r)
    args$k = ncol(args$R)
  })

  try({
    args$m = sum(args$phi == 1)
    if(length(args$c1) != args$m){
      args$c1 = rep(args$c1, args$m)
      args$c2 = rep(args$c2, args$m)
    }
    })

  y = eval(parse(text = family))(y)
  X = eval(parse(text = family))(X)

  paramsAtual = chuteInicial(y, X, args)
  medias = estimaMedia(X, paramsAtual$params, args)
  crit = 1
  it = 0
  veroAnt = 0
  vero0 = vero(y, medias, paramsAtual, args)

  while(((crit > tol) & (it < max_iter)) | (it < min_iter)){

    # Etapa E
    U = etapaE(y, X, paramsAtual, medias, args)

    # Etapa M
    paramsNovo = etapaM(y, X, U, paramsAtual, args)

    # Estimando valores esperados
    medias = estimaMedia(X, paramsNovo$params, args)

    paramsAtual = paramsNovo

    # Calculando critério
    veroAtual = vero(y, medias, paramsAtual, args)

    #crit = abs((veroAtual-vero0)/(vero0))

    # Aitken Acceleration
    ck = (veroAtual-vero0)/(vero0-veroAnt)
    veroInf = vero0+(veroAtual-vero0)/(1-ck)
    crit = abs(veroAtual-veroInf)

    veroAnt = vero0
    vero0 = veroAtual

    if(verbose){
      print(paramsNovo$params)
      cat('Loglikelihood =', veroAtual, '\n')
    }
    it = it+1
  }

  conv = T
  if((it == max_iter) && (it != min_iter) && verbose){
    print("Warning: O algortimo parou pelo máximo de itereções, e não convergiu")
    conv = F
  }

  veroAtual = vero(y, medias, paramsAtual, args)

  Par = c(paramsAtual$params[,!colnames(paramsAtual$params) %in% c("delta", "gama")])
  nPar = length(Par[!is.na(Par)])-length(args$lambda)
  aic = -2*veroAtual + 2*nPar
  bic = -2*veroAtual + log(args$n)*nPar

  gruposEM = apply(U$Z, 1, which.max)
  r = rank(-table(gruposEM))
  grupos_ordem = numeric(args$n)
  for(i in 1:length(r)){
    grupos_ordem[gruposEM == names(r[i])] = r[i]
  }
  gruposEM = grupos_ordem

  rownames(paramsAtual$params) = 1:nrow(paramsAtual$params)
  if(showSE) se = estimaSe(y, X, paramsAtual, args = args, U = U) else se = NULL

  resultados = list(
    Iteracoes = it,
    Convergiu = conv,
    g = g,
    l = veroAtual,
    AIC = aic,
    BIC = bic,
    Parametros = t(paramsAtual$params),
    U = U,
    se = se,
    P = paramsAtual$P
  )
  class(resultados) = c("resultadosEM", family)
  return(resultados)
}


