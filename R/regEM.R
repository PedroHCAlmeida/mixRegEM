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

    if(ll[1] == ll[2] == ll[3]) crit = 0

    crit = abs(llInf - ll[3])

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
  nPar = length(Par[!is.na(Par)])-length(args$lambda)
  aic = -2*ll[3] + 2*nPar
  bic = -2*ll[3] + log(args$n)*nPar

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
    l = ll[3],
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


