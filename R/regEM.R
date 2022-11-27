#' @include classes.R

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
                 grupoReal = NULL, max_iter = 1000, min_iter = 5, verbose = T){

  args = list(...)
  args$n = length(y)
  args$g = g

  X = cbind(rep(1, args$n), x)
  args$p = ncol(X)
  try({
    args$R = cbind(rep(1, args$n), args$r)
  })

  y = eval(parse(text = family))(y)
  X = eval(parse(text = family))(X)

  paramsAtual = chuteInicial(y, X, args)
  medias = estimaMedia(X, paramsAtual$params, args)
  crit = 1
  it = 0

    while(((crit > tol) & (it < max_iter)) | (it < min_iter)){
      # Calculando Verossimilhança
      vero0 = vero(y, X, paramsAtual)
      print(vero0)
      # Etapa E
      U = etapaE(y, X, paramsAtual, medias, args)

      # Etapa M
      paramsNovo = etapaM(y, X, U, paramsAtual, args)

      # Estimando valores esperados
      medias = estimaMedia(X, paramsAtual$params, args)

      paramsAtual = paramsNovo

      # Calculando critério
      crit = abs(vero(y, X, paramsAtual)/vero0 - 1)
      if(verbose) print(crit)
      it = it+1
    }

    gruposEM = apply(U$Z, 1, which.max)

    if(!is.null(grupoReal)){

      ordemReg = order(table(gruposEM))
      ordemReal = order(table(grupoReal))

      tabela = table(grupoReal, replace(gruposEM, ordemReg, ordemReal),
                     dnn = list('Real', 'Modelo'))

      #tabela = table(grupoReal, gruposEM, dnn = list('Real', 'Modelo'))
      }
    else tabela = NA

    rownames(paramsNovo$params) = 1:nrow(paramsNovo$params)
    eps = estimaEp(X, paramsNovo, args = args, U = U)
    metricas = lapply(1:g, function(j) calculaMetricas(y[gruposEM == j],
                                                       medias[gruposEM == j, j]))

    resultados = list(
      Iteracoes = it,
      l = vero(y, X, paramsNovo),
      Parametros = t(paramsNovo$params),
      U = U,
      eps = eps,
      tabela = tabela,
      medias = medias,
      metricas = metricas,
      params = paramsNovo,
      metricasTotais = colSums(do.call(rbind,
                                       lapply(metricas,
                                              function(x) sapply(2:4, function(i) x$n*x[[i]]/n)))),
      residuos = lapply(1:g, function(j) medias[gruposEM == j, j] - y[gruposEM == j])
    )
  class(resultados) = c("resultadosEM", family)
  return(resultados)
}

