#' @title search
#' @param ... Parâmetros da função regEM
#' @param g número de grupos
#' @param tol tolerância
#' @param family família
#' @param criteria critério de seleção
#' @param verbose mostrar a verossimilhança passo a passo
#' @param best retornar apenas o melhor modelo
#' @return resultados finais
#' @export
search = function(..., g = 2:8, family = "MixNormal", criteria = "BIC", verbose = F, best = F){

  results = expand.grid(g, family)
  colnames(results) = c("g", "family")

  models =
    mapply(
      function(g, family) try(do.call(function(...) regEM(..., g = g, family = family, verbose = verbose), list(...))),
      g = as.integer(results$g),
      family = as.character(results$family),
      SIMPLIFY = F
    )

  names(models) = paste0(results$family, results$g)

  metrics = models |>
    sapply(function(x) try({x[[criteria]]}, silent = T))

  results[,criteria] = metrics
  print(results)

  if(best) best_model = tryCatch({best_model = models[[which.min(metrics)]]}, error = function(e) NULL, silent = T)
  else best_model = models
  list(best = best_model, results = results)
}

