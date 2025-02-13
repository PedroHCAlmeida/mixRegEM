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
search = function(..., g = 2:8, family = "MixNormal", criteria = "BIC", verbose = F, best = F, max_iters_search = 5){

  results = expand.grid(g, family)
  colnames(results) = c("family", "g")

  models = list()
  for(f in family){
    for(gi in g){
      print(sprintf("Start g = %s; family = %s", gi, f))
      cont = 0
      print(sprintf("Tentativas: %s", cont))
      x = list()
      while(length(x)==0 && cont < max_iters_search){
        x = tryCatch({
          x = do.call(function(...) regEM(..., g = gi, family = f, verbose = verbose), list(...))
          },
          error = function(e) return(list())
          )
        cont = cont+1
        print(x)
      }
      if(is.null(x)){print("Fail")}
      else{print(sprintf("End g = %s; family = %s", gi, f))}

      models[[f]][[gi]] = x
    }
  }

    # mapply(
    #   function(g, family){
    #     cont = 0
    #     x = NULL
    #     while(is.null(x) && cont < max_iters_search){
    #       x = try(do.call(function(...) regEM(..., g = g, family = family, verbose = verbose), list(...)))
    #       cont = cont+1
    #     }
    #     return(x)
    #     },
    #   g = as.integer(results$g),
    #   family = as.character(results$family),
    #   SIMPLIFY = F
    # )

  #names(models) = paste0(results$family, results$g)

  metrics = models |>
    sapply(
      function(x){
        x |>
          sapply(
            function(y){
               bic = try({y[[criteria]]}, silent = T)
               if(is.null(bic)) return(NA)
               return(bic)
            }
          )
      }
    )

  results[,criteria] = unlist(c(metrics))
  print(results)

  if(best) best_model = tryCatch({best_model = models[[floor(which.min(metrics)/max(g))]][[which.min(metrics)%%max(g)]]}, error = function(e) NULL, silent = T)
  else best_model = models
  list(best = best_model, results = results)
}

