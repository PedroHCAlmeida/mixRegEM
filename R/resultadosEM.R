#' @export
print.resultadosEM = function(resultadosEM){
  print(resultadosEM$`Parametros`)
}
.S3method("print", "resultadosEM", print.resultadosEM)

summary.resultadosEM = function(resultadosEM){

  se = resultadosEM$se
  params = resultadosEM$`Parametros`[!rownames(resultadosEM$`Parametros`) %in% c("delta", "gama"),]
  theta_na = do.call(rbind, lapply(1:nrow(params), function(i) cbind(params[i,])))
  theta = theta_na[!is.na(theta_na)]

  sig = ifelse(se[,2] > 0.05, ' ', ifelse(se[,2] > 0.01, '*', ifelse(se[2,] > 0.001, '**', '***')))
  tabela = data.frame(Estimate = theta, `Std error` = se[,1], `p-value` = paste0(se[,2], sig))

  cat('------------------------------------------------------------------------ \n')
  cat('  \n')
  print(tabela)
  cat('  \n')
  cat('Loglikelihood =', resultadosEM$l, '   AIC =', resultadosEM$AIC, '   BIC =', resultadosEM$BIC, '\n')
  cat('  \n')
  cat('  \n')
  cat('Ps: the results take into account the normality of the ML estimators for large n. \n')
  cat('------------------------------------------------------------------------')
}
.S3method("summary", "resultadosEM", summary.resultadosEM)
