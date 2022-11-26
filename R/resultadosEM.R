print.resultadosEM = function(resultadosEM){
  print(resultadosEM$`Parâmetros`)
}
.S3method("print", "resultadosEM", print.resultadosEM)
