setClassUnion("Y",   c("numeric",   "matrix", "array"))

setClass("Normal",
         representation(Y = "Y"))

setClass("MixNormal",
         representation(Y = "Y"),
         contains="Normal")

MixNormal = function(x){
  class(x) = c("MixNormal", "Normal", class(x))
  return(x)
}

setClass("MoENormal",
         representation(Y = "Y"),
         contains = c("MixNormal"))

MoENormal = function(x){
  class(x) = c("MoENormal", "MixNormal", class(x))
  return(x)
}

setClass("MoEKernelNormal",
         representation(Y = "Y"),
         contains = c("MixNormal"))

MoEKernelNormal = function(x){
  class(x) = c("MoEKernelNormal", "MoENormal", class(x))
  return(x)
}

setClass("MixT",
         representation(Y = "Y"))

MixT = function(x){
  class(x) = c("MixT", class(x))
  return(x)
}

setClass("resultadosEM",
         representation(resultados = "list"))

