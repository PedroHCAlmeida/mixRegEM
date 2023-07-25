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

setClass("MoET",
         representation(Y = "Y"))

MoET = function(x){
  class(x) = c("MoET", "MixT", class(x))
  return(x)
}

setClass("MixSN",
         representation(Y = "Y"))

MixSN = function(x){
  class(x) = c("MixSN", "MixNormal", class(x))
  return(x)
}

setClass("MoESN",
         representation(Y = "Y"))

MoESN = function(x){
  class(x) = c("MoESN", "MoECenSN", "MixSN", "MixNormal", class(x))
  return(x)
}

setClass("MixST",
         representation(Y = "Y"))

MixST = function(x){
  class(x) = c("MixST", "MixCenST", "MixT", class(x))
  return(x)
}


MoEST = function(x){
  class(x) = c("MixST", "MixCenST", "MixT", class(x))
  return(x)
}

setClass("MoEST",
         representation(Y = "Y"))

MoEST = function(x){
  class(x) = c("MoEST", "MoECenST", "MoET", class(x))
  return(x)
}

setClass("MixCenSN",
         representation(Y = "Y"))

MixCenSN = function(x){
  class(x) = c("MixCenSN", "MixSN", "MixNormal", class(x))
  return(x)
}

setClass("MoECenSN",
         representation(Y = "Y"))

MoECenSN = function(x){
  class(x) = c("MoECenSN", "MixSN", "MixNormal", class(x))
  return(x)
}

setClass("MixCenST",
         representation(Y = "Y"))

MixCenST = function(x){
  class(x) = c("MixCenST", "MixST", "MixT", class(x))
  return(x)
}

setClass("MixCenST",
         representation(Y = "Y"))

MoECenST = function(x){
  class(x) = c("MoECenST", "MixSN", "MixNormal", class(x))
  return(x)
}

setClass("resultadosEM",
         representation(resultados = "list"))

