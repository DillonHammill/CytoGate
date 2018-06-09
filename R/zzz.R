#' Automatically register DrawGate and ManualGate with openCyto upon loading...
.onLoad <- function(libname,pkgname){
  openCyto::registerPlugins(fun = gate_draw, methodName = "DrawGate")
  openCyto::registerPlugins(fun = gate_manual, methodName = "ManualGate")
}
