library(openCyto)
library(flowDensity)


DrawGate <- function(fr, channels, adjust = 1, ...){
  
  # Check that length(channels) %in% c(1,2)
  if(!length(channels) %in% c(1:2)){
    stop("Please supply either 1 or 2 fluorescent channels for gating.")
  }
  
  # Determine whether R is being run in RStudio
  if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    # if TRUE we need to open X11() interactive graphics device
    X11()
  }
  
  if(length(channels) == 1){
    message("Select 2 points on the plot to define the limits of the 1D gate.")
    
    dens <- density(exprs(fr)[,channels], adjust = adjust)
    plot(dens, main = " ")
    
    pts <- locator(type = "o", lwd = 2, pch = 16)
    abline(v = pts$x, col = "red")
    
    pts <- matrix(pts$x, ncol = 1, dimnames = list(c("min","max"),channels))
  
  } else if(length(channels) == 2){
  # Plot the data for gating use flowDensity::plotDens - locator() only works for base graphics
  cat("Draw 2D polygon gate around population. \n")
  
  flowDensity::plotDens(fr, channels = channels, cex = 3, ...)
  
  # Extract points of drawn gate
  pts <- locator(type = "o", lwd = 2, pch = 16)
  
  if (length(pts$x) < 3) stop("A minimum of 3 points is required to construct a polygon gate.")
  lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
  
  pts <- as.data.frame(pts)
  colnames(pts) <- channels
  }
  
  return(pts)
}


gate_draw <- function(fr, pp_res, channels, filterId = "", gate_range = NULL, min = NULL, max = NULL, ...){
  
  # Two fluorescent channels must be supplied
  if(missing(channels) | !length(channels) %in% c(1:2)){
    stop("Two fluorescent channels must be specified to draw polygon gate.")
  }
  
  # Truncate flowFrame if min and max arguments are supplied
  if(!(is.null(min) && is.null(max))){
    fr <- openCyto::.truncate_flowframe(fr, channels = channels, min = min, max = max)
  }
  
  # Determine vertices of polygon using DrawGate
  pts <- DrawGate(fr, channels)
  
  if(length(channels) == 1){
    rectangleGate(.gate = pts)
  
  }else if (length(channels) == 2){
  
  # Construct polygonGate
    polygonGate(.gate = pts)
  }
}

#### Register DrawGate with openCyto
registerPlugins(fun = gate_draw, methodName = "DrawGate")

##### Load in Samples ----------------------
fs <- read.flowSet(path = "Samples", pattern = ".fcs")
gs <- GatingSet(fs)

trans <- estimateLogicle(gs[[1]], colnames(fs)[-c(1,14)])
gs <- transform(gs,trans)

template <- add_pop(
  gs, alias = "Cells", pop = "+", parent = "root", dims = "FSC-A", gating_method = "DrawGate",
  collapseDataForGating = TRUE, groupBy = 2, gating_args = "subSample=20000"
)
Rm("Cells",gs)
