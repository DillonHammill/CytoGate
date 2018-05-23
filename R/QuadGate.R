#' Quadrant Gates for Flow Cytometry Data
#' 
#' Simply select a single point designating the center of the 4 rectangleGates. Gate will be constructed in the
#' following order: bottom left, bottom right, top right and top left.
#'
#' \code{QuadGate} facilitates manual drawing of quadrant gates onto flow cytometry plots.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector of length 2 indicating the fluorescent channels to be used to construct 2D plot
#' and \code{QuadGates}.
#' @param ... additional arguments for plotDens.
#'
#' @return a \code{dataframe} object containing the quadrant gate filters.
#'
#' @keywords manual, gating, draw, FlowJo, polygonGate, openCyto
#' @import flowDensity
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
QuadGate <- function(fr, channels, ...){
  
  # Check that 2 channels have been supplied
  if(!length(channels) == 2){
    stop("Please supply 2 fluorescent channels for gating")
  }
  
  # Extract data for plotting and gating
  x <- exprs(fr[,channels])
  
  # Determine whether R is being run in RStudio
  if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    # if TRUE we need to open X11() interactive graphics device
    X11()
  }
  
  # Plot the data for gating use flowDensity::plotDens - locator() only works for base graphics
  cat("Select a single point designating the center of the QuadGates. \n")
  
  flowDensity::plotDens(fr, channels = channels, cex = 3, ...)
  
  # Extract points of drawn gate
  pts <- locator(type = "o", lwd = 2, pch = 16)
  
  if (length(pts$x) > 1) stop("Only a single point is required to construct the rectangleGates.")
  lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
  
  pts <- as.data.frame(pts)
  colnames(pts) <- channels
  
  return(pts)
}

#' DrawGate plugin for openCyto
#'
#' \code{DrawGate} allows the user to draw polygon gates directly onto plots of flow cytometry data.
#' Simply left click to gate and right click to close the gate.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for gating.
#' @param pp_res output of preprocessing function.
#' @param channels a vector of length 2 indicating the channels used to construct the 2D plot.
#' @param filterId gate name assigned by openCyto from the \code{gatingTemplate}.
#' @param gate_range range in which gate should be constructed (only needed for autogating functions)
#' @param min argument passed to \code{truncate_flowFrame} to restrict data to values > \code{min}.
#' @param max argument passed to \code{truncate_flowFrame} to restrict data to values < \code{max}.
#' @param ... additional arguments passsed to \code{DrawGate}.
#'
#' @return a \code{polygonGate} constructed from coordinates supplied by \code{DrawGate}.
#'
#' @keywords manual, gating, polygon, polygonGate
#' @import flowDensity openCyto
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
#' 
#' @examples{
#' library(openCyto)
#' 
#' registerPlugins(fun = gate_draw, methodName = "DrawGate")
#' listgtMethods()   # check plugin has been registered with openCyto
#' 
#' fs <- read.flowSet(path = "Samples", pattern = ".fcs") # load in .fcs files
#' 
#' gs <- GatingSet(fs) # add flowSet to GatingSet
#' 
#' template <- add_pop(
#' gs, alias = "Lymphocytes", pop = "+", parent = "root", dims = "FSC-A,SSC-A", gating_method = "DrawGate",
#' gating_args = "subSample=10000", collapseDataForGating = TRUE, groupBy = 2
#' )
#' 
#' # gating window will open to construct gate left click vertices on plot and close gate by right click and selecting "stop".
#' 
#' ggcyto(gs[[1]], subset = "root", aes(x = "FSC-A",y = "SSC-A")) + geom_hex(bins = 100) + geom_stats()
#' 
#' }
gate_quad <- function(fr, pp_res, channels, filterId = "", gate_range = NULL, min = NULL, max = NULL, ...){
  
  # Two fluorescent channels must be supplied
  if(missing(channels) | length(channels) != 2){
    stop("Two fluorescent channels must be specified to draw polygon gate.")
  }
  
  # Truncate flowFrame if min and max arguments are supplied
  if(!(is.null(min) && is.null(max))){
    fr <- openCyto::.truncate_flowframe(fr, channels = channels, min = min, max = max)
  }
  
  # Determine vertices of polygon using DrawGate
  pts <- QuadGate(fr, channels)
  
  # Construct 4 RectangleGates
    
    # Q1 <- Bottom Left
    q1.gate <- data.frame(x = c(-Inf,pts[1,1]), y = c(-Inf, pts[1,2]))
    colnames(q1.gate) <- channels
    q1 <- rectangleGate(.gate = q1.gate)
    
    # Q2 <- Bottom Right
    q2.gate <- data.frame(x = c(pts[1,1], Inf), y = c(-Inf, pts[1,2]))
    colnames(q2.gate) <- channels
    q2 <- rectangleGate(.gate = q2.gate)
    
    # Q3 <- Top Right
    q3.gate <- data.frame(x = c(pts[1,1], Inf), y = c(pts[1,2], Inf))
    colnames(q3.gate) <- channels
    q3 <- rectangleGate(.gate = q3.gate)
    
    # Q4 <- Top Left
    q4.gate <- data.frame(x = c(-Inf, pts[1,1]), y = c(pts[1,2], Inf))
    colnames(q4.gate) <- channels
    q4 <- rectangleGate(.gate = q4.gate)
    
    filters(list(q1,q2,q3,q4))
 
  
}
