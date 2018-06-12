#' Updated version of DrawGate
#'
#' \code{DrawGate} implements an interactive manual gating routine for flow cytometry data. Users can easily
#' select gate coordinates on plots of flow cytometry data using a mouse click. Based on the user input, \code{DrawGate}
#' can construct many different types of gates, including \code{polygon}, \code{rectangle}, \code{interval}, \code{threshold}, \code{ellipse},
#' and \code{quadrant}. The type of gate to be constructed must be supplied as the \code{gate_type} argument
#' which by default is set to a \code{polygonGate}. Each \code{gate_type} has specific gating instructions which are printed to
#' the console during gating. The selection of multiple gates is supported for rectangle and interval gate types.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param gate_type a character string of length 1 indicating the type of gate to be constructed. Supported gates are \code{"polygon"},
#' \code{"rectangle"}, \code{"interval"}, \code{"threshold"}, \code{"ellipse"} and \code{"quadrant"}.
#' @param N an integer indicating the number of gates to construct.
#' @param axis indicates which axis should be gated for \code{gate_type="interval"} with 2 fluorescent channel supplied.
#' @param adjust numeric smoothing factor used for 1D density plots.
#' @param ... additional arguments for plotDens.
#'
#' @return a \code{dataframe} object containing the coordinates required to construct the gate.
#'
#' @keywords manual, gating, draw, FlowJo, polygonGate, openCyto
#' @import flowDensity
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
cytoGate <- function(fr, channels, alias, subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  # Supported gate types
  gate_types <- c("polygon", "Polygon", "p", "P","rectangle", "Rectangle", "r", "R","interval", "Interval", "i", "I","threshold", "Threshold", "t", "T", "boundary", "Boundary", "b", "B","ellipse", "Ellipse", "e", "E","quadrant", "Quadrant", "q", "Q")
  
  # Check that gate_type has been supplied or default to polygon type
  if(missing(gate_type)){
    
    message("No gate type supplied - gate type set to polygon.")
    gate_type <- "polygon"
    
  }else if(anyNA(match(gate_type, gate_types))){
    
    message("Invalid gate type supplied - gate type set to polygon. Supported gate types include polygon, rectangle, interval, threshold, boundary, ellipse and quadrant")
    gate_type[which(is.na(match(gate_type,gate_types)))] <- "polygon"
    
  }
  
  gates <- mapply(function(gate_type, alias){
    
    if(dev.cur() == 2){ # Graphics device not open
      plot <- TRUE
    }else{
      plot <- FALSE
    }
    
    if(gate_type %in% c("polygon", "Polygon", "p", "P")){
      
      drawPolygon(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
    
    }else if(gate_type %in% c("rectangle", "Rectangle", "r", "R")){
      
      drawRectangle(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("interval", "Interval", "i", "I")){
      
      drawInterval(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, axis = axis, labs = labs,...)
      
    }else if(gate_type %in% c("threshold", "Threshold", "t", "T")){
      
      drawThreshold(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("boundary", "Boundary", "b", "B")){ 
      
      drawBoundary(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("ellipse", "Ellipse", "e", "E")){
      
      drawEllipse(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }else if(gate_type %in% c("quadrant", "Quadrant", "q", "Q")){
      
      drawQuadrants(fr = fr, channels = channels, alias = alias, subSample = subSample, plot = plot, labs = labs,...)
      
    }
  }, gate_type, alias)
  
  names(gates) <- alias
  return(gates)
  
}
