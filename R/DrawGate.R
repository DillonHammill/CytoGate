#' Interactive Gate Drawing for Flow Cytometry Data.
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
DrawGate <- function(fr, channels, gate_type, N = 1, axis = "x", adjust = 1.5,...){
  
  # Check that length(channels) %in% c(1,2)
  if(!length(channels) %in% c(1,2) | missing(channels)){
    stop("Please supply fluorescent channel(s) for gating.")
  }
  
  # Determine whether R is being run in RStudio
  if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    # if TRUE we need to open X11() interactive graphics device
    X11()
  }
  
  # Check that gate_type has been supplied or default to polygon type
  if(missing(gate_type)){
    
    message("No gate type supplied - gate type set to polygon.")
    gate_type <- "polygon"
    
  }else if(length(gate_type) == 1 & !gate_type %in% c("polygon","rectangle", "interval", "threshold", "ellipse", "quadrant")){
    
    message("Invalid gate type supplied - gate type set to polygon. Supported gate types include polygon, rectangle, interval, threshold and quadrant")
    gate_type <- "polygon"
    
  }
  
  # Construct Plot
  if(length(channels) == 1){
    
    # Density Plot
    d <- density(exprs(fr)[,channels], adjust = adjust)
    plot(d, main=paste(channels))
    polygon(d, col="blue", border="black")
    
  }else if(length(channels) == 2){
    
    # 2D Plot using flowDensity::plotDens
    flowDensity::plotDens(fr, channels = channels, cex = 3, ...)
  }
  
  
  if(gate_type == "polygon"){
    
    if(length(channels) != 2) stop("Two fluorescent channels are required to construct polygonGate.")
    
    # Co-ordinates of gate(s)
    pts <- list()
    for(i in 1:N){
      
      cat("Select at least 3 points to construct a polygon gate around the population. \n")
      
      # Extract gate coordinates
      coords <- locator(type = "o", lwd = 2, pch = 16)
      
      if (length(coords$x) < 3) stop("A minimum of 3 points is required to construct a polygon gate.")
      lines(x = coords$x[c(1, length(coords$x))], y = coords$y[c(1, length(coords$x))], lwd = 2)
      
      coords <- as.data.frame(coords)
      colnames(coords) <- channels
      pts[[i]] <- coords
    }
    
    gates <- lapply(pts, function(pts){
      pts <- data.frame(pts)
      colnames(pts) <- channels
      polygonGate(.gate = pts)
    })
    
    gates <- filters(gates)
    
  }else if(gate_type == "rectangle"){
    
    if(length(channels) != 2) stop("Two fluorescent channels are required to construct rectangleGate.")
    
    # Coordinates of gate(s)
    pts <- list()
    for(i in 1:N){
      
      cat("Select 2 diagonal points to construct a rectangle gate around the population. \n")
      
      # Extract gate coordinates
      coords <- locator(n = 2, type = "p", lwd = 2, pch = 16, col = "red")
      coords <- data.frame(coords)
      colnames(coords) <- channels
      pts[[i]] <- coords
      
      rect(xleft = min(coords[,1]), ybottom = min(coords[,2]), xright = max(coords[,1]), ytop = max(coords[,2]), border = "red", lwd = 2)
    }
    
    gates <- lapply(pts, function(pts){
      pts <- data.frame(pts)
      colnames(pts) <- channels
      rectangleGate(pts)
    })
    
    gates <- filters(gates)
    
  }else if(gate_type == "interval"){
    
    # Coordinates of gates(s)
    pts <- list()
    for(i in 1:N){
      
      cat("Select 2 points to define the lower and upper bounds of the population. \n")
      
      # Extract gate coordinates
      coords <- locator(n=2, type = "o", lwd = 2, pch = 16, col = "red")
      coords <- data.frame(coords)
      
      if(length(channels) == 1){
        colnames(coords) <- c(channels[1],"Density")
      }else{
        colnames(coords) <- channels
      }
      pts[[i]] <- coords
      
      if(axis == "x"){
        abline(v = coords[,1], lwd = 2, col = "red")
      }else if(axis == "y"){
        abline(h = coords[,2], lwd = 2, col = "red")
      }
      
    }
    
    if(axis == "x"){
      gates <- lapply(pts, function(pts){
        
        if(length(channels) == 1){
          pts <- data.frame(x = pts[,1])
          colnames(pts) <- channels[1]
        }else if(length(channels) == 2){
          pts <- data.frame(x = pts[,1], y = c(-Inf,Inf))
          colnames(pts) <- channels
        }
        
        rectangleGate(.gate = pts)
      })
    }else if(axis == "y"){
      if(length(channels) == 1) stop("Cannot gate y axis if a single fluorescent channel is provided.")
      gates <- lapply(pts, function(pts){
        pts <- data.frame(x = c(-Inf,Inf), y = pts[,2])
        colnames(pts) <- channels
        
        rectangleGate(.gate = pts)
      })
    }
    
    gates <- filters(gates)
    
  }else if(gate_type == "threshold"){
    
    cat("Select 1 point indicating the lower bound of the threshold gate. \n")
    
    if(N > 1){
      message("Multiple threhold gates are not supported - a single threshold will be returned")
    }
    
    # Extract gate coordinates
    coords <- locator(n=1, type = "o", lwd = 2, pch = 16, col = "red")
    
    if(length(channels) == 1){
      pts <- data.frame(x = c(coords$x,Inf))
      colnames(pts) <- channels[1]
      abline(v = coords$x, lwd = 2, col = "red")
    }else if(length(channels) == 2){
      pts <- data.frame(x = c(coords$x,Inf), y = c(coords$y,Inf))
      colnames(pts) <- channels
      rect(xleft = min(coords$x), ybottom = min(coords$y), xright = max(exprs(fr)[,channels[1]]), ytop = max(exprs(fr)[, channels[2]]), border = "red", lwd = 2)
    }
    
    gates <- rectangleGate(.gate = pts)
    
  }else if(gate_type == "ellipse"){
    
    if(length(channels) != 2) stop("Two fluorescent channels are required to construct ellipsoidGate.")
    
    gates <- list()
    for(i in 1:N){
      
      cat("Select at 4 points to define limits of the ellipsoidGate. \n")
      
      # Extract gate coordinates
      coords <- locator(n=4, type = "p", lwd = 2, pch = 16, col = "red")
      coords <- data.frame(coords)
      
      # Find which points are on major axis
      dst <- as.matrix(stats::dist(coords))
      mj.pts <- coords[which(dst == max(dst), arr.ind = TRUE)[1,],]
      
      # Find which points are on minor axis
      mr.pts <- coords[!coords$x %in% mj.pts$x & !coords$y %in% mj.pts$y,]
      
      # Find center of the major axis
      center <- c((sum(mj.pts$x)/nrow(mj.pts)), (sum(mj.pts$y)/nrow(mj.pts)))
      points(x = center[1], y = center[2], col = "red", pch = 16)
      
      # Find major point which lies above center
      max.pt <- mj.pts[mj.pts$y > center[2] ,]
      
      # Radius of the major axis
      a <- stats::dist(mj.pts)/2
      
      # Radius of the minor axis
      b <- stats::dist(mr.pts)/2
      
      # Angle between horizontal line through center and max.pt
      if(max.pt[1] > center[1]){           # angle < pi/2
        mj.pt.ct <- cbind(max.pt[1],center[2])
        colnames(mj.pt.ct) <- c("x","y")
        adj <- stats::dist(rbind(center,mj.pt.ct))
        angle <- acos(adj/a)
      }else if(max.pt[1] < center[1]){     # angle > pi/2
        mj.pt.ct <- cbind(center[1], max.pt[2])
        colnames(mj.pt.ct) <- c("x","y")
        opp <- stats::dist(as.matrix(rbind(max.pt,mj.pt.ct)))
        angle <- pi/2 + asin(opp/a)
      }
      
      # Covariance matrix
      cinv <- matrix(c(0,0,0,0), nrow = 2, ncol = 2)
      cinv[1,1] <- (((cos(angle)*cos(angle))/(a*a)) + ((sin(angle)*sin(angle))/(b*b)))
      cinv[2,1] <- sin(angle)*cos(angle)*((1/(a*a))-(1/(b*b)))
      cinv[1,2] <- cinv[2,1]
      cinv[2,2] <- (((sin(angle)*sin(angle))/(a*a)) + ((cos(angle)*cos(angle))/(b*b)))
      
      cvm <- solve(cinv)
      
      dimnames(cvm) <- list(channels,channels)
      
      DescTools::DrawEllipse(x = center[1], y = center[2], radius.x = a, radius.y = b, rot = angle, border = "red", lwd = 2)
      
      gates[[i]] <- ellipsoidGate(.gate = cvm, mean = center)
      
    }
    
    gates <- filters(gates)
    
  }else if(gate_type == "quadrant"){
    
    if(length(channels) != 2) stop("Two fluorescent channels are required to construct quadrant gates.")
    
    if(length(N > 1)){
      message("For quadrant gates only N = 1 - ignoring argument N > 1")
    }
    
    # Construct quadrant gates
    
    cat("Select a single point designating the center of the quadrant gates. \n")
    
    # Extract points of drawn gate
    pts <- locator(n=1, type = "o", lwd = 2, pch = 16)
    
    lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
    abline(v = pts$x, h = pts$y, lwd = 2)
    
    pts <- as.data.frame(pts)
    colnames(pts) <- channels
    
    # Construct quadrant gates
    
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
    
    gates <- filters(list(q1,q2,q3,q4))
    
  }
  
  return(gates)
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
#' @param gate_range range in which gate should be constructed (only needed for autogating functions).
#' @param gate_type type of gate to be constructed, supported types include 
#' \code{c("polygon", "rectangle", "interval", "threshold", "ellipse", "quadrant")}.
#' @param min argument passed to \code{truncate_flowFrame} to restrict data to values > \code{min}.
#' @param max argument passed to \code{truncate_flowFrame} to restrict data to values < \code{max}.
#' @param N an integer indicating the number of gates to construct, set to 1 by default.
#' @param axis indicates the axis to use for gating for \code{gate_type="interval"} when 2 fluorescent channel are supplied.
#' @param adjust numeric smoothing factor used for 1D density plots.
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
#' gs, alias = "Lymphocytes", pop = "+", parent = "root", dims = "FSC-A,SSC-A", gating_method = "DrawGate", gaating_args = "gate_type='ellipse'",
#' gating_args = "subSample=10000", collapseDataForGating = TRUE, groupBy = 2
#' )
#' 
#' # gating window will open to construct gate left click vertices on plot and close gate by right click and selecting "stop".
#' 
#' ggcyto(gs[[1]], subset = "root", aes(x = "FSC-A",y = "SSC-A")) + geom_hex(bins = 100) + geom_stats()
#' 
#' }
gate_draw <- function(fr, pp_res, channels, filterId = "", gate_range = NULL, min = NULL, max = NULL, gate_type = c("polygon", "rectangle", "interval", "threshold", "ellipse", "quadrant"), N = 1, axis = c("x","y"), adjust = 1.5,...){
  
  gate_type <- match.arg(gate_type)
  
  if(missing(axis)){
    axis <- "x"
  }
  
  axis <- match.arg(axis)
  
  # Two fluorescent channels must be supplied
  if(missing(channels) | !length(channels) %in% c(1,2)){
    stop("Two fluorescent channels must be specified to draw polygon gate.")
  }
  
  # Truncate flowFrame if min and max arguments are supplied
  if(!(is.null(min) && is.null(max))){
    fr <- openCyto::.truncate_flowframe(fr, channels = channels, min = min, max = max)
  }
  
  # Determine vertices of polygon using DrawGate
  gates <- DrawGate(fr, channels, gate_type = gate_type, N = N, axis = axis, adjust = adjust)
  
  return(gates)
}
