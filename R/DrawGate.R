#' Interactive Gate Drawing for Flow Cytometry Data.
#'
#' \code{DrawGate} implements an interactive manual gating routine for flow cytometry data. Users can easily
#' select gate coordinates on plots of flow cytometry data using a mouse click. Based on the user input, \code{DrawGate}
#' can construct many different types of gates, including \code{polygon}, \code{rectangle}, \code{interval}, \code{multiinterval}, \code{threshold},
#' and \code{quadrant}. The type of gate to be constructed must be supplied as the \code{gate_type} argument
#' which by default is set to a \code{polygonGate}. Each \code{gate_type} has specific gating instructions which are printed to
#' the console during gating.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector indicating the fluorescent channel(s) to be used for gating. If a single channel is supplied, a histogram of
#' of the kernel density will be constructed.
#' @param gate_type a character string of length 1 indicating the type of gate to be constructed. Supported gates are \code{"polygon"},
#' \code{"rectangle"}, \code{"interval"}, \code{multiinterval}, \code{"threshold"} and \code{"quadrant"}.
#' @param ... additional arguments for plotDens.
#'
#' @return a \code{dataframe} object containing the coordinates required to construct the gate.
#'
#' @keywords manual, gating, draw, FlowJo, polygonGate, openCyto
#' @import flowDensity
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
DrawGate <- function(fr, channels, gate_type, ...){
  
  # Check that length(channels) %in% c(1,2)
  if(!length(channels) %in% c(1,2) | missing(channels)){
    stop("Please supply fluorescent channel(s) for gating.")
  }

  # Extract data for plotting and gating
  x <- exprs(fr[,channels])

  # Determine whether R is being run in RStudio
  if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    # if TRUE we need to open X11() interactive graphics device
    X11()
  }

  # Check that gate_type has been supplied or default to polygon type
  if(missing(gate_type)){
    message("No gate type supplied - gate type set to polygon.")
    gate_type <- "polygon"
  }else if(length(gate_type) == 1 & !gate_type %in% c("polygon","rectangle", "interval", "multiinterval","threshold", "quadrant")){
    message("Invalid gate type supplied - gate type set to polygon. Supported gate types include polygon, rectangle, interval, multiinterval, threshold and quadrant")
    gate_type <- "polygon"
  }
  
  
  if(gate_type == "polygon"){
    
    # Construct polygon gate 
    
    cat("Select at least 3 points to construct a polygon gate around the population. \n")
    
    # Create plot for gating
    flowDensity::plotDens(fr, channels = channels, cex = 3, ...)
    
    # Extract gate coordinates
    pts <- locator(type = "o", lwd = 2, pch = 16)
    
    if (length(pts$x) < 3) stop("A minimum of 3 points is required to construct a polygon gate.")
    lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
    
    pts <- as.data.frame(pts)
    colnames(pts) <- channels
    
  }else if(gate_type == "rectangle"){
  
    # Construct rectangle gate

    cat("Select 4 points to construct a rectangle gate around the population. \n")
    
    # Create plot for gating
    flowDensity::plotDens(fr, channels = channels, cex = 3, ...)
    
    # Extract gate coordinates
    pts <- locator(type = "o", lwd = 2, pch = 16)
    
    if (!length(pts$x) == 4) stop("Exactly 4 points are required to construct a rectangle gate.")
    lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
    rect(xleft = min(pts$x), ybottom = min(pts$y), xright = max(pts$x), ytop = max(pts$y), border = "red", lwd = 2)
    
    pts <- data.frame(x = c(min(pts$x), max(pts$x)), y = c(min(pts$y), max(pts$y)))
    colnames(pts) <- channels
  
  }else if(gate_type == "interval"){
  
    # Construct interval gate 
    
    cat("Select 2 points indicating the lower and upper bounds of the interval gate. \n")
    
    if (length(channels) != 1) stop("A single fluorescent channel is required to construct an interval gate")
    
    d <- density(exprs(fr)[,channels])
    plot(d, main=paste(channels))
    polygon(d, col="red", border="black")
    
    # Extract gate coordinates
    pts <- locator(type = "o", lwd = 2, pch = 16)
    
    if (!length(pts$x) == 2) stop("Exactly 2 points are required to construct an interval gate.")
    lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
    abline(v = pts$x, lwd = 2)
    
    pts <- data.frame(x = c(pts$x))
    colnames(pts) <- channels
  
  }else if(gate_type == "multiinterval"){
   
     # Construct multiple interval gates
    
    cat("Select 2 points per gate to indicate lower and upper bounds of each interval gate. \n")
    
    if (length(channels) != 1) stop("A single fluorescent channel is required to construct an interval gates.")
    
    d <- density(exprs(fr)[,channels])
    plot(d, main=paste(channels))
    polygon(d, col="red", border="black")
    
    # Extract gate coordinates
    pts <- locator(type = "p", pch = 16)
    pts <- data.frame(x = pts$x, y = pts$y)
    
    ord <- order(pts$x)
    pts <- pts[ord,]
    
    channel <- channels
    colnames(pts) <- c(channel,"Density")
    
    if(!nrow(pts) %% 2 == 0) stop("2 points are required per gate to define lower and upper bounds.")
    
    sp.pts <- split(pts, rep(1:(nrow(pts)/2),each=2))
    
    for(i in 1:length(sp.pts)){
    lines(x = sp.pts[[i]][,1], y = sp.pts[[i]][,2], lwd = 2, col = "black")
    abline(v = sp.pts[[i]][,1], lwd = 2)
    }
    
    pts <- sp.pts
    print(pts)
    print(pts[[1]][,1])
    
  }else if(gate_type == "threshold"){
  
    # Construct threshold gate

    cat("Select 1 points indicating the lower bound of the threshold gate. \n")
    
    if (length(channels) != 1) stop("A single fluorescent channel is required to construct an threshold gate")
    
    d <- density(exprs(fr)[,channels])
    plot(d, main=paste(channels))
    polygon(d, col="red", border="black")
    
    # Extract gate coordinates
    pts <- locator(type = "o", lwd = 2, pch = 16)
    
    if (!length(pts$x) == 1) stop("Exactly 1 point is required to define the lower bound of the threshold gate.")
    lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
    abline(v = pts$x, lwd = 2)
    
    pts <- data.frame(x = c(pts$x,Inf))
    colnames(pts) <- channels
  
  }else if(gate_type == "quadrant"){
  
    # Construct quadrant gates

    cat("Select a single point designating the center of the quadrant gates. \n")
    
    flowDensity::plotDens(fr, channels = channels, cex = 3, ...)
    
    # Extract points of drawn gate
    pts <- locator(type = "o", lwd = 2, pch = 16)
    
    if (length(pts$x) > 1) stop("Only a single point is required to construct the quadrant gates.")
    lines(x = pts$x[c(1, length(pts$x))], y = pts$y[c(1, length(pts$x))], lwd = 2)
    
    pts <- as.data.frame(pts)
    colnames(pts) <- channels
  }
  
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
#' @param gate_range range in which gate should be constructed (only needed for autogating functions).
#' @param gate_type type of gate to be constructed, supported types include 
#' \code{c("polygon", "rectangle", "interval", "multiinterval","threshold", "quadrant")}.
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
gate_draw <- function(fr, pp_res, channels, filterId = "", gate_range = NULL, min = NULL, max = NULL, gate_type = c("polygon", "rectangle", "interval", "multiinterval","threshold", "quadrant"), ...){
  
  gate_type <- match.arg(gate_type)
  
  # Two fluorescent channels must be supplied
  if(missing(channels) | !length(channels) %in% c(1,2)){
    stop("Two fluorescent channels must be specified to draw polygon gate.")
  }

  # Truncate flowFrame if min and max arguments are supplied
  if(!(is.null(min) && is.null(max))){
    fr <- openCyto::.truncate_flowframe(fr, channels = channels, min = min, max = max)
  }

  # Determine vertices of polygon using DrawGate
  pts <- DrawGate(fr, channels, gate_type = gate_type)

  
  if(gate_type == "polygon"){
 
  # Construct polygonGate
    
  gate <- polygonGate(.gate = pts)
  
  }else if(gate_type == "rectangle"){
  
  # Construct rectangle gate  

  gate <- rectangleGate(.gate = pts)
  
  }else if(gate_type == "interval"){
  
  # Construct interval gate

  gate <- rectangleGate(.gate = pts)
  
  }else if(gate_type == "multiinterval"){
    
  # Construct multiple interval gates
    
  channel <- channels
    
  gate <- lapply(pts, function(pts){
    pts <- data.frame(x = pts[,1])
    colnames(pts) <- channel
    rectangleGate(.gate = pts)
  })
  gate <- filters(gate)
  
  }else if(gate_type == "threshold"){
  
  # Construct recatngle gate
    
  gate <- rectangleGate(.gate = pts)
  
  }else if(gate_type == "quadrant"){
  
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
    
    gate <- filters(list(q1,q2,q3,q4))
  }
    
  return(gate)
}
