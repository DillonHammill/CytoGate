#' Draw 2D Gates on Plots of Flow Cytometry Data
#'
#' \code{DrawGate} facilitates manual drawing of polygon gates onto flow cytometry plots.
#'
#' @param fr a \code{flowFrame} object containing the flow cytometry data for plotting and gating.
#' @param channels a vector of length 2 indicating the fluorescent channels to be used to construct 2D plot
#' and \code{polygonGate}.
#'
#' @return a \code{dataframe} object containing the vertices of the polygon gate.
#'
#' @keywords manual, gating, draw, FlowJo, polygonGate, openCyto
#' @import flowDensity
#' @export
#'
#' @author Dillon Hammill (Dillon.Hammill@anu.edu.au)
DrawGate <- function(fr, channels){

  # Check that length(channels) %in% c(1,2)
  if(!length(channels) == 2){
    stop("Please supply either 2 fluorescent channels for gating")
  }

  # Extract data for plotting and gating
  x <- exprs(fr[,channels])

  # Determine whether R is being run in RStudio
  if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    # if TRUE we need to open X11() interactive graphics device
    X11()
  }

  # Plot the data for gating use flowDensity::plotDens - locator() only works for base graphics
  cat("Draw 2D polygon gate around population. \n")

  flowDensity::plotDens(fr, channels = channels, cex = 3)

  # Extract points of drawn gate
  pts <- locator(type = "l", lwd = 2)

  if (length(pts$x) < 3) stop("A minimum of 3 points is required to construct a polygon gate.")
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
gate_draw <- function(fr, pp_res, channels, filterId = "", gate_range = NULL, min = NULL, max = NULL, ...){

  # Two fluorescent channels must be supplied
  if(missing(channels) | length(channels) != 2){
    stop("Two fluorescent channels must be specified to draw polygon gate.")
  }

  # Truncate flowFrame if min and max arguments are supplied
  if(!(is.null(min) && is.null(max))){
    fr <- openCyto::.truncate_flowframe(fr, channels = channels, min = min, max = max)
  }

  # Determine vertices of polygon using DrawGate
  pts <- DrawGate(fr, channels)

  # Construct polygonGate
  polygonGate(.gate = pts)

}
