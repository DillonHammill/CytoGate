drawGate <- function(x) UseMethod("drawGate")

drawGate.default <- function(x, ...){
  
  warning(paste("drawGate does not know how to handle object of class ", 
                class(x), 
                "and can only be used on classes flowFrame, flowSet, and GatingSet"))
  
}

drawGate.flowFrame <- function(x, channels, alias, subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  fr <- x
  
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
    
    if(dev.cur() %in% c(1,2)){ # Graphics device not open
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

  gates <- filters(gates)
  return(gates)
  
}

drawGate.flowSet <- function(x, channels, alias, subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  fr <- as(x, "flowFrame")
  
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
    
    if(dev.cur() %in% c(1,2)){ # Graphics device not open
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
  
  gates <- filters(gates)
  return(gates)
  
}

drawGate.GatingSet <- function(x, channels, parent = "root", alias, file = "GatingTemplate.csv",subSample = NULL, gate_type = NULL, axis = "x", adjust = 1.5, plot = TRUE, labs = TRUE, ...){
  
  gs <- x
  fs <- getData(x, parent)
  fr <- as(fs, "flowFrame")
  
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
    
    if(dev.cur() %in% c(1,2)){ # Graphics device not open
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
  
  gates <- filters(gates)
  
  # Prepare arguments for GatingTemplate Entry
  if(length(channels) == 2){
    channels <- paste(channels, collapse = ",")
  }
  
  if(length(alias) > 1){
    alias <- paste(alias, collapse = ",")
  }
  
  if(length(alias) > 1){
    pop <- "*"
  }else{
    pop = "+"
  }
  
  
  # Check if GatingTemplate exists - if not create one
  if ("template" %in% ls(envir = .GlobalEnv)) {

    template <- rbind(get("template", envir = .GlobalEnv), add_pop(
      gs, alias = alias, parent = parent, pop = pop, dims = channels, gating_method = "manualGate",
      gating_args = list(gate = gates)
    ))
    
  } else {
    
    template <- add_pop(
      gs, alias = alias, parent = parent, pop = pop, dims = channels, gating_method = "manualGate",
      gating_args = list(gate = gates)
    )
    
  }
  
  template <<- template
  write.csv(template, file)
  
  return(gs)
  
}
