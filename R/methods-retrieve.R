if(!isGeneric("BH")){
  setGeneric(name="BH", 
             def=function(x,...){standardGeneric("BH")})
}

setMethod("BH",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    if(is.null(x@BH)){
      cat("BH method not calculated for provided object.\n")
      return(NULL)
    } else {
      return(x@BH)
    }
  }
)

if(!isGeneric("blocks")){
  setGeneric(name="blocks", 
             def=function(x,...){standardGeneric("blocks")})
}

setMethod("blocks",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    if( length(x@chgPts) > 1.5 ) {
      blk <- list()
      j <- 1L
      for( i in seq(from=1L, to=length(x@chgPts), by=2L) ) {
        blk[[j]] <- x@chgPts[i]:{x@chgPts[i+1L]-1L}
        j <- j + 1L
      }
    } else {
      blk <- NULL
    }

    return(blk)
  }
)

if(!isGeneric("CW")){
  setGeneric(name="CW", 
             def=function(x,...){standardGeneric("CW")})
}

setMethod("CW",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    return(x@CW)
  }
)

if(!isGeneric("changePts")){
  setGeneric(name="changePts", 
             def=function(x,...){standardGeneric("changePts")})
}

setMethod("changePts",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    return(x@chgPts)
  }
)

if(!isGeneric("FDRL")){
  setGeneric(name="FDRL", 
             def=function(x,...){standardGeneric("FDRL")})
}

setMethod("FDRL",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    if(is.null(x@FDRL)){
      cat("FDR_L method not calculated for provided object.\n")
      return(NULL)
    } else {
      return(x@FDRL)
    }
  }
)

if(!isGeneric("critical")){
  setGeneric(name="critical", 
             def=function(x,...){standardGeneric("critical")})
}

setMethod("critical",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    return(x@gammaStar)
  }
)



if(!isGeneric("numAlt")){
  setGeneric(name="numAlt", 
             def=function(x,...){standardGeneric("numAlt")})
}

setMethod("numAlt",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    return(x@num_alt)
  }
)

if(!isGeneric("piAlt")){
  setGeneric(name="piAlt", 
             def=function(x,...){standardGeneric("piAlt")})
}

setMethod("piAlt",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    return(x@pi_alt)
  }
)

if(!isGeneric("sigmaSq")){
  setGeneric(name="sigmaSq", 
             def=function(x,...){standardGeneric("sigmaSq")})
}

setMethod("sigmaSq",
  signature = c(x="changePoint"),
  definition = function(x, ...){
    return(x@sigmaSq)
  }
)


