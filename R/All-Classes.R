setClassUnion("numeric or NULL", members = c("numeric","NULL"))

setClass("changePoint", 
  representation(CW = "numeric or NULL", 
                 chgPts = "numeric or NULL",
                 pi_alt = "numeric",
                 num_alt = "numeric",
                 FDRL = "numeric or NULL", 
                 BH = "numeric or NULL",
                 gammaStar = "numeric",
                 sigmaSq = "numeric",
                 pVals = "numeric"),
  prototype(CW = 0, 
            chgPts = 0,
            pi_alt = -1,
            num_alt = -1,
            FDRL = 0,
            BH = 0,
            gammaStar = -1,
            sigmaSq = -1,
            pVals = 0)
)


setMethod(f="show",
          signature=c(object="changePoint"),
          definition = function(object){
                         result <- matrix(data = 0.0,
                                          nrow = 1L,
                                          ncol = 2L,
                                          dimnames = list("changePoint",
                                                          c("# Rejected","% Rejected")))

                         result[1L,] <- c(object@num_alt, object@pi_alt)

                         if( !is(object@FDRL,"NULL") ) {
                           FDR_L <- c(sum(object@FDRL), 
                                      sum(object@FDRL)/length(object@FDRL))
                           result <- rbind(result, FDR_L)
                         }
                         if( !is(object@BH,"NULL") ) {
                           BH <- c(sum(object@BH), 
                                      sum(object@BH)/length(object@BH))
                           result <- rbind(result, BH)
                         }
                         print(result)
                         cat("gamma*:  ", object@gammaStar, ".\n", sep = "")
                         cat("sigma^2: ", object@sigmaSq, ".\n", sep = "")
                       } )

if (!isGeneric("plot"))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod(f="plot",
          signature=c(x="changePoint", y = "missing"),
          definition = function(x, y, logp = FALSE, ...){

                         blk <- blocks(x)
                         if( length(blk) < 0.5 ) {
                           cat("None detected.\n")
                           return
                         }

                         np <- length(blk)
                         nr <- ceiling(np / 4)

                         hold <- par(mfrow=c(nr, min(np,4)))

                         if( !logp ) {
                           y <- x@pVals
                           ny <- "p-value"
                         } else {
                           y <- -log(x@pVals)
                           ny <- "-log(p)"
                         }

                         for( i in 1L:np ) {
                           rng <- blk[[i]]
                           plot(x = rng, 
                                y = y[rng],
                                xlab = "position",
                                ylab = ny, ...)
                         }
                         
                         par(hold)
                       } )


setMethod(f="print",
          signature=c(x="changePoint"),
          definition = function(x){

                         result <- matrix(data = 0.0,
                                          nrow = 1L,
                                          ncol = 2L,
                                          dimnames = list("changePoint",
                                                          c("# Rejected","% Rejected")))

                         result[1L,] <- c(x@num_alt, x@pi_alt)

                         if( !is(x@FDRL,"NULL") ) {
                           FDR_L <- c(sum(x@FDRL), 
                                      sum(x@FDRL)/length(x@FDRL))
                           result <- rbind(result, FDR_L)
                         }
                         if( !is(x@BH,"NULL") ) {
                           BH <- c(sum(x@BH), 
                                   sum(x@BH)/length(x@BH))
                           result <- rbind(result, BH)
                         }
                         print(result)
                         cat("gamma*:  ", x@gammaStar, ".\n", sep = "")
                         cat("sigma^2: ", x@sigmaSq, ".\n", sep = "")
                       } )

