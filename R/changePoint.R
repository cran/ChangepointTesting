#**********************************************************************#
# changePoint : public function implementing change point method       #
#**********************************************************************#
#                                                                      #
# Input                                                                #
#                                                                      #
# pvalues  : an object of class numeric.                               #
#            a vector of pvalues                                       #
#                                                                      #
# alpha    : an object of class numeric                                #
#            level of significance for determining the critical value  #
#                                                                      #
# km       : an object of class numeric                                #
#            size of the window defining the neighborhood in left and  #
#            right distances                                           #
#                                                                      #
# lm       : an object of class numeric                                #
#            size of window defining the neighborhood in variance est  #
#                                                                      #
# compare  : an object of class character                              #
#            one of {Both, FDRL, BH} indicating what methods are to be #
#            compared to the changePoint method.                       #
#                                                                      #
# fdrlWindow : an object of class numeric                              #
#            size of the window defining the neighborhood              #
#                                                                      #
# fdrlNStep  : an object of class numeric                              #
#            the number of threshold values to consider                #
#                                                                      #
# fdrlLambda : an object of class numeric                              #
#            tuning constant                                           #
#                                                                      #
# Output                                                               #
#                                                                      #
#  An object of class changePoint is returned                          #
#                                                                      #
#**********************************************************************#
changePoint <- function(pvalues, 
                        alpha, 
                        km, 
                        lm,
                        compare="BOTH",
                        fdrlWindow = 3,
                        fdrlNStep = 300,
                        fdrlLambda = 0.1){

  compare <- toupper(compare)

  if("BOTH" %in% compare) compare <- c("BH","FDRL")

  if(!all(compare %in% c("BH","FDRL","NONE"))) {
    stop("method specified for comparison is not available")
  }

  #------------------------------------------------------------------#
  # Determine the number of p-values passed in.                      #
  #------------------------------------------------------------------#
  m <- length(pvalues)

  #------------------------------------------------------------------#
  # Calculate the variance estimator                                 #
  #------------------------------------------------------------------#
  mean_pvalue <- mean(pvalues)
  pvalueCumSum <- cumsum(pvalues)

  #------------------------------------------------------------------#
  # j = 1 term                                                       #
  #------------------------------------------------------------------#
  sigmasq <- (pvalueCumSum[lm]/lm - mean_pvalue)^2

  #------------------------------------------------------------------#
  # j = 2, ..., m-lm+1 terms                                         #
  #------------------------------------------------------------------#
  temp <- (pvalueCumSum[{lm+1L}:m] - pvalueCumSum[1L:{m-lm}])/lm - 
          mean_pvalue
  sigmasq <- sigmasq + sum(temp*temp)
  sigmasq <- sigmasq * lm / {m - lm + 1.0}

  cat("km=",km,"lm=",lm,"sigmasq=",sigmasq,"\n")

  #------------------------------------------------------------------#
  # Simulation assisted approach to determine the critical value     #
  #------------------------------------------------------------------#
  iter <- 10000L

  tempFunc <- function(x) {

    cumP <- cumsum(runif(n = m, min = 0, max = 1))
    maxL <- max( abs( c( (0.5 - cumP[km]/km), 
                         (0.5 - (cumP[(km+1L):m] - 
                                 cumP[1L:(m-km)])/km ) ) ))

    return( maxL )
  }

  dist <- sapply(1L:iter, tempFunc)

  gamma <- quantile(dist,(1.0-alpha))

  critical <- gamma*sqrt(12.0*sigmasq)
  names(critical) <- NULL

  rm(dist)

  cat("gamma=",gamma,"critical=",critical,"\n")

  #------------------------------------------------------------------#
  # Calculate left and right differences                             #
  #------------------------------------------------------------------#
  lbnds <- 1L:{m-1L} - km
  lbnds[lbnds < 1L] <- 1L
  rbnds <- 1L:{m-1L} + km
  rbnds[rbnds > m] <- m

  Li <- c(0.0, 
          abs((pvalueCumSum[1L:{m-1L}] - pvalueCumSum[lbnds])/km - 0.5))
  Ri <- c(abs(pvalueCumSum[km]/km - 0.5), 
          abs((pvalueCumSum[rbnds] - pvalueCumSum[1L:{m-1L}])/km - 0.5))

  Qi <- ({Li > critical} + {Ri > critical})*2.0

  #------------------------------------------------------------------#
  # Within each window, determine dominant critical test value       #
  #------------------------------------------------------------------#
  most <- numeric(m)
  halfWin <- floor(km/2)

  for( i in 1L:m ) {

    bnds <- max(1L, i-halfWin) : min(i+halfWin, m)

    tests <- c(sum(Qi[bnds] < 1L),
               sum(Qi[bnds] < 3L),
               sum(Qi[bnds] < 5L))

    most[i] <- which.max( c(tests[1L], 
                            tests[2L]-tests[1L], 
                            tests[3L]-tests[2L]) ) - 1L

  }

  #------------------------------------------------------------------#
  # For transition of 1 -> 2, determine change point                 #
  #------------------------------------------------------------------#
  count12 <- which( (most[1:{m-1L}] > 0.5 & most[1:{m-1L}] < 1.5) & 
                    (most[2:m] > 1.5) )

  cp <- NULL
  for( i in count12 ) {

    k <- i
    while( most[k] == 1L ) {
      k <- k - 1L
      if( k < 1L ) break
    }

    k <- k + 1L

    if( (i - k) <= km/2) next

    cut <- which.max(Li[k:i]*(Li[k:i] > critical) + 
                     Ri[k:i]*(Ri[k:i] > critical))

    #--------------------------------------------------------------#
    # If the first element is the max, shift cutoff to avoid       #
    # moving out of the window.                                    #
    #--------------------------------------------------------------#
    if( cut == 1L ) cut <- 2L

    cp <- c(cp, k + cut - 2L)

  }

  #------------------------------------------------------------------#
  # For transition of 2 -> 1, determine change point                 #
  #------------------------------------------------------------------#
  count21 <- which( (most[1:{m-1L}] > 1.5) & 
                    ((most[2:m] > 0.5) & (most[2:m] < 1.5)) )
  for( i in count21 ) {

    j <- i + 1L

    k <- j
    while( most[k] == 1L ) {
      k <- k + 1L
      if( k > m ) break
    }

    k <- k - 1L
    
    if( (k - j) <= km/2) next

    cut <- which.max(Li[j:k]*(Li[j:k] > critical) + 
                     Ri[j:k]*(Ri[j:k] > critical))

    #--------------------------------------------------------------#
    # If the first element is the max, shift cutoff to avoid       #
    # moving out of the window.                                    #
    #--------------------------------------------------------------#
    if( cut == 1L ) cut <- 2L

    #--------------------------------------------------------------#
    # if this change point was identified for the 0 -> 1 transition#
    # indicates that there are no 2s. Remove change point.         #
    #--------------------------------------------------------------#
    dup <- cp == (j + cut - 2L)

    if( any(dup) ) {
      cp <- cp[!dup]
      next
    }

    cp <- c(cp, j + cut - 2L)

  }


  Qi[] <- 0L

  if( length(cp) > 1.0 ) {
    #--------------------------------------------------------------#
    # Sort change points                                           #
    #--------------------------------------------------------------#
    cp <- sort(cp)

    #--------------------------------------------------------------#
    # Generate indicator function                                  #
    #--------------------------------------------------------------#
    for( i in seq(1L, length(cp),2) ) {
      Qi[ {cp[i]+1L}:{cp[i+1L]-1L} ] <- 1L
    } 
  }

  numAlt <- sum(Qi)
  propAlt <- numAlt/m

  if("BH" %in% compare){
    indicator.bh <- BenjaminiHochberg(pvalues = pvalues, alpha = alpha)
  } else {
    indicator.bh <- NULL
  }

  if("FDRL" %in% compare){
    indicator.fdrl <- FDRLMethod(pvalues = pvalues, 
                                 window = fdrlWindow, 
                                 alpha = alpha, 
                                 nstep = fdrlNStep, 
                                 lambda = fdrlLambda)
  } else {
    indicator.fdrl <- NULL
  }

  obj <- new("changePoint",
             CW = Qi,
             chgPts = cp,
             pi_alt = propAlt,
             num_alt = numAlt,
             FDRL = indicator.fdrl$ind, 
             BH = indicator.bh$ind,
             gammaStar = critical,
             sigmaSq = sigmasq,
             pVals = pvalues)

  show(obj)

  return( obj )

}
