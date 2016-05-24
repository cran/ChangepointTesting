#******************************************************************************#
# FDRL : public function implementing FDR_L method of Zhang et al.             #
#******************************************************************************#
#                                                                              #
# Input                                                                        #
#                                                                              #
# pvalues  : an object of class numeric.                                       #
#            a vector of pvalues                                               #
#                                                                              #
# window   : an object of class numeric                                        #
#            size of the window defining the neighborhood                      #
#                                                                              #
# alpha    : an object of class numeric                                        #
#            the level of significance for determining the critical value      #
#                                                                              #
# nstep    : an object of class numeric                                        #
#            the number of threshold values to consider                        #
#                                                                              #
# lambda   : an object of class numeric                                        #
#            tuning constant                                                   #
#                                                                              #
# Output                                                                       #
#                                                                              #
#  A list is returned                                                          #
#                                                                              #
#    ind : a vector of indicator functions. (1) rejected (0) not-rejected      #
#                                                                              #
#    threshold : critical value                                                #
#                                                                              #
#    numAlt : number of rejected hypotheses                                    #
#                                                                              #
#    propAlt : proportion of hypotheses rejected.                              #
#                                                                              #
#******************************************************************************#
FDRLMethod <- function(pvalues, window, alpha, nstep = 300, lambda = 0.1) {

  if( !is(window,"integer") ) window <- as.integer(round(window,0))

  if( lambda < 0.0 ) stop("lambda must be [0,1]")

  tol <- 1.5e-8

  m <- length(pvalues)
  pstar <- numeric(length = m)

  for( i in 1L:m ) {
    inx <- max(1L,(i-window)):min((i+window),m)
    pstar[i] <- median(pvalues[inx])
  }

  #--------------------------------------------------------------------------#
  # Number of non-rejections with a threshold of lambda                      #
  # p > lambda                                                               #
  #--------------------------------------------------------------------------#
  Wlambda <- sum( (pstar - lambda) > tol )

  #--------------------------------------------------------------------------#
  # Denominator of Ghat Eq 3.3                                               #
  # tst1 = p >= 0.5                                                          #
  # tst2 = p >  0.5                                                          #
  # tst1 - tst2 = p == 0.5                                                   #
  # 2(p>0.5) + (p==0.5) = 2*tst2 + (tst1 - tst2) = tst2 + tst1               #
  #--------------------------------------------------------------------------#
  tst1 <- sum((pstar - 0.5) > -tol)
  tst2 <- sum((pstar - 0.5) >  tol)
  denom <- 1.0/(tst2 + tst1)


  #--------------------------------------------------------------------------#
  # Ghat for threshold lambda                                                #
  #--------------------------------------------------------------------------#
  if( (lambda - 0.5) < tol ) {
    gHatLambda <- sum((pstar - (1.0 - lambda)) > -tol) * denom
  } else {
    gHatLambda <- 1.0 - sum((pstar - lambda) > tol) * denom
  }

  tVec <- (1L:nstep)/nstep

  thresh <- 0.0

  lim <- floor(nstep/2)

  for(i in 1L:lim) {

    #----------------------------------------------------------------------#
    # Number of rejections with a threshold of t                           #
    # p <= t                                                               #
    #----------------------------------------------------------------------#
    Rt <- max(1.0, sum(pstar <= (tVec[i]+tol)))

    #----------------------------------------------------------------------#
    # Ghat for threshold t <= 0.5                                          #
    #----------------------------------------------------------------------#
    gHatt <- sum( (pstar - (1.0-tVec[i])) > - tol) * denom

    if( ( Wlambda*gHatt / (Rt*(1.0-gHatLambda)) ) < alpha ) thresh <- i/nstep

  }

  for(i in (lim+1L):nstep) {

    #----------------------------------------------------------------------#
    # Number of rejections with a threshold of t                           #
    # p <= t                                                               #
    #----------------------------------------------------------------------#
    Rt <- max(1.0, sum(pstar <= (tVec[i]+tol)))

    #----------------------------------------------------------------------#
    # Ghat for threshold t > 0.5                                           #
    #----------------------------------------------------------------------#
    gHatt <- 1.0 - sum( (pstar - tVec[i]) > tol) * denom

    if( ( Wlambda*gHatt / (Rt*(1.0-gHatLambda)) ) < alpha ) thresh <- i/nstep

  }

  testresult <- as.numeric(pstar <= thresh)

  numAlt <- sum(testresult)
  propAlt <- numAlt/m

  return( list("ind" = testresult, 
               "threshold" = thresh, 
               "numAlt" = numAlt, 
               "propAlt" = propAlt) )

}

