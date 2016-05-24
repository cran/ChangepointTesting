#******************************************************************************#
# Benjamini Hochberg method                                                    #
#******************************************************************************#
#                                                                              #
#  Inputs :                                                                    #
#                                                                              #
#   pvalues : vector of p-values                                               # 
#                                                                              #
#   alpha   : single significant level at which the FDR is to be 'controlled'  #
#                                                                              #
#  Outputs :                                                                   #
#                                                                              #
#  A list is returned                                                          #
#                                                                              #
#    ind : a vector of indicator functions. (1) rejected (0) not-rejected      #
#                                                                              #
#    numAlt : number of rejected hypotheses                                    #
#                                                                              #
#    propAlt : proportion of hypotheses rejected.                              #
#                                                                              #
#******************************************************************************#
BenjaminiHochberg <- function(pvalues, alpha) {

  adjusted.p_value <- p.adjust(pvalues, method="BH")

  indicator.bh <- adjusted.p_value <= alpha

  numAlt <- sum(indicator.bh)
  propAlt <- numAlt/length(pvalues)

  return( list("ind" = as.numeric(indicator.bh), 
               "numAlt" = numAlt, 
               "propAlt" = propAlt) )

}
