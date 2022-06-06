\name{changePoint-class}
\Rdversion{1.0}
\docType{class}
\alias{changePoint-class}
\alias{BH,changePoint-method}
\alias{blocks,changePoint-method}
\alias{CW,changePoint-method}
\alias{changePts,changePoint-method}
\alias{critical,changePoint-method}
\alias{numAlt,changePoint-method}
\alias{piAlt,changePoint-method}
\alias{FDRL,changePoint-method}
\alias{sigmaSq,changePoint-method}
\alias{print,changePoint-method}
\alias{plot,changePoint,missing-method}
\alias{show,changePoint-method}

\title{Class \code{"changePoint"}}

\description{Value object returned by call to \code{changePoint()}. 
}
\section{Objects from the Class}{This object should not be created by users.}

\section{Slots}{
  \describe{
    \item{\code{CW}:        }{Object of class \code{numeric} or NULL. 
                              A vector of 1/0 values; 1 indicates  
                              that hypothesis was rejected by the  
                              Cao-Wu method.}
    \item{\code{chgPts}:    }{Object of class \code{numeric} or NULL.  
                              The vector of change points identified  
                              by the Cao-Wu method. If no change  
                              points are identified, NULL.}
    \item{\code{pi_alt}:    }{Object of class \code{numeric}. 
                              The estimated proportion of alternative  
                              hypotheses calculated using the Cao-Wu  
                              method.}
    \item{\code{num_alt}:   }{Object of class \code{numeric}. 
                              The estimated number of alternative hypotheses
                              calculated using the Cao-Wu method.}
    \item{\code{FDRL}:      }{Object of class \code{numeric} or NULL.  
                              A vector of 1/0 values; 1 indicates that  
                              hypothesis was rejected by the FDR_L  
                              method.}
    \item{\code{BH}:        }{Object of class \code{numeric} or NULL.  
                              A vector of 1/0 values; 1 indicates that  
                              hypothesis was rejected by the FDR_L  
                              method.}
    \item{\code{gammaStar}: }{Object of class \code{numeric}.  The  
                              estimated critical value for testing  
                              used by the Cao-Wu method.}
    \item{\code{sigmaSq}:   }{Object of class \code{numeric}.  The  
                              estimated variance used to determine  
                              the critical value of the Cao-Wu method.}
    \item{\code{pVals}:     }{Object of class \code{numeric}.  The  
                              original p-values provided as input.}
  }
}
\section{Methods}{
  \describe{
    \item{BH}{\code{signature(x = "changePoint")}: 
              Retrieves a vector of integer values. 
              An elements is 1 if the null hypothesis is rejected
              by the Benjamini-Hochberg (1995) method.}

    \item{blocks}{\code{signature(x = "changePoint")}: 
              Retrieves a list, each element of which is 
              a vector of integer values. 
              Each vector contains the indices of 
              an alternative hypothesis block.}

    \item{CW}{\code{signature(x = "changePoint")}: 
              Retrieves a vector of integer values. 
              An element is 1 if the null hypothesis is rejected 
              by the Cao-Wu change point (2015) method.}

    \item{changePts}{\code{signature(x = "changePoint")}:
              Retrieves a vector of integer values.
              The vector of change points identified by the
              Cao-Wu (2015) method. If no change points are
              identified, NULL is returned.}

    \item{FDRL}{\code{signature(x = "changePoint")}:  
              Retrieves a vector of integer values. 
              Elements are 1 if the null hypothesis is rejected 
              by the FDR_L (Zhang et al. 2011) method. }

    \item{critical}{\code{signature(x = "changePoint")}: Retrieves the
        estimated critical value for testing used by
        the Cao-Wu (2015) method.}

    \item{numAlt}{\code{signature(x = "changePoint")}: Retrieves the
        estimated number of alternative hypotheses
        obtained by the Cao-Wu (2015) method. }

    \item{piAlt}{\code{signature(x = "changePoint")}: Retrieves the
        estimated proportion of alternative hypotheses
        obtained by the Cao-Wu (2015) method. }

    \item{plot}{\code{signature(x = "changePoint", y = "missing", 
                                logp = FALSE, ...)}: 
        Generates x-y plots
        of -log(p) vs position or p-value vs position for
        each alternative hypothesis block obtained
        by the Cao-Wu (2015) method. logp is TRUE/FALSE
        indicating if -log(p)/p-values are plotted on the y-axis.}

    \item{sigmaSq}{\code{signature(x = "changePoint")}: Retrieves the
        estimated variance used to determine the critical value of
        the Cao-Wu (2015) method.}

}}
\references{
Benjamini, Y. and Hochberg, Y. (1995).
Controlling the false discovery rate: A practical and powerful approach to 
multiple testing.
{Journal of the Royal Statistical Society: Series B}, 57, 289--300.

Cao, H. and Wu, W. B. (2015)
Changepoint estimation: Another look at multiple testing problems.
Biometrika, 102, 974--980.

Zhang, C., Fan, J., and Yu, T. (2011).
Multiple testing via FDRL for large-scale imaging data.
Anals of Statistics, 39, 613--642.
}
\author{
Hongyuan Cao, Wei Biao Wu, and Shannon T.  Holloway
Maintainer: Shannon T. Holloway <shannon.t.holloway@gmail.com>
}
\examples{
showClass("changePoint")
}
\keyword{classes}
