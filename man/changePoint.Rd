\name{changePoint}
\alias{changePoint}
\alias{plot,changePoint}
\alias{BH,changePoint}
\alias{CW,changePoint}
\alias{FDRL,changePoint}
\alias{critical,changePoint}
\alias{numAlt,changePoint}
\alias{pAlt,changePoint}
\alias{sigmaSq,changePoint}
\alias{BH}
\alias{blocks}
\alias{CW}
\alias{FDRL}
\alias{critical}
\alias{numAlt}
\alias{piAlt}
\alias{sigmaSq}
\alias{changePts}
\title{ Change Point Estimation for Clustered Signals }
\description{ A multiple testing procedure for clustered alternative 
              hypotheses. It is assumed that the p-values under the 
              null hypotheses follow U(0,1) and that the distributions
              of p-values from the alternative hypotheses are
              stochastically smaller than U(0,1).
              By aggregating information, this method is more
              sensitive to detecting signals of low magnitude than
              standard methods. Additionally,
              sporadic small p-values appearing within a
              null hypotheses sequence are avoided 
              by averaging on the neighboring p-values. }
\usage{
changePoint(pvalues, alpha, km, lm, compare = "BOTH", fdrlWindow = 3, 
fdrlNStep = 300, fdrlLambda = 0.1)
}
\arguments{
  \item{pvalues   }{an object of class numeric.
                    A vector of p-values.}
  \item{alpha     }{an object of class numeric. The significant level for
                    the estimation of the critical value, gamma*.}
  \item{km        }{an object of class numeric. The size of the window 
                    defining the neighborhood in left and right distances.}
  \item{lm        }{an object of class numeric. The size of the 
                    window defining the 
                    neighborhood in the long-run variance estimation.}
  \item{compare   }{one of ("FDRL", "BH", "Both", "None"). 
                    In addition to the Cao-Wu method, 
                    obtain significance indicators using
                    the FDR_L method (FDRL) 
                    (Zhang et al., 2011),
                    the Benjamini-Hochberg method (BH), 
                    (Benjamini andHochberg, 1995),
                    "both" the FDRL and the BH methods, or 
                    do not consider alternative methods (none).}
  \item{fdrlWindow}{an object of class numeric.
                    If FDR_L method requested, the size of the window 
                    defining the neighborhood.}
  \item{fdrlNStep }{an object of class numeric.
                    If FDR_L method requested, the number of threshold
                    values to consider.}
  \item{fdrlLambda}{and object of class numeric.
                    If FDR_L method requested, the tuning constant.}
             
}
\details{
The comparison capability is included only for convenience and 
reproducibility 
of the original manuscript. The Benjamini-Hochberg and FDR_L
methods cannot be accessed outside of the changePoint function.

The following methods retrieve individual results from a changePoint
object, x: 

    \code{BH(x)}: 
              Retrieves a vector of integer values. 
              An element is 1 if the null hypothesis is rejected
              by the Benjamini-Hochberg (1995) method.

    \code{blocks(x)}: 
              Retrieves a list, each element of which is 
              a vector of integer values. 
              Each vector contains the indices of 
              an alternative hypothesis block.

    \code{CW(x)}: 
              Retrieves a vector of integer values. 
              An element is 1 if the null hypothesis is rejected 
              by the Cao-Wu change point (2015) method.

    \code{changePts(x)}:
              Retrieves a vector of integer values.
              The vector of change points identified by the
              Cao-Wu (2015) method. If no change points are
              identified, NULL is returned.

    \code{FDRL(x)}:  
              Retrieves a vector of integer values. 
              Elements are 1 if the null hypothesis is rejected 
              by the FDR_L (Zhang et al. 2011) method. 

    \code{critical(x)}: Retrieves the
        estimated critical value for testing used by
        the Cao-Wu (2015) method.

    \code{numAlt(x)}: Retrieves the
        estimated number of alternative hypotheses
        obtained by the Cao-Wu (2015) method. 

    \code{piAlt(x)}: Retrieves the
        estimated proportion of alternative hypotheses
        obtained by the Cao-Wu (2015) method. 

    \code{plot(x, y, logp, ...)}: Generates plots
        of -log(p) vs position or p-value vs position for
        each alternative hypothesis block obtained
        by the Cao-Wu (2015) method. logp is TRUE/FALSE
        indicating if -log(p)/p-values are plotted on the y-axis.

    \code{sigmaSq(x)}: Retrieves the
        estimated variance used to determine the critical value of
        the Cao-Wu (2015) method.

}
\value{
Returns an object of class \code{changePoint}.
}
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
Maintainer: Shannon T. Holloway <sthollow@ncsu.edu>
}
\examples{

m <- 5000

T <- c(rep(0.1,  75), rep( 0.8, 75), rep(1.8, 100), 
       rep(0.0,2250), rep(-1.5,250), rep(0.0,2250)) + 
     rnorm(m, mean=0.0, sd = 1.0)

pv <- 2.0*(1.0-pnorm(abs(T)))

res <- changePoint(pvalues = pv, 
                   alpha = 0.05, 
                   km = {log(m)}^2, 
                   lm = m^{1/4}, 
                   compare = "Both")

print(changePts(res))

print(head(cbind(BH(res),FDRL(res),CW(res))))

}

