\name{SSF}
\alias{SSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation function to assess power of mixed models }
\description{ Given a specific total number of observations and variance-covariance structure for random
  effect, the function simulates different association of number of
  group and replicates (giving the specified sample size) and assess p-values and power of random
  intercept and random slope }
\usage{
SSF(numsim, tss, nbstep = 10, randompart, fixed = c(0, 1, 0), intercept = 0, exgr = NA, exrepl = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numsim}{ number of simulation for each step }
  \item{tss}{ total sample size (nb group * nb replicates) }
  \item{nbstep}{ number of group*replicates associations to simulate }
  \item{randompart}{ vector of lenght 4 or 5 with 1) variance component
    of intercept (VI); 2) variance component of slope (VS); 3) residual
    variance (VR); 4) relation between random intercept and random
    slope; 5) "cor" or "cov" determine id the relation between I ans S is correlation or covariance (set to correlation by default) }
  \item{fixed}{ vector of lenght 3 with mean, variance and estimate of fixed effect to simulate }
  \item{intercept}{a numeric value giving the expected intercept value. Default:0 }
  \item{exgr}{ a vector specifying minimum and maximum value for number of group. Default:c(2,tss/2) }
  \item{exrepl}{ a vector specifying minimum and maximum value for number of replicates. Default:c(2,tss/2) }
  }
\details{
 P-values for random effects are estimated using a log-likelihood ratio
 test between two models with and without the effect. Power represent
 the percentage of simulations providing a significant p-value for a
 given random structure
}
\value{
  data frame reporting estimated P-values and power with CI for random
  intercept and random slope  
 }
\references{ ... }
\author{ Julien Martin }
\note{ }
\section{Warning }{the simulation is based on a balanced data set with
  unrelated group}
\seealso{\code{\link{PAMM}}, \code{\link{EAMM}}, \code{\link{plot.SSF}} }
\examples{
   \dontrun{
  ours<- SSF(100,200,10,c(0.4,0.1,0.6,0))
  plot(ours)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc }% __ONLY ONE__ keyword per line
