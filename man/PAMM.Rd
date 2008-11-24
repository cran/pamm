\name{PAMM}
\alias{PAMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation function to assess power of mixed models }
\description{ Given a specific varaince-covariance structure for random
  effect, the function simulate different group size and assess p-values and power of random
  intercept and random slope }
\usage{
PAMM(numsim, group, repl, randompart, fixed = c(0, 1, 0), intercept = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numsim}{ number of simulation for each step }
  \item{group}{ number of group (individuals). Could be specified as a vector }
  \item{repl}{ number of replicates (observations) per group . Could be specified as a vector }
  \item{randompart}{ vector of lenght 4 or 5, with 1) variance component
    of intercept (VI); 2) variance component of slope (VS); 3) residual
    variance (VR); 4) relation between random intercept and random
    slope; 5) "cor" or "cov" determine if the relation 4) between I ans S is a correlation or a covariance (set to correlation by default) }
  \item{fixed}{ vector with mean, variance and estimate of fixed effect to simulate. c(0,1,0) by default }
  \item{intercept}{a numeric value giving the expected intercept value. Default:0 }
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
\seealso{\code{\link{EAMM}}, \code{\link{SSF}}, \code{\link{plot.PAMM}} }
\examples{
   \dontrun{
  ours=PAMM(numsim=10,group=c(seq(10,50,10),100),repl=c(2,4,6),
           randompart=c(0.4,0.1,0.5,0.1),fixed=c(0,1,0.7))       
  plot(ours,"both")
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc }% __ONLY ONE__ keyword per line
