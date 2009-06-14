\name{EAMM}
\alias{EAMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SIMULATION FUNCTION FOR EXPLORATORY POWER ANALYSIS FOR RANDOM EFFECTS}
\description{
 Given a specific sample size (fixed number of group and replicates per group), the function simulate different variance-covariance structure and assess p-values and power of random intercept and random slope using lmer.  
}
\usage{
EAMM(numsim, group, repl, fixed = c(0, 1, 0), VI = seq(0.05, 0.95,0.05),
       VS = seq(0.05, 0.5, 0.05), CoIS = 0, relIS = "cor")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numsim}{number of simulation for each step}
  \item{group}{number of group (individuals)}
  \item{repl}{ number of replicates (observations) per group}
  \item{fixed}{vector of lenght 3 with mean, variance and estimate of fixed effect to simulate. default value: c(0,1,0)}
  \item{VI}{variance component of intercept (ID). Could be specified as a vector. default value: seq(0.05,0.95,0.05)}
  \item{VS}{variance component of slope (ID*fixed effect interaction). Could be specified as a vector. default value :seq(0.05,0.5,0.05))}
  \item{CoIS}{value of correlation or covariance between random intercept and random slope}
  \item{relIS}{"cor" or "cov" set the type of relation give in CoIS. By default the relation is set to correlation}
}
\details{
 P-values for random effects are estimated using a log-likelihood ratio
 test between two models with and without the effect. Power represent
 the percentage of simulations providing a significant p-value for a
 given random structure.  
Residual variance (e) is calculted as 1-VI.  
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
\seealso{ \code{\link{PAMM}},\code{\link{SSF},\code{\link{plot.EAMM}}}
}
\examples{
## Not run:
#  ours=EAMM(numsim=10,group=100,repl=4,fixed=c(0,1,1),VI=seq(0.05,0.3,0.05),
#             VS=seq(0.05,0.2,0.05))
#  plot(ours, "both")
## End(Not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
