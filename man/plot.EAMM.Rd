\name{plot.EAMM}
\alias{plot.EAMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graphic output of the EAMM function }
\description{
  provide graphic interpretation of the simulation results
}
\usage{
\method{plot}{EAMM}(x, graphtype = "both", vi, vs,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an EAMM object }
  \item{graphtype}{ "VI", "VS","both" or "both.dyn"
		"VI" give graphs with varying variance component of intercept and with a fixed variance component for slope (should be specified in vs)
              	"VS" give graphs with varying variance component for slope and with a fixed variance component of intercept (should be specified in vi)
  		"both" 3-D plot using lattice package
	        "both.dyn" 3D plot using rgl package}
  \item{vi}{ VI for which plots the output (necessary for "VS" type of graph) }
  \item{vs}{ VS for which plots the output (necessary for "VI" type of graph) }
  \item{\dots}{ potentially further arguments to pass to methods }
}
\details{ }
\value{ }
\references{ }
\author{ Julien Martin }
\note{ }
\seealso{\code{\link{EAMM}}, \code{\link{plot.PAMM}}, \code{\link{plot.SSF}} }
\examples{
   \dontrun{
  ours=EAMM(numsim=10,group=10,repl=4,VI=seq(0.05,0.95,0.05),VS=c(0.05,0.1))
  plot(ours, "both")
  plot(ours, "VI",vi=0.4)
  plot(ours,"VS",vs=0.1) 
   } 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc }% __ONLY ONE__ keyword per line
