\name{plot.SSF}
\alias{plot.SSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graphic output of the PAMM function }
\description{
  provide graphic interpretation of the simulation results
}
\usage{
\method{plot}{SSF}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an SSF object }
  \item{\dots}{ potentially further arguments to pass to methods }
}
\details{ }
\value{ }
\references{ }
\author{ Julien Martin }
\note{ }
\seealso{\code{\link{SSF}}, \code{\link{plot.PAMM}}, \code{\link{plot.EAMM}} }
\examples{
   \dontrun{
   oursSSF=SSF(50,100,10,c(0.4,0.1,0.6,0))
   plot(oursSSF)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
\keyword{ misc }% __ONLY ONE__ keyword per line
