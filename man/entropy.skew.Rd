\name{entropy.skew}
\alias{entropy.skew}

\title{
Shannon Entropy for Skewed Distributions
}
\description{Compute the Shannon Entropy for Skew-Normal and Skew-t Distributions}
\usage{
entropy.skew(X, family = "ST")
}
\arguments{
  \item{X}{a numerical vector or a numerical matrix with complete data}
  \item{family}{string to select the skewed distribution. Its available for: \code{ST} (Skew-t, by default) and \code{SN} (Skew-Normal)}
}
\details{
Shannon Entropy and Mutual information Index for Skewed Distributions are provided in Arellano-Valle et al. (2012). 
Estimation of parameters are provided by the functions \code{\link{mst.fit}} (Skew-t case) and \code{\link{msn.fit}} 
(Skew-Normal case) from \code{sn} package
}
\value{
\item{H}{a numerical value of the entropy}
\item{family}{a string of the used distribution}
}
\references{
Arellano-Valle, R. B., Contreras-Reyes, J. E., Genton, M. G. (2012). Shannon Entropy and Mutual Information for Multivariate
Skew-Elliptical Distributions. \emph{Scand. J. Stat.}, doi: 10.1111/j.1467-9469.2011.00774.x.

Javier, W. R., Gupta, A. K. (2009). Mutual information for certain multivariate distributions. \emph{Far East
J. Theor. Stat.}, 29, 39-51.
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{\link{mst.fit}}, \code{\link{msn.fit}}
}
\examples{
data(OzoneStgo)
data1=OzoneStgo[,3:9]

entropy.skew(X=data1[1:100,], family = "ST")
entropy.skew(X=data1[1:100,], family = "SN")

## Mutual Information:

Hx = entropy.skew(data1[1:100,1:6],"ST")$H 
Hy = entropy.skew(data1[1:100,7],"ST")$H 
Hxy = entropy.skew(data1[1:100,1:7],"ST")$H
Hx + Hy - Hxy

Hx = entropy.skew(data1[1:100,-5],"SN")$H 
Hy = entropy.skew(data1[1:100,5],"SN")$H 
Hxy = entropy.skew(data1[1:100,],"SN")$H
Hx + Hy - Hxy

Hx = entropy.skew(data1[1:100,-5],"ST")$H 
Hy = entropy.skew(data1[1:100,5],"ST")$H 
Hxy = entropy.skew(data1[1:100,],"ST")$H
Hx + Hy - Hxy
}
\keyword{entropy}
\keyword{information}
\keyword{skewtools}
