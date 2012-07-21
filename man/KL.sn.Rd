\name{KL.sn}
\alias{KL.sn}

\title{
Kullback-Leibler divergence for univariate Skew-normal distributions
}
\description{
Compute the Kullback-Leibler divergence measure between two univariate Skew-Normal distributions 
given by its respective parameters
}
\usage{
KL.sn(beta1, beta2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta1}{a numerical vector of length = 3 with location, dispertion and shape parameters of the first distribution}
  \item{beta2}{a numerical vector of length = 3 with location, dispertion and shape parameters of the second distribution}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{Dkl}{a numeric value corresponding to Kullback-Leibler divergence}
\item{max.entropy}{a numeric value corresponding to upper bound of Shannon entropy given by maximum entropy principle}
}

\references{
Contreras-Reyes, J.E., Arellano-Valle, R.B. (2012). Kullback-Leibler divergence measure for Multivariate
Skew-Normal Distributions. Pre-print.

Arellano-Valle, R.B., Contreras-Reyes, J.E., Genton, M.G. (2012). Shannon Entropy and Mutual Information for Multivariate
Skew-Elliptical Distributions. \emph{Scand. J. Stat.}, doi: 10.1111/j.1467-9469.2011.00774.x.
}
\author{
Javier E. Contreras-Reyes
}

\seealso{
\code{\link{Jsn.test}}, \code{\link{entropy.skew}} 
}
\examples{
data(Mag27F)
length(Mag27F)
data1=Mag27F[Mag27F <=6][1:100]
hist(data1)
data2=Mag27F
hist(data2)

modelX=msn.fit(y=data1, plot.it=FALSE)$dp
modelY=msn.fit(y=data2, plot.it=FALSE)$dp
beta1=c(modelX$beta, modelX$Omega, modelX$alpha)
beta2=c(modelY$beta, modelY$Omega, modelY$alpha)

KL.sn(beta1, beta2)

entropy.skew(X=data1, family = "SN")
entropy.skew(X=data2, family = "SN")
}

\keyword{divergence}
\keyword{information}
\keyword{skewtools}
