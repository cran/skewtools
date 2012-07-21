\name{Jsn.test}
\alias{Jsn.test}

\title{
Kupperman-test to comparing two Skew-normal samples
}
\description{
Performs two sample Kupperman-tests on vectors of data distributed univariate or/and multivariate Skew-normal
}
\usage{
Jsn.test(X, Y)
}

\arguments{
  \item{X}{a numerical vector or a numerical matrix with complete data corresponding to the first distribution}
  \item{Y}{a numerical vector or a numerical matrix with complete data corresponding to the second distribution}
}

\details{
This test is based on J-distance \eqn{J(X,Y)=D_{KL}(X,Y)+D_{KL}(Y,X)} and it is extractable from \code{STATISTIC} value using
\eqn{\frac{N_{1}N_{2}}{N_{1}+N_{2}} J(X,Y)}. If \eqn{dim(X)=dim(Y)=1}, \code{Jsn.test} consider \code{\link{KL.sn}} function.
}

\value{
\item{STATISTIC}{the value of the test statistic}
\item{p.value}{the p-value for the test}
}

\references{
Contreras-Reyes, J.E., Arellano-Valle, R.B. (2012). Kullback-Leibler divergence measure for Multivariate Skew-Normal Distributions. Pre-print.

Salicru, M., Menendez, M.L., Pardo, L., Morales, D. (1994). On the applications of divergence type measures in testing statistical hypothesis. \emph{ J. Multivariate Anal.}, 51, 372-391.

Kupperman, M. (1957). \emph{ Further applications of information theory to multivariate analysis and statistical Inference}. PhD dissertation.
}

\author{
Javier E. Contreras-Reyes
}

\seealso{
\code{\link{KL.sn}}
}

\examples{
data(Mag27F)
length(Mag27F)
data1=Mag27F[Mag27F <=6][1:100]
hist(data1)
data2=Mag27F
hist(data2)

Jsn.test(data1, data2)
}

\keyword{divergence}
\keyword{information}
\keyword{skewtools}
