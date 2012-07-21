\name{infl.growth}
\alias{infl.growth}

\title{
Local Influence Analysis
}
\description{
Diagnostic the influential observations of Heteroscedastic Nonlinear Regression Growth models using Local Influence Analysis of Cook
with curvature correction
}
\usage{
infl.growth(y, x, model, plot.it = TRUE)
}
\arguments{
\item{y}{a numeric vector of lengths}
  \item{x}{a numeric vector of ages}
  \item{model}{a object related to fitted model of type \code{\link{nls}}}
  \item{plot.it}{a logical value to plot. By default, is \code{TRUE}}
}
\value{
\item{residuals}{a numeric vector of model residuals}
\item{outliers}{a numeric vector with positions of outliers}
}
\references{
Cook, R. D., (1986). Assessment of local influence (with discussion). \emph{J. R. Stat. Soc. Ser. B}, 48, 2, 133-169.

Poon, W., Poon, Y. S., (1999). Conformal normal curvature and assessment of local influence. \emph{ J. R. Stat. Soc. Ser. B}, 61, 1, 51-61.
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{\link{conf.im}}, \code{\link{HNL.skew}}
}
\examples{
data(merluzaChile)
x <- merluzaChile$edad
y <- merluzaChile$long
plot(x, y, main="Age-Length", ylab="Length (cm)", xlab="Age (years)")

beta <- c(80, 0.08, -0.187)
rho <- -0.1
sigma2 <- 3.2726
shape <- 0.1698
nu <- 11

modelVB1 <- HNL.skew(y, x, beta, rho, sigma2, shape, nu, loglik = TRUE, 
model = "VB", type = "T", m.type = "power", error = 0.00001)
modelVB2 <- HNL.skew(y, x, beta, rho, sigma2, shape, nu, loglik = TRUE, 
model = "VB", type = "N", m.type = "power", error = 0.00001)

I1=infl.growth(y, x, model=modelVB1, plot.it = TRUE)
I2=infl.growth(y, x, model=modelVB2, plot.it = TRUE)

plot(x, y, main="Age-Length", ylab="Length (cm)", xlab="Age (years)")
points(x[I1$outliers],y[I1$outliers],col="red")
legend(10,30,c("T"),col=c("red"),lty=1)
points(x[I2$outliers],y[I2$outliers],col="yellow")
legend(10,20,c("N: tau = 3"),col=c("yellow"),lty=1)
}
\keyword{Diagnostic}
\keyword{skewtools}
