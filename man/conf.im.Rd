\name{conf.im}
\alias{conf.im}
\title{
Confidence Intervals for Heteroscedastic Nonlinear Regression Growth models  
}
\description{
Plot the Heteroscedastic nonlinear Regression Growth model and its confidence intervals
}
\usage{
conf.im(x, y, model, alpha = 0.05, plot.it = TRUE)
}
\arguments{
  \item{y}{a numeric vector of lengths}
  \item{x}{a numeric vector of ages}
  \item{model}{a object related to fitted model of type \code{\link{nls}}}
  \item{alpha}{confidence level of the interval}
  \item{plot.it}{a logical value to plot. By default, is \code{TRUE}}
}

\value{
\item{yeval}{a numeric vector about the evaluation of \code{x} on growth model}
\item{ICS}{a numeric vector related to upper band of estimation confidence intervals}
\item{ICI}{a numeric vector related to low band of estimation confidence intervals}
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{\link{HNL.skew}}
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

conf.im(x, y, model=modelVB1, alpha = 0.05, plot.it = TRUE)
conf.im(x, y, model=modelVB2, alpha = 0.05, plot.it = TRUE)
}
\keyword{Diagnostic}
\keyword{skewtools}
