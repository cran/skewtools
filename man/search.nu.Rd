\name{search.nu}
\alias{search.nu}

\title{
Profile Maximum Log-likelihood based on heavy tailed distributions
}
\description{
Computes the maximum log-likelihood profiles for a given nu value. 
}
\usage{
search.nu(y, x, betas, rho, sigma2, shape, nu = c(2:15), 
model = "VB", type = "ST", m.type = "power", plot.it = TRUE)
}

\arguments{
\item{y}{a numeric vector of lengths}
  \item{x}{a numeric vector of ages}
  \item{betas}{a initial numeric vector of Growth model parameters}
  \item{rho}{a initial numeric value of mean residual parameter}
  \item{sigma2}{a initial numeric value of variance residual parameter}
  \item{shape}{a initial numeric value of shape residual parameter}
  \item{nu}{a initial numeric value of freedom degree residual parameter}
  \item{model}{a string related to growth model. It is available for:  \code{VB} (Von Bertalanffy, by default), 
 \code{Gompertz},  \code{Logistic} and \code{Richards}}
  \item{type}{a string related to residual distribution. It is available for:  \code{ST} (Skew-t, by default),  \code{SN} (Skew-Normal),
 \code{T} (T-Student) and \code{N} (Normal)}
  \item{m.type}{a string related to heteroscedastic residual modelation. It is available for:  \code{power} (power type,
by default) and  \code{exp} (exponencial type). See below for details}
  \item{plot.it}{a logical value to plot the values \code{nu} vs \code{loglik}. By default is \code{TRUE}}
}
\value{
\item{table.res}{a numeric matrix with results. Include the values \code{betas}, \code{rho}, \code{sigma2}, \code{shape} and \code{loglik} for each \code{nu}}
\item{pos.nu}{a numeric value of \code{nu} associated to maximum likelihood}
}
\references{
Lange, K. L., Sinsheimer, J. S., (1993). Normal/independent distributions and their applications in robust Regression Growth. \emph{J. Comput. Graph. Stat}, 2, 2, 175-198.
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

S1 <- search.nu(y, x, beta, rho, sigma2, shape, 
nu = c(25:30), model = "VB", type = "T", m.type="power", 
plot.it = TRUE)
S1$table.res
}
\keyword{Estimation}
\keyword{skewtools}
