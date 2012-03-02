\name{boot.var}
\alias{boot.var}
\title{
Bootstrap estimator for the variance of the Heteroscedastic Nonlinear Regression Growth model parameters
}
\description{
This function computes the estimated variance of the Heteroscedastic 
Nonlinear Regression Growth Model using Bootstrap estimator

}
\usage{
boot.var(data, betas, rho, sigma2, shape, nu, B, size, method = "EM", 
type = "N", model = "VB", m.type = "power", replace = TRUE, 
category = TRUE, prop = TRUE)
}
\arguments{
  \item{data}{a numeric matrix of \code{data.frame} type with labels \code{x} (ages) and \code{y} (lengths)}
  \item{betas}{a initial numeric vector of Growth model parameters}
  \item{rho}{a initial numeric value of mean residual parameter}
  \item{sigma2}{a initial numeric value of variance residual parameter}
  \item{shape}{a initial numeric value of shape residual parameter}
  \item{nu}{a initial numeric value of freedom degree residual parameter}
  \item{B}{a integer value of number of Bootstrap simulations}
  \item{size}{a integer value of \code{B} Bootstrap simulation sample size}
  \item{method}{a string to select the method of estimation. For maximum likelihood: \code{MLE}, for the \code{\link{nls}} funcion; 
  \code{EM}, for the \code{\link{HNL.skew}} funcion}
  \item{model}{a string related to growth model. It is available for:  \code{VB} (Von Bertalanffy, by default), 
 \code{Gompertz},  \code{Logistic} and \code{Richards}}
  \item{type}{a string related to residual distribution. It is available for:  \code{ST} (Skew-t, by default),  \code{SN} (Skew-Normal),
 \code{T} (T-Student) and \code{N} (Normal)}
  \item{m.type}{a string related to heteroscedastic residual modelation. It is available for:  \code{power} (power type,
by default) and  \code{exp} (exponencial type). See \code{\link{HNL.skew}} function for details}
  \item{replace}{a logical value to select \code{replace} in Bootstrap. By default is \code{TRUE} (with \code{replace})}
  \item{category}{a logical value to select ages class. By default is \code{TRUE} (use the class)}
  \item{prop}{a logical value to select proportion at age, i.e, maintaining the original proportion at when   
Bootstrapping. By default is \code{TRUE}}
}
\value{
\item{V}{a numeric matrix of estimated variance}
}
\references{
Wasserman, L. (2006). \emph{ All of Nonparametric Statistics}. Springer-Verlag, New York Inc, USA.

Efron, B., Tibshirani, R., J. (1993). \emph{ An Introduction to the Bootstrap}. Chapman & Hall/CRC, Monographs on Statistics & Applied Probability.
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{\link{jn.var}}, \code{\link{HNL.skew}}, \code{\link{nls}}
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

X = data.frame(list(x=merluzaChile$edad, y=merluzaChile$long))
p=round(nrow(X)/2,0)

# Run:
# bvar <-  boot.var(data = X, beta, rho, sigma2, shape, nu, B=50, size=p, 
#	   method = "EM", type = "N", model = "VB", m.type="power", 
#	   replace = TRUE, category = TRUE, prop = TRUE)
# bvar$V

modelVB <- HNL.skew(y, x, beta, rho, sigma2, shape, nu, 
		loglik = FALSE, model = "VB", type = "N", 
		m.type = "power", error = 0.00001)

sqrt(diag(modelVB$V))
}
\keyword{Variance}
\keyword{Estimation}
\keyword{skewtools}
