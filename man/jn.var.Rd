\name{jn.var}
\alias{jn.var}

\title{
Jackknife estimator for the variance of the Heteroscedastic Nonlinear
Regression Growth model parameters
}
\description{
This function computes the estimated variance of the Heteroscedastic Nonlinear Regression Growth
Model using Jackknife and Weighted Jackknife described in Wu (1986)
}
\usage{
jn.var(data, betas, rho, sigma2, shape, nu, weight = FALSE, 
method = "MLE", model = "VB", type = "ST", m.type = "power")
}
\arguments{
   \item{data}{a numeric matrix of \code{data.frame} type with labels \code{x} (ages) and \code{y} (lengths)}
  \item{betas}{a initial numeric vector of Growth model parameters}
  \item{rho}{a initial numeric value of mean residual parameter}
  \item{sigma2}{a initial numeric value of variance residual parameter}
  \item{shape}{a initial numeric value of shape residual parameter}
  \item{nu}{a initial numeric value of freedom degree residual parameter}
  \item{weight}{a logical value to use the weighted Jackknife method. By default is \code{FALSE} (usual method)}
  \item{method}{a string to select the method of estimation. For maximum likelihood: \code{MLE}, for the \code{\link{nls}} funcion; 
  \code{EM}, for the \code{\link{HNL.skew}} funcion}
  \item{model}{a string related to growth model. It is available for:  \code{VB} (Von Bertalanffy, by default), 
 \code{Gompertz},  \code{Logistic} and \code{Richards}}
  \item{type}{a string related to residual distribution. It is available for:  \code{ST} (Skew-t, by default),  \code{SN} (Skew-Normal),
 \code{T} (T-Student) and \code{N} (Normal)}
  \item{m.type}{a string related to heteroscedastic residual modelation. It is available for:  \code{power} (power type,
by default) and  \code{exp} (exponencial type). See \code{\link{HNL.skew}} function for details}
}
\value{
\item{V}{a numeric matrix of estimated variance}
}
\references{
Wu, C., F., J. (1986). Jackknife, Bootstrap and Other Resampling Methods in Regression Growth Analysis. \emph{Annals of Statistics}, 14, 4, 1261-1295.
}
\author{
Javier E. Contreras-Reyes
}

\seealso{
\code{\link{boot.var}}, \code{\link{HNL.skew}}, \code{\link{nls}}
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
# bvar <-  jn.var(data = X, beta, rho, sigma2, shape, nu, weight=FALSE, 
#	   method = "MLE", model = "VB", type = "T", m.type="power")
# bvar$V

st = c(b1=beta[1], b2=beta[2], b3=beta[3])
ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3))), data=X, start=st)
sqrt(diag(vcov(ml)))

# Run:
# bvar2 <-  jn.var(data = X, beta, rho, sigma2, shape, nu, weight=FALSE, 
# 	   method = "MLE", model = "VB", type = "T", m.type="power")
# bvar2$V

modelVB <- HNL.skew(y, x, beta, rho, sigma2, shape, nu, 
		loglik = FALSE, model = "VB", type = "T", 
		m.type = "power", error = 0.00001)

sqrt(diag(modelVB$V))
}
\keyword{Variance}
\keyword{Estimation}
\keyword{skewtools}
