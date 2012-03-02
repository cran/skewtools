\name{AICgrowth}
\alias{AICgrowth}
\title{
Information Criterion for Heteroscedastic Nonlinear Regression Growth models
 
}
\description{
This function computes AIC, BIC and AICc (corrected AIC) information criteria for the Heteroscedastic nonlinear Regression Growth models
}
\usage{
AICgrowth(model)
}

\arguments{
  \item{model}{a fitted model object obtained by the function \code{HNL.skew}}
}
\details{
The Akaike Information Criterion (AIC) value is computed as follow: 

\deqn{AIC(\widehat{\theta})=-2(L(\widehat{\theta})-k).}

where \eqn{theta} is a parameter vector and \eqn{k} is the number of parameters. The corrected AIC value is

\deqn{AICc(\widehat{\theta})=AIC(\widehat{\theta})+\frac{2k(k+1)}{n-k-1}.}

The Schwarz Information criterion (BIC) is

\deqn{BIC(\widehat{\theta})=-2L(\widehat{\theta})+k\log(n).}
}
\value{
\item{n}{a numeric value of the sample size}
\item{k}{a numeric value of the parameters number of the model}
\item{L}{a numeric value with the corresponding log-likelihood of the model}
\item{AIC}{a numeric value with the corresponding AIC}
\item{AICc}{a numeric value with the corresponding corrected AIC}
\item{BIC}{a numeric value with the corresponding BIC}
}
\references{
Akaike, H. (1974). A new look at the statistical model identification. \emph{IEEE Transactions on Automatic Control}, 19, 6, 716-723.

Burnham, K. P., Anderson, D. R. (2002). \emph{Model Selection and Multimodel Inference, A Practical Information-Theoretic Approach}, 2nd ed., Springer-Verlag.

Schwarz, G. E. (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, 6, 2, 461-464.

Hurvich, C. M., Tsai, C. (1989). Regression Growth and time series model selection in small samples. \emph{Biometrika}, 76, 297-307.
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

# If loglik = FALSE, AICm not run

AICgrowth(modelVB1)
AICgrowth(modelVB2)
}

\keyword{Diagnostic}
\keyword{skewtools}
