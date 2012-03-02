\name{HNL.skew}
\alias{HNL.skew}
\title{
Expected Conditional Maximization Estimation of Heteroscedastic Nonlinear Regression Growth model parameters
}
\description{
Parameters estimation of the von Bertalanffy, Gompertz, Logistic and Richards growth models us-
ing the robust Expected Conditional Maximization Estimation (ECME) method in heteroscedastic
Nonlinear Regression Growth models on scale mixtures of skew-normal distributions. 
The distributions considered are: : Normal, T-Student, Skew-Normal and Skew-t
}
\usage{
HNL.skew(y, x, betas, rho, sigma2, shape, nu, loglik = FALSE, 
model = "VB", type = "ST", m.type = "power", error = 1e-05)
}
\arguments{
  \item{y}{a numeric vector of lengths}
  \item{x}{a numeric vector of ages}
  \item{betas}{a initial numeric vector of Growth model parameters}
  \item{rho}{a initial numeric value of mean residual parameter}
  \item{sigma2}{a initial numeric value of variance residual parameter}
  \item{shape}{a initial numeric value of shape residual parameter}
  \item{nu}{a initial numeric value of freedom degree residual parameter}
  \item{loglik}{a logical value to obtain the log-likelihood. By default is \code{FALSE}}
  \item{model}{a string related to growth model. It is available for: \code{VB} (Von Bertalanffy, by default), 
 \code{Gompertz},  \code{Logistic} and \code{Richards}}
  \item{type}{a string related to residual distribution. It is available for:  \code{ST} (Skew-t, by default),  \code{SN} (Skew-Normal),
 \code{T} (T-Student) and \code{N} (Normal)}
  \item{m.type}{a string related to heteroscedastic residual modelation. It is available for:  \code{power} (power type,
by default) and  \code{exp} (exponencial type). See below for details}
  \item{error}{a numeric value for convergence tolerance of ECME method. By default is \code{1e-05}}
}
\details{
This function is based on the ECME method of the scripts written by Lachos et al. (2012) 
(open access: \url{http://www.ime.unicamp.br/~hlachos/codeSMSN-HNLM.rar})

The \code{HNL.skew} function consider a Heteroscedastic Nonlineal Regression Growth model 
\eqn{L_t=\eta(x_t,\beta)+\varepsilon_t} and, are available the following growth models:


von Bertalanffy (von Bertalanffy, 1941)

\deqn{L_t=L_{\infty}(1-e^{-k(t-t_0)}),}

Gompertz (Gompertz, 1825)

\deqn{L_t=L_{\infty}(e^{-e^{-k(t-t_0)}}),}

Logistic (Verhulst, 1838)

\deqn{L_t=\frac{L_{\infty}}{1+e^{-k(t-t_0)}},}

and Richards (Richards, 1959)

\deqn{L_t=L_{\infty}(1-e^{-k(t-t_0)})^M}

The residuals of the model are represented by the scale mixture  of skew-normal, where the first and second moment depends
on \eqn{m(x_t)} function. The \code{m.type} option allow to select the power function \eqn{m(x_t)=x_{t}^{\rho}} through
\code{power} or, the exponential function \eqn{m(x_t)=e^{\rho x_{t}}} using \code{exp}.

An incorrect input of the initial values for the parameters may induce an error in convergence. Therefore, is recommended 
to assess visually how accurate the initial parameters are by plotting the chosen growth function model and the data. 
In the example (see below) it is shown the initial values  of \code{beta}, \code{rho}, \code{sigma2}, \code{lambda} 
and \code{nu}. The ECME algorithm need, at least, a good vector of \code{beta}; the residuals parameters shown in 
example, should be used for estimation of others species
}
\value{
\item{y}{a numeric vector of lengths used in the estimation}
\item{x}{a numeric vector of ages used in the estimation}
\item{betas}{a numeric vector of estimated parameters of the growth model by ECME method}
\item{rho}{a numeric vector of mean residual parameter of the growth model by ECME method}
\item{sigma2}{a numeric vector of variance residual parameter of the growth model by ECME method}
\item{shape}{a numeric vector of shape residual parameter of the growth model by ECME method}
\item{loglik}{a numeric value of log-likelihood of the obtained model}
\item{iter}{a numeric value of total iterations made by ECME method}
\item{res}{a numeric vector of model residuals}
\item{n}{a numeric value of sample size}
\item{V}{a numeric matrix of variance-covariance of the growth model parameters}
\item{model}{a call of the fitted model of type \code{\link{nls}}}
\item{curve}{a string related to the used growth model}
\item{distr}{a string related to the used residual distribution}
\item{m.type}{a string related to the used heteroscedastic residual modeling}
}
\references{
Basso, R. M., Lachos, V. H., Barbosa, C. R., Ghosh, P., (2010). Robust mixture modeling based on scale mixtures of skew-normal distributions. \emph{Comput. Stat. Data An.}, { 54}, { 12}, 2926-2941.

Gompertz, B. (1825). On the Nature of the Function Expressive of the Law of Human Mortality, and on a New Mode of Determining the Value of Life Contingencies. 
\emph{Philosophical Transactions of the Royal Society of London}, 115, 513-585.

Lachos, V. H., Bandyopadhyay, D., Garay,  A. M., (2011). Heteroscedastic nonlinear Regression Growth models based on scale mixtures of skew-normal distributions. \emph{ Stat. Probabil. Lett.}, { 81}, { 8}, 1208-1217.

Lachos, V. H., Garay, A. M., Ortega, E. M., Vilca, L. F. (2012). Estimation and diagnostics for heteroscedastic nonlinear Regression Growth models based on scale mixtures of skew-normal distributions. \emph{J. Stat. Plann. Inf.}, to appear.

Schnute, J. T., Richards, L. J. (1990). A Unified Approach to the Analysis of Fish Growth, Maturity, and Survivorship Data. \emph{Can. J. Fish. Aquat. Sci.}, 47, 24-37. 

Verhulst, P. -F. (1838). Notice sur la loi que la population poursuit dans son accroissement. \emph{Correspondance mathematique et physique}, 10, 113-121.

von Bertalanffy, L. (1938).Aquantitative theory of organic growth (Inquiries on growth laws. II). Human Biology 10, 181-213.
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{\link{nls}}
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

modelVB <- HNL.skew(y, x, beta, rho, sigma2, shape, nu, loglik = FALSE, 
model = "VB", type = "T", m.type = "power", error = 0.00001)
modelVB
}

\keyword{Estimation}
\keyword{skewtools}
