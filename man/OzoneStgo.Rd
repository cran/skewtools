\name{OzoneStgo}
\alias{OzoneStgo}
\docType{data}
\title{
Ozone Time Series Data
}
\description{
Ozone Time Series Data set of hourly observations (Santiago de Chile)
}
\usage{data(OzoneStgo)}
\format{
  A data frame with 744 observations on the following 9 variables.
  \describe{
    \item{\code{date}}{a numeric vector of dates}
    \item{\code{hour}}{a numeric vector of hours}
    \item{\code{F}}{a numeric vector of \emph{Independencia} Station Ozone concentrations}
    \item{\code{L}}{a numeric vector of \emph{La Florida} Station Ozone concentrations}
    \item{\code{M}}{a numeric vector of \emph{Las Condes} Station Ozone concentrations}
    \item{\code{N}}{a numeric vector of \emph{Parque O'Higgins} Station Ozone concentrations}
    \item{\code{O}}{a numeric vector of \emph{Pudahuel} Station Ozone concentrations}
    \item{\code{P}}{a numeric vector of \emph{Cerrillos} Station Ozone concentrations}
    \item{\code{Q}}{a numeric vector of \emph{El Bosque} Station Ozone concentrations}
  }
}
\details{
This Dataset correspond to MACAM network in Seremi de Salud (March, 2006) of time series of ozone concentrations at 7 monitoring stations 
(\code{F}, \code{L}, \code{M}, \code{N}, \code{O}, \code{P}, \code{Q}) with \code{n=5208} hourly observations
}
\source{
Data base available at Seremi de Salud (2006), \emph{Indices de Calidad del Aire, Santiago de Chile},
\url{http://www.seremisaludrm.cl/sitio/pag/aire/indexjs3aireindgasesdemo-prueba.asp} 

}
\references{
Arellano-Valle, R. B., Contreras-Reyes, J. E., Genton, M. G. (2012). Shannon Entropy and Mutual Information for Multivariate
Skew-Elliptical Distributions. \emph{Scandinavian Journal of Statistics}, doi: 10.1111/j.1467-9469.2011.00774.x.

Silva, C., Quiroz, A. (2003). Optimization of the atmospheric pollution monitoring network at Santiago de Chile. 
\emph{Atmospheric Environment}, 37, 2337-2345.
}
\examples{
data(OzoneStgo)
str(OzoneStgo) 

data1=OzoneStgo[,3:9]

## Movile Average s-Function 
mmf <- function(X, s=1){
m=dim(X)[2]; t=dim(X)[1]
n=dim(X)[1]/s; datos.a=matrix(NA,t,m)
for(i in (s+1):t) for(j in 1:m){ 
aux=X[(i-s+1):i,j]
datos.a[i,j]=mean(aux[is.na(aux)==FALSE]) 
}
datos.a[(s+1):t,]
}

est=2 # L Station
dat0=mmf(data1,s=8)[,est]
dat1=mmf(data1,s=16)[,est]
dat2=mmf(data1,s=24)[,est]
dat3=mmf(data1,s=32)[,est]

plot(data1[1:(24*7),est],type="l",lty=2,
ylab=expression(paste(O[3] (mu ,g/m^3))), 
xlab="time (hours)", main="Ozone - Station L")
lines(dat0[1:(24*7)],type="l",lwd=2,col="yellow") 
lines(dat1[1:(24*7)],type="l",lwd=2,col="red") 
lines(dat2[1:(24*7)],type="l",lwd=2,col="green") 
lines(dat3[1:(24*7)],type="l",lwd=2,col="blue") 
legend(122,90,c("Original",expression(MA[8]),expression(MA[16]),
expression(MA[24]),expression(MA[32])), 
lty=c(2,1,1,1,1),lwd=c(1,2,2,2,2),
col=c("black","yellow","red","green","blue"))
}
\keyword{datasets}
\keyword{skewtools}
