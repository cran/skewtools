\name{merluzaChile}
\alias{merluzaChile}
\docType{data}
\title{
Merluccius gayi gayi Data
}
\description{
Length at age database of Common Hake (Merluccius gayi gayi, males fished off Chilean coast during 2005)

}
\usage{data(merluzaChile)}
\format{
  A data frame containing 1638 observations on the following 2 variables.
  \describe{
    \item{\code{edad}}{a numeric vector of ages (years)}
    \item{\code{long}}{a numeric vector of lengths (cm)}
  }
}

\source{
Instituto de Fomento Pesquero, Valparaiso, Chile. URL \url{http://www.ifop.cl}
}
\examples{
data(merluzaChile)
str(merluzaChile)
dim(merluzaChile)

x <- merluzaChile$edad
y <- merluzaChile$long
plot(x, y, main="Age-Length", ylab="Length (cm)", xlab="Age (years)")
}
\keyword{datasets}
\keyword{skewtools}