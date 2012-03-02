\name{Mag27F}
\alias{Mag27F}
\docType{data}
\title{
Maule 2010 Earthquake Data
}
\description{
Dataset of Local Richter Magnitudes of 27 Fabruary 2010 Earthquake (Maule, Chile)
}
\usage{data(Mag27F)}
\format{
  A numeric vector containing 629 observations. The format is:
 num [1:629] 2.3 4.1 2.5 8.8 6 5.6 6 5.2 4.9 4.9 ...
}
\source{
Servicio Sismologico Nacional, Universidad de Chile. URL \url{http://ssn.dgf.uchile.cl/}
}
\references{
Lange, D., Tilmann, F., Barrientos, S., Contreras-Reyes, E., Methe, P., Moreno, M., Heit, B.,
Bernard, P., Vilotte, J., Beck, S. (2012). Aftershock seismicity of the 27 February 2010 Mw 8.8
Maule earthquake rupture zone. \emph{Earth and Planetary Science Letters}, doi:10.1016/j.epsl.2011.11.034.
}
\examples{
data(Mag27F)
str(Mag27F)
require(sn)

hist(Mag27F, main="27 February 2010 Mw 8.8 Maule earthquake")
msn.fit(y=Mag27F)
mst.fit(y=Mag27F)
}
\keyword{datasets}
\keyword{skewtools}