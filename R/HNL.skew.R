HNL.skew <-
function(y, x, betas, rho, sigma2, shape, nu, loglik = FALSE, model = "VB", type = "ST", m.type = "power", error = 0.00001){

z=x
options(warn=-1)

d.rho <-
function(x, rho, type="power", der=0)
{
if(type=="power") {
if(der==0) d=x^rho
if(der==1) d=x^rho * log(x)
if(der==2) d=x^rho * log(x)^2
}
if(type=="exp") { 
if(der==0) d=exp(x * rho)
if(der==1) d=exp(x * rho) * x
if(der==2) d=exp(x * rho) * x * x
}
return(d)
}

CGM <-
function(betas, x, type="VB", der=0)
{
n <- length(x)
if(type=="Logistic") {
if(der==0) d <- betas[1]/(1+exp(-betas[2]*(x-betas[3])))
if(der==1) {
db1 <- 1/(1 + exp(-betas[2] * (x - betas[3])))
db2 <- betas[1] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]))/(1 + exp(-betas[2] * (x - betas[3])))^2
db3 <- -(betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2])/(1 + exp(-betas[2] * (x - betas[3])))^2)
}
if(der==2) {
db1b1 <- 0
db2b2 <- -(betas[1] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]) * (x - betas[3]))/(1 + exp(-betas[2] *     (x - betas[3])))^2 - betas[1] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3])) * (2 *     (exp(-betas[2] * (x - betas[3])) * (x - betas[3]) * (1 + exp(-betas[2] * (x - betas[3])))))/((1 + exp(-betas[2] * (x - betas[3])))^2)^2)
db3b3 <- -(betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2] * betas[2])/(1 + exp(-betas[2] * (x - betas[3])))^2 -     betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2]) * (2 * (exp(-betas[2] * (x - betas[3])) * betas[2] * (1 + exp(-betas[2] * (x - betas[3])))))/((1 + exp(-betas[2] * (x - betas[3])))^2)^2)
db1b2 <- exp(-betas[2] * (x - betas[3])) * (x - betas[3])/(1 + exp(-betas[2] * (x - betas[3])))^2
db1b3 <- -(exp(-betas[2] * (x - betas[3])) * betas[2]/(1 + exp(-betas[2] * (x - betas[3])))^2)
db2b3 <- betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2] * (x - betas[3]) - exp(-betas[2] * (x - betas[3])))/(1 +     exp(-betas[2] * (x - betas[3])))^2 - betas[1] * (exp(-betas[2] * (x - betas[3])) * (x -     betas[3])) * (2 * (exp(-betas[2] * (x - betas[3])) * betas[2] * (1 + exp(-betas[2] * (x -     betas[3])))))/((1 + exp(-betas[2] * (x - betas[3])))^2)^2
}
}

if(type=="VB") {
if(der==0) d <- betas[1]*(1-exp(-betas[2]*(x-betas[3])))
if(der==1) {
db1 <- (1 - exp(-betas[2] * (x - betas[3])))
db2 <- betas[1] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]))
db3 <- -(betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2]))
}
if(der==2) {
db1b1 <- 0
db2b2 <- -(betas[1] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]) * (x - betas[3])))
db3b3 <- -(betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2] * betas[2]))
db1b2 <- exp(-betas[2] * (x - betas[3])) * (x - betas[3])
db1b3 <- -(exp(-betas[2] * (x - betas[3])) * betas[2])
db2b3 <- betas[1] * (exp(-betas[2] * (x - betas[3])) * betas[2] * (x - betas[3]) - exp(-betas[2] * (x - betas[3])))
}
}

if(type=="Gompertz") {
if(der==0) d <- betas[1]*(exp(-exp(-betas[2]*(x-betas[3]))))
		if(der==1) {
		db1 <- (exp(-exp(-betas[2] * (x - betas[3]))))
		db2 <- betas[1] * (exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * (x -     betas[3])))
		db3 <- -(betas[1] * (exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * betas[2])))
		}
		if(der==2) {
		db1b1 <- 0
		db2b2 <- betas[1] * (exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * (x -     betas[3])) * (exp(-betas[2] * (x - betas[3])) * (x - betas[3])) - exp(-exp(-betas[2] *     (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]) * (x - betas[3])))
		db3b3 <- -(betas[1] * (exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * betas[2] *     betas[2]) - exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) *     betas[2]) * (exp(-betas[2] * (x - betas[3])) * betas[2])))
	      db1b2 <- exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]))
     		db1b3 <- -(exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * betas[2]))
     		db2b3 <- betas[1] * (exp(-exp(-betas[2] * (x - betas[3]))) * (exp(-betas[2] * (x - betas[3])) * betas[2] *     (x - betas[3]) - exp(-betas[2] * (x - betas[3]))) - exp(-exp(-betas[2] * (x - betas[3]))) *     (exp(-betas[2] * (x - betas[3])) * betas[2]) * (exp(-betas[2] * (x - betas[3])) * (x -     betas[3])))
		}
	}

	if(type=="Richards") {
		if(der==0) d <- betas[1]*(1-exp(-betas[2]*(x-betas[3])))^betas[4]
		if(der==1) {
		db1 <- (1 - exp(-betas[2] * (x - betas[3])))^betas[4]
		db2 <- betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (betas[4] * (exp(-betas[2] * (x -     betas[3])) * (x - betas[3]))))
		db3 <- -(betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (betas[4] * (exp(-betas[2] *     (x - betas[3])) * betas[2]))))
		db4 <- betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^betas[4] * log((1 - exp(-betas[2] * (x - betas[3])))))
		}
		if(der==2) {
		db1b1 <- 0
		db2b2 <- betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^((betas[4] - 1) - 1) * ((betas[4] - 1) *     (exp(-betas[2] * (x - betas[3])) * (x - betas[3]))) * (betas[4] * (exp(-betas[2] * (x -     betas[3])) * (x - betas[3]))) - (1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) *     (betas[4] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]) * (x - betas[3]))))
		db3b3 <- -(betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (betas[4] * (exp(-betas[2] *     (x - betas[3])) * betas[2] * betas[2])) - (1 - exp(-betas[2] * (x - betas[3])))^((betas[4] -     1) - 1) * ((betas[4] - 1) * (exp(-betas[2] * (x - betas[3])) * betas[2])) * (betas[4] *     (exp(-betas[2] * (x - betas[3])) * betas[2]))))
		db4b4 <- betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^betas[4] * log((1 - exp(-betas[2] * (x -     betas[3])))) * log((1 - exp(-betas[2] * (x - betas[3])))))
	      db1b2 <- (1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (betas[4] * (exp(-betas[2] * (x - betas[3])) *     (x - betas[3])))
     		db1b3 <- -((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (betas[4] * (exp(-betas[2] * (x -     betas[3])) * betas[2])))
     		db1b4 <- (1 - exp(-betas[2] * (x - betas[3])))^betas[4] * log((1 - exp(-betas[2] * (x - betas[3]))))
		db2b3 <- betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (betas[4] * (exp(-betas[2] * (x -     betas[3])) * betas[2] * (x - betas[3]) - exp(-betas[2] * (x - betas[3])))) - (1 - exp(-betas[2] *     (x - betas[3])))^((betas[4] - 1) - 1) * ((betas[4] - 1) * (exp(-betas[2] * (x - betas[3])) *     betas[2])) * (betas[4] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]))))
     		db2b4 <- betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * log((1 - exp(-betas[2] *     (x - betas[3])))) * (betas[4] * (exp(-betas[2] * (x - betas[3])) * (x - betas[3]))) + (1 -     exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * (exp(-betas[2] * (x - betas[3])) * (x -     betas[3])))
     		db3b4 <- -(betas[1] * ((1 - exp(-betas[2] * (x - betas[3])))^(betas[4] - 1) * log((1 - exp(-betas[2] *     (x - betas[3])))) * (betas[4] * (exp(-betas[2] * (x - betas[3])) * betas[2])) + (1 - exp(-betas[2] *     (x - betas[3])))^(betas[4] - 1) * (exp(-betas[2] * (x - betas[3])) * betas[2])))
		}
	}

	if(der==0) L=d
	if(der==1) {
		if(type=="Logistic" | type=="VB" | type=="Gompertz") deriv1 <- cbind(db1,db2,db3)
		if(type=="Richards") deriv1 <- cbind(db1,db2,db3,db4)
		L=deriv1
	} 
	if(der==2) {
		if(type=="Logistic" | type=="VB" | type=="Gompertz") {
			deriv2 <- matrix(NA,3,3)		
			diag(deriv2) <- c(db1b1, db2b2, db3b3)
			deriv2[1,2] <- deriv2[2,1] <- db1b2; deriv2[1,3] <- deriv2[3,1] <- db1b3; deriv2[2,3] <- deriv2[3,2] <- db2b3 
		}
		if(type=="Richards" | type=="GVB") {
			deriv2 <- matrix(NA,4,4)		
			diag(deriv2) <- c(db1b1, db2b2, db3b3, db4b4)
			deriv2[1,2] <- deriv2[2,1] <- db1b2; deriv2[1,3] <- deriv2[3,1] <- db1b3; deriv2[1,4] <- deriv2[4,1] <- db1b4 
			deriv2[2,3] <- deriv2[3,2] <- db2b3; deriv2[2,4] <- deriv2[4,2] <- db2b4; deriv2[3,4] <- deriv2[4,3] <- db3b4
		}
		L=deriv2
	}

	return(L)	
}

dSN <- function(y, mu , sigma2 = 1, shape=1){
  dens <- 2*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*((y - mu)/sqrt(sigma2)))
  return(dens)
}

dSNH <- function(rho, y, z, mu, sigma2 = 1, shape=1, type){
  sigma2i<-sigma2*d.rho(z, rho, type, der=0)
  dens <- 2*dnorm(y, mu, sqrt(sigma2i))*pnorm(shape*((y - mu)/sqrt(sigma2i)))
  return(dens)
}


dt.ls <- function(x, loc , sigma2 = 1,shape=1, nu = 4){
  d <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}


dtH.ls <- function(rho, y, z, loc, sigma2 = 1,shape=1, nu = 4, type){

  sigma2i <- sigma2*d.rho(z, rho, type, der=0)
  d <- (y - loc)/sqrt(sigma2i)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2i)
  return(dens)
}

dSNC <- function(y, mu, sigma2, shape, nu){
    dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-1/2)*(y-mu)))
    return(dens)
  }

dSNCH <- function(rho, y, z, mu, sigma2, shape, nu, type){
    sigma2i<-sigma2*d.rho(z, rho, type, der=0)
    dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2i/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2i^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2i))*pnorm(shape*sigma2i^(-1/2)*(y-mu)))
    return(dens)
  }


dSS <- function(y, mu, sigma2, shape,nu){
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu[i],sqrt(sigma2[i]/u))*pnorm(u^(1/2)*shape*(sigma2[i]^(-1/2))*(y[i]-mu[i]))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}


dSSH <- function(rho, y, z, mu, sigma2, shape, nu, type){
  sigma2i<-sigma2*d.rho(z, rho, type, der=0)
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu[i],sqrt(sigma2i[i]/u))*pnorm(u^(1/2)*shape*(sigma2i[i]^(-1/2))*(y[i]-mu[i]))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}



rmix <- function(n, p1, rF1, arg1=NULL, arg2=NULL) {
	x1 <- vector(mode = "numeric", length = n)
    for (i in 1:n){
      u <- runif(1)
      if (u < p1) x1[i] <- do.call("rF1", c(list(1), arg1))
      if (u > p1) x1[i] <- do.call("rF1", c(list(1), arg2))
    }
    return(x1)
    }

geraSkewNormal <- function(mu, sigma2, shape){
   n=length(sigma2)
    delta <- shape / sqrt(1 + shape^2)
    y <- mu+ sqrt(sigma2)*(delta*abs(rnorm(n)) + (1 - delta^2)^(1/2)*rnorm(n))
    return(y)
  }

geraSNIt <- function(mu, sigma2, shape, nu ){
  n=length(sigma2)
  y <- mu + (rgamma(n, nu/2, nu/2))^(-1/2)*geraSkewNormal(0, sigma2, shape)
}


geraSNC <- function(mu, sigma2, shape, nu, gama){
   n=length(sigma2)
    y <- vector(mode = "numeric", length = n)
    for(i in 1:n){
    y[i]<-rmix(1, nu, geraSkewNormal, c(mu,sigma2[i]/gama,shape), c(mu,sigma2[i],shape))
                  }
                  return(y)
}

geraSS <- function(mu, sigma2, shape, nu){
   n=length(sigma2)
    u2 <- rbeta(n,nu,1)   # formula 10 do artigo e método da inversão
    ys <- mu + (u2)^(-1/2)*geraSkewNormal(0, sigma2, shape)
    return(ys)
}

##### Algoritmo EMCM


  if (type == "T"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
      k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-d.rho(z, rho, type=m.type, der=0)
      teta <- c(betas,rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old <- rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error){
      count <- count + 1
      ####print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

	    mu.1<-CGM(betas, x, type=model, der=0)
          mu<-mu.1+b*Delta*sqrt(wi)

          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / (Mtij)

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)

          S1<- u/wi
          S2  <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S22 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###

	    m1 <- nls(ymod ~ CGM(betas, x, type=model, der=0), start=list(betas=betas.old),trace=TRUE,weights=S1)
          betas <- as.vector(coef(m1)) #sum(S1*y - Delta.old*S2) / sum(S1)
	    V <- vcov(m1)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <-0
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.ST <- function(rho) sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
          rho <- optimize(logvero.ST, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas,rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        rho.old <- rho
        Delta.old <- Delta
        Gama.old <- Gama
        wi<-d.rho(z, rho, type=m.type, der=0)
        if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (loglik == TRUE){
        lk <- sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho,sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      } else {
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas,rho=rho, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      }
      class(obj.out) <- type
      return(obj.out)
  }

  if (type == "ST"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
      k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-d.rho(z, rho, type=m.type, der=0)
      teta <- c(betas,rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old <- rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error){
      count <- count + 1
      ####print(count)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

	    mu.1<-CGM(betas, x, type=model, der=0)
          mu<-mu.1+b*Delta*sqrt(wi)

          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / (Mtij)

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S22 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###

	    m1 <- nls(ymod ~ CGM(betas, x, type=model, der=0),start=list(betas=betas.old),trace=TRUE,weights=S1)
	    betas <- as.vector(coef(m1)) #sum(S1*y - Delta.old*S2) / sum(S1)
	    V <- vcov(m1)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <-sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.ST <- function(rho) sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
          rho <- optimize(logvero.ST, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas,rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        rho.old <- rho
        Delta.old <- Delta
        Gama.old <- Gama
        wi<-d.rho(z, rho, type=m.type, der=0)
        if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (loglik == TRUE){
        lk <- sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho,sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      } else {
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas,rho=rho, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      }
      class(obj.out) <- type
      return(obj.out)
  }

  if (type == "SCN"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-nu[1]/nu[2]^(1/2)+1-nu[1]
      k2<-nu[1]/nu[2]+1+nu[1]
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-d.rho(z, rho, type=m.type, der=0)
      teta <- c(betas,rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error){
      count <- count + 1
      ###print(count)

        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

	    mu.1<-CGM(betas, x, type=model, der=0)
          mu<-mu.1+b*Delta*sqrt(wi)

          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(wi*sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          u=(2/dSNC(y,mu,sigma2*wi,shape,nu))*(nu[1]*nu[2]*dnorm(y,mu,sqrt(wi*sigma2/nu[2]))*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2*wi))*pnorm(A,0,1))
          E=(2/dSNC(y,mu,sigma2*wi,shape,nu))*(nu[1]*sqrt(nu[2])*dnorm(y,mu,sqrt(sigma2*wi/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2*wi))*dnorm(A,0,1))

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S22 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###

	    m1 <- nls(ymod ~ CGM(betas, x, type=model, der=0),start=list(betas=betas.old),trace=TRUE,weights=S1)
          betas <-as.vector(coef(m1)) #sum(S1*y - Delta.old*S2) / sum(S1)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
    	    V <- vcov(m1)
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SCN <- function(rho) sum(log(dSNCH(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
          rho <- optimize(logvero.SCN, c(-5,-0.00001), tol = 0.0000001, maximum = TRUE)$maximum

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old <- rho
        wi<-d.rho(z, rho, type=m.type, der=0)
        if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (loglik == TRUE){
        lk <- sum(log(dSNCH(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas =betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      } else {                
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      }

      class(obj.out) <- type
      return(obj.out)
  }

  if (type == "SSL"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-2*nu/(2*nu-1)
      k2<-2*nu/(2*nu-2)
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      teta <- c(betas, rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      wi<-d.rho(z, rho, type=m.type, der=0)
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error){
      count <- count + 1
      ###print(count)
          u <- vector(mode="numeric",length=n)
          E <- vector(mode="numeric",length=n)
          S1 <- matrix(0, n)
          S2 <- matrix(0, n)
          S3 <- matrix(0, n)

          mu.1<-CGM(betas, x, type=model, der=0)
          mu<-mu.1+b*Delta*sqrt(wi)

          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(wi*sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b <-mutij+b
          A <- mutij / Mtij

          for(i in 1:n){
            E[i] <- (((2^(nu + 1))*nu*gamma(nu + 1))/(dSS(y[i], mu[i], sigma2*wi[i], shape, nu)*pi*sqrt(sigma2*wi[i])))* ((dj[i]+A[i]^2)^(-nu-1))*pgamma(1,nu+1,(dj[i]+A[i]^2)/2)
            faux <- function(u) u^(nu+0.5)*exp(-u*dj[i]/2)*pnorm(u^(1/2)*A[i])
            aux22 <- integrate(faux,0,1)$value
            u[i] <- ((sqrt(2)*nu) / (dSS(y[i], mu[i], sigma2*wi[i], shape, nu)*sqrt(pi)*sqrt(sigma2*wi[i])))*aux22
          }

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S22 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)


          ### M-step: atualizar mu, Delta, Gama, sigma2 ###

	    m1 <- nls(ymod ~ CGM(betas, x, type=model, der=0),start=list(betas=betas.old),trace=TRUE,weights=S1)
          betas <-as.vector(coef(m1)) #sum(S1*y - Delta.old*S2) / sum(S1)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
    	    V <- vcov(m1)
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SS <- function(rho) sum(log(dSSH(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
          rho <- optimize(logvero.SS, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old <- rho
        wi<-d.rho(z, rho, type=m.type, der=0)
        if (count==200) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (loglik == TRUE){
        lk <- sum(log(dSSH(rho, y, z, mu , sigma2 ,shape, nu, type=m.type)))
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu,loglik = lk, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      } else {
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu,iter = count, res=residuos, n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      }

      class(obj.out) <- type
      return(obj.out)
  }

  if (type == "SN"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-1
      k2=1
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-d.rho(z, rho, type=m.type, der=0)
      teta <- c(betas, rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error){
      count <- count + 1
      ###print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-CGM(betas, x, type=model, der=0)
          mu<-mu.1++b*Delta*sqrt(wi)

          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E = dnorm(mutij/Mtij) / pnorm(mutij/Mtij)
          u = rep(1, n)

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S22 <- (mutij_b*u + Mtij*E)
          S3 <- ((mutij_b)^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###

	    m1 <- nls(ymod ~ CGM(betas, x, type=model, der=0),start=list(betas=betas.old),trace=TRUE,weights=S1)
          betas <-as.vector(coef(m1)) #sum(S1*y - Delta.old*S2) / sum(S1)
    	    V <- vcov(m1)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SN <- function(rho) sum(log(dSNH(rho, y, z, mu , sigma2 ,shape, type=m.type)))
          rho <- optimize(logvero.SN, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old<-rho
        wi<-d.rho(z, rho, type=m.type, der=0)
        if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (loglik == TRUE){
        lk <- sum(log( dSNH (rho, y, z, mu , sigma2 ,shape, type=m.type) ))
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, loglik = lk, iter = count, res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      } else {
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, iter = count,res=residuos,n = length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      }

      class(obj.out) <- type
      return(obj.out)
  }



 if (type == "N"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-1
      k2=1
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-d.rho(z, rho, type=m.type, der=0)
      teta <- c(betas, rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error){
      count <- count + 1
      ###print(count)

        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-CGM(betas, x, type=model, der=0)
          mu<-mu.1+b*Delta*sqrt(wi)

          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E = dnorm(mutij/Mtij) / pnorm(mutij/Mtij)
          u = rep(1, n)

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S22 <- (mutij_b*u + Mtij*E)
          S3 <- ((mutij_b)^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###

	    m1 <- nls(ymod ~ CGM(betas, x, type=model, der=0),start=list(betas=betas.old),trace=TRUE,weights=S1)
          betas <-as.vector(coef(m1)) #sum(S1*y - Delta.old*S2) / sum(S1)
    	    V <- vcov(m1)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- 0
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SN <- function(rho) sum(log(dSNH(rho, y, z, mu , sigma2 ,shape, type=m.type)))
          rho <- optimize(logvero.SN, c(-5,-0.00001), tol = 0.0000001, maximum = TRUE)$maximum

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old<-rho
        wi<-d.rho(z, rho, type=m.type, der=0)
        if (count==200) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (loglik == TRUE){
        lk <- sum(log( dSNH (rho, y, z, mu , sigma2 ,shape, type=m.type) ))
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, loglik = lk, iter = count, res=residuos, n=length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      } else {
        obj.out <- list(y=y, x=x, u=u, uti=S22, uti2=S3, betas = betas, rho=rho, sigma2 = sigma2, shape = shape, iter = count,res=residuos, n=length(y), mah=dj, V=V, model=m1, curve = model, distr = type, m.type = m.type)
      }

      class(obj.out) <- type
      return(obj.out)
  }
}

