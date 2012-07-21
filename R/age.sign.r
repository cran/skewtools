
age.sign <- function(y, model) {

distr = model$distr
funcA = model$curve
n=length(y)

geraSkewNormal <- function(mu, sigma2, shape){
n=length(sigma2)
delta <- shape / sqrt(1 + shape^2)
y <- mu+ sqrt(sigma2)*(delta*abs(rnorm(n)) + (1 - delta^2)^(1/2)*rnorm(n))
return(y)
}

geraSNIt <- function(mu, sigma2, shape, nu){
n=length(sigma2)
y <- mu + (rgamma(n, nu/2, nu/2))^(-1/2)*geraSkewNormal(0, sigma2, shape)
}

Age <- function(y) {
Linf=model$betas[1]; k=model$betas[2]; t0=model$betas[3]
edad=t0-log(1-y/Linf)/k
return(edad)
}

edad = error = rep(NA, n)

for(i in 1:n){
if(distr=="N" | distr=="SN") error[i] = geraSkewNormal(model$rho, model$sigma2, model$shape)
if(distr=="T" | distr=="ST") error[i] = geraSNIt(model$rho, model$sigma2, model$shape, model$nu)
edad[i] = Age(y[i]) + error[i]
}

return(list(edad=edad, error=error))
}