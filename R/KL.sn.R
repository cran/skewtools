
KL.sn <-
function(beta1,beta2){

C1=is.numeric(dim(beta1)); C2=is.numeric(dim(beta2))
if(C1 & C2) stop("This function is only available for univariate case")

M1=beta1[1]; O1=beta1[2]; E1=beta1[3]; delta1=O1*E1/sqrt(1+O1*E1^2)
M2=beta2[1]; O2=beta2[2]; E2=beta2[3]; delta2=O2*E2/sqrt(1+O2*E2^2)
V1=O1-(2/pi)*delta1^2; V2=O2-(2/pi)*delta2^2
ME=0.5*(1+log(2*pi)+c(log(V1),log(V2)))   

Emint <- function(beta1, beta2){
delta=beta1[2]*beta1[3]/sqrt(1+beta1[2]*beta1[3]^2)
muY=beta2[3]*(beta1[1]-beta2[1])
sigmaY=beta1[2]*beta2[3]^2
lambdaY=beta2[3]*delta/sqrt(beta1[2]*beta2[3]^2-(beta2[3]*delta)^2)
f1 <- function(x) {
2*dnorm(x, muY, sqrt(sigmaY))*pnorm(lambdaY*x, muY, sqrt(sigmaY))*log(2*pnorm(x))
} 
W <- seq(-50,0,1)
a=try(integrate(f1,W[1],Inf, subdivisions=100)$value, TRUE)
k=2
while(!is.numeric(a)) {
a=try(integrate(f1,W[k],0, subdivisions=100)$value + integrate(f1,0,Inf, subdivisions=100)$value, TRUE) 
k=k+1
}
return(a)
}

DKL0 = 0.5*(log(O2/O1) + O1/O2 + (M1-M2)^2/O2 - 1)
KL = DKL0 + sqrt(2/pi)*(M1-M2)*delta1/O2 + Emint(beta1, beta1) - Emint(beta1, beta2)
entropy = 0.5*(1+log(2*pi)) + 0.5*c(log(O1), log(O2)) - c(Emint(beta1, beta1), Emint(beta2, beta2))

return(list(Dkl=as.numeric(KL), max.entropy=ME))
}

