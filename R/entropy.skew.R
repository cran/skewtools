entropy.skew <-
function(X, family="ST"){

norm <- function(wYY,wZZ,wYZ,wZY,eY,eZ){
wZZy=wZZ-wZY%*%solve(wYY)%*%wYZ 		
wYYz=wYY-wYZ%*%solve(wZZ)%*%wZY
eYz=as.vector(eY+solve(wYY)%*%wYZ%*%eZ)/sqrt(1+t(eZ)%*%wZZy%*%eZ)
eZy=as.vector(eZ+solve(wZZ)%*%wZY%*%eY)/sqrt(1+t(eY)%*%wYYz%*%eY)
NeY=sqrt(t(eYz)%*%wYY%*%eYz)
NeZ=sqrt(t(eZy)%*%wZZ%*%eZy)
NeYZ=sqrt(NeY^2 + NeZ^2 + 2*t(eYz)%*%wYZ%*%eZy)
list(NeY=NeY, NeZ=NeZ, NeYZ=NeYZ)	
}

Emint <- function(a){
f1 <- function(x) 2*pnorm(a*x)*log(2*pnorm(a*x))*dnorm(x) 
x=37.5/a
integrate(f1,-x,100, subdivisions=100)$value
}

Emint2 <- function(a,nu,k){
h01= log(gamma(0.5*(nu+k)))-log(gamma(0.5*nu))-0.5*k*log(nu*pi)
h02= 0.5*(nu+k)*(digamma(0.5*(nu+k))-digamma(0.5*nu)) 
f1 <- function(x) log(2*pt(sqrt(k+nu)*a*x/sqrt(nu+k-1+x^2),nu+k))*dst(x,0,1,a,nu+k-1)
-h01 + h02 - integrate(f1,-100,100, subdivisions=100)$value
}

if(family=="ST") {
model=mst.fit(y=X,plot.it=FALSE)$dp
O=model$Omega; E=model$alpha; n.eta=sqrt(t(E)%*%O%*%E); k=1
if(!is.null(dim(X))) k=dim(X)[2]
H = 0.5*log(det(O)) + Emint2(n.eta,nu=model$df,k=k)
}
	
if(family=="SN") {
model=msn.fit(y=X, plot.it=FALSE)$dp
O=model$Omega; E=model$alpha; n.eta=sqrt(t(E)%*%O%*%E); k=1
if(!is.null(dim(X))) k=dim(X)[2]
H = 0.5*log(det(O)) + 0.5*k*(1+log(2*pi)) - Emint(n.eta)
}

return(list(H=H, family=family))
}

