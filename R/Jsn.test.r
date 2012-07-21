
Jsn.test <-
function(X,Y){

C1=is.numeric(dim(X)); C2=is.numeric(dim(Y))
if(C1 & C2) stop("This function is only available for univariate case")
N1=length(X); N2=length(Y)

modelX=msn.fit(y=X, plot.it=FALSE)$dp; beta1=c(modelX$beta, modelX$Omega, modelX$alpha)
modelY=msn.fit(y=Y, plot.it=FALSE)$dp; beta2=c(modelY$beta, modelY$Omega, modelY$alpha)
J = KL.sn(beta1,beta2)$Dkl + KL.sn(beta2,beta1)$Dkl

STATISTIC = N1*N2*J/(N1+N2)
PVAL = pchisq(STATISTIC, df=3, lower.tail=FALSE)
structure(list(STATISTIC=STATISTIC, p.value=PVAL))
}





