search.nu <-
function(y, x, betas, rho, sigma2, shape, nu=c(2:15), model="VB", type="ST", m.type="power", plot.it=TRUE)
{

loglik <- nu
B <- matrix(NA, length(nu), length(betas)+4)
rownames(B)=nu
if(model=="Richards") colnames(B)=c("b1", "b2", "b3", "b4", "rho", "sigma2", "lambda", "loglik")
if(model!="Richards") colnames(B)=c("b1", "b2", "b3", "rho", "sigma2", "lambda", "loglik")

for(i in 1:length(nu)){ 
	par <- HNL.skew(y, x, betas, rho, sigma2, shape, nu[i], loglik=TRUE, model, type, m.type, error=0.00001)
	loglik[i] = par$loglik
	B[i,] = c(par$betas, par$rho, par$sigma2, par$shape, loglik[i])
}

l.def=loglik[2:length(nu)]
nu.max=nu[loglik==max(loglik, na.rm = TRUE)]

if(plot.it) {
plot(nu, loglik, type="l", main="Profile log-likelihood")
abline(v=nu.max, lty=3, col="blue")
}

return(list(table.res=B, pos.nu=nu.max))
}

