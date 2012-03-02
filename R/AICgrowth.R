AICgrowth <-
function(model)
{
L=model$loglik
type=model$distr
n=model$n
k=6
if(type=="N" | type=="T") k=5
AIC = -2 * (L - k)
AICc = AIC + 2 * k * (k + 1) / (n - k - 1)
BIC = -2 * L + k * log(n)
return(list(n=n, k=k, L=L, AIC=AIC, AICc=AICc, BIC=BIC))
}

