conf.im <-
function(x, y, model, alpha=0.05, plot.it=TRUE)
{

betas=model$betas
rho=model$rho
sigma2=model$sigma2
lambda=model$shape
curva=model$curve
type=model$distr

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

	if(type=="NIST") {
		if(der==0) d <- exp(-betas[1]*x)/(betas[2]+betas[3]*x)
		if(der==1) {
		db1 <- -(exp(-betas[1] * x) * x/(betas[2] + betas[3] * x))
		db2 <- -(exp(-betas[1] * x)/(betas[2] + betas[3] * x)^2)
		db3 <- -(exp(-betas[1] * x) * x/(betas[2] + betas[3] * x)^2)
		}
		if(der==2) {
		db1b1 <- exp(-betas[1] * x) * x * x/(betas[2] + betas[3] * x)
		db2b2 <- exp(-betas[1] * x) * (2 * (betas[2] + betas[3] * x))/((betas[2] + betas[3] * x)^2)^2
		db3b3 <- exp(-betas[1] * x) * x * (2 * (x * (betas[2] + betas[3] * x)))/((betas[2] + betas[3] * x)^2)^2
		db1b2 <- exp(-betas[1] * x) * x/(betas[2] + betas[3] * x)^2
		db1b3 <- exp(-betas[1] * x) * x * x/(betas[2] + betas[3] * x)^2
		db2b3 <- exp(-betas[1] * x) * (2 * (x * (betas[2] + betas[3] * x)))/((betas[2] + betas[3] * x)^2)^2
		}
	}

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

	if(type=="GVB") {
		if(der==0) d <- betas[1]*(1-exp(-betas[2]*(1-betas[4])*(x-betas[3]))^(1/(1-betas[4])))
		if(der==1) {
		db1 <- (1 - exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(1/(1 - betas[4])))
		db2 <- betas[1] * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * ((1 - betas[4]) * (x -     betas[3])))))
		db3 <- betas[1] * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (-betas[2] * (1 - betas[4])))))
		db4 <- -(betas[1] * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (betas[2] * (x - betas[3])))) +     exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(1/(1 - betas[4])) * (log(exp(-betas[2] *         (1 - betas[4]) * (x - betas[3]))) * (1/(1 - betas[4])^2))))
		}
		if(der==2) {
		db1b1 <- 0
		db2b2 <- -(betas[1] * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * ((1 - betas[4]) * (x -     betas[3])) * ((1 - betas[4]) * (x - betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x -     betas[3]))^(((1/(1 - betas[4])) - 1) - 1) * (((1/(1 - betas[4])) - 1) * (exp(-betas[2] *     (1 - betas[4]) * (x - betas[3])) * ((1 - betas[4]) * (x - betas[3])))) * ((1/(1 -     betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * ((1 - betas[4]) * (x -     betas[3]))))))
		db3b3 <- -(betas[1] * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (-betas[2] * (1 - betas[4])) *     (-betas[2] * (1 - betas[4])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(((1/(1 -     betas[4])) - 1) - 1) * (((1/(1 - betas[4])) - 1) * (exp(-betas[2] * (1 - betas[4]) *     (x - betas[3])) * (-betas[2] * (1 - betas[4])))) * ((1/(1 - betas[4])) * (exp(-betas[2] *     (1 - betas[4]) * (x - betas[3])) * (-betas[2] * (1 - betas[4]))))))
		db4b4 <- -(betas[1] * ((exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(((1/(1 - betas[4])) - 1) - 1) * (((1/(1 - betas[4])) - 1) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     (betas[2] * (x - betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 -     betas[4])) - 1) * (log(exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))) * (1/(1 -     betas[4])^2))) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     (betas[2] * (x - betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 -     betas[4])) - 1) * (1/(1 - betas[4])^2 * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     (betas[2] * (x - betas[3]))) + (1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x -     betas[3])) * (betas[2] * (x - betas[3])) * (betas[2] * (x - betas[3])))) + ((exp(-betas[2] *     (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) *     (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (betas[2] * (x - betas[3])))) + exp(-betas[2] *     (1 - betas[4]) * (x - betas[3]))^(1/(1 - betas[4])) * (log(exp(-betas[2] * (1 - betas[4]) *     (x - betas[3]))) * (1/(1 - betas[4])^2))) * (log(exp(-betas[2] * (1 - betas[4]) *     (x - betas[3]))) * (1/(1 - betas[4])^2)) + exp(-betas[2] * (1 - betas[4]) * (x -     betas[3]))^(1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (betas[2] *     (x - betas[3]))/exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (1/(1 - betas[4])^2) +     log(exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))) * (2 * (1 - betas[4])/((1 -         betas[4])^2)^2)))))	      
     		db1b2 <- exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * ((1 - betas[4]) * (x -     betas[3]))))
		db1b3 <- exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (-betas[2] * (1 - betas[4]))))
     		db1b4 <- -(exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (betas[2] * (x - betas[3])))) +     exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(1/(1 - betas[4])) * (log(exp(-betas[2] *         (1 - betas[4]) * (x - betas[3]))) * (1/(1 - betas[4])^2)))
		db2b3 <- -(betas[1] * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 - betas[4])) - 1) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (1 - betas[4]) + exp(-betas[2] *     (1 - betas[4]) * (x - betas[3])) * (-betas[2] * (1 - betas[4])) * ((1 - betas[4]) * (x -     betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(((1/(1 - betas[4])) -     1) - 1) * (((1/(1 - betas[4])) - 1) * (exp(-betas[2] * (1 - betas[4]) * (x -     betas[3])) * (-betas[2] * (1 - betas[4])))) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 -     betas[4]) * (x - betas[3])) * ((1 - betas[4]) * (x - betas[3]))))))
     		db2b4 <- betas[1] * ((exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(((1/(1 - betas[4])) - 1) - 1) * (((1/(1 - betas[4])) - 1) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (betas[2] *         (x - betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 -     betas[4])) - 1) * (log(exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))) * (1/(1 -     betas[4])^2))) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     ((1 - betas[4]) * (x - betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 -     betas[4])) - 1) * (1/(1 - betas[4])^2 * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     ((1 - betas[4]) * (x - betas[3]))) + (1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) *     (x - betas[3])) * (betas[2] * (x - betas[3])) * ((1 - betas[4]) * (x - betas[3])) - exp(-betas[2] *     (1 - betas[4]) * (x - betas[3])) * (x - betas[3]))))
     		db3b4 <- betas[1] * ((exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^(((1/(1 - betas[4])) - 1) - 1) * (((1/(1 - betas[4])) - 1) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) * (betas[2] *         (x - betas[3])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 -     betas[4])) - 1) * (log(exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))) * (1/(1 -     betas[4])^2))) * ((1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     (-betas[2] * (1 - betas[4])))) + exp(-betas[2] * (1 - betas[4]) * (x - betas[3]))^((1/(1 -     betas[4])) - 1) * (1/(1 - betas[4])^2 * (exp(-betas[2] * (1 - betas[4]) * (x - betas[3])) *     (-betas[2] * (1 - betas[4]))) + (1/(1 - betas[4])) * (exp(-betas[2] * (1 - betas[4]) *     (x - betas[3])) * (betas[2] * (x - betas[3])) * (-betas[2] * (1 - betas[4])) + exp(-betas[2] *     (1 - betas[4]) * (x - betas[3])) * betas[2])))		
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
		if(type=="Logistic" | type=="VB" | type=="Gompertz" | type=="NIST") deriv1 <- cbind(db1,db2,db3)
		if(type=="Richards" | type=="GVB") deriv1 <- cbind(db1,db2,db3,db4)
		L=deriv1
	} 
	if(der==2) {
		if(type=="Logistic" | type=="VB" | type=="Gompertz" | type=="NIST") {
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

if (type == "T" | type == "ST") {
	nu=model$nu
	k1=sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
	k2=(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
}
if (type == "SN" | type == "N") k1=k2=1
if (type == "SSL") {
	nu=model$nu	
	k1=2*nu/(2*nu-1)
	k2=2*nu/(2*nu-2)
}

SE <- function(x){
V=(k2-(-sqrt(2/pi)*k1*lambda/(sqrt(1+lambda^2)))^2)*sigma2*d.rho(x, rho, type=model$m.type, der=0)
1.96*sqrt(V)
}

yeval=CGM(betas, x, curva, der=0)
ICS=yeval+SE(x)
ICI=yeval-SE(x)

if(plot.it) {
plot(x, y, col="gray", ylab="Longitude", xlab="Age", main=paste("",type,"- EM fit"),
xlim=c(0, max(x)), ylim=c(-5, max(ICS,y)))
x.points=seq(0,max(x),0.5)
yeval1=CGM(betas, x.points, curva, der=0)
SE1=SE(x.points)
points(x.points, CGM(betas, x.points, curva, der=0), lwd=2, type="l")
points(x.points, yeval1 + SE1, col="red", lty=2, type="l")
points(x.points, yeval1 - SE1, col="red", lty=2, type="l")
}

return(list(yeval=yeval, ICS=ICS, ICI=ICI))
}

