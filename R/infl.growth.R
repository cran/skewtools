infl.growth <-
function(y, x, model, plot.it=TRUE)
{

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
	
H.SNI <-
function(y, x, betas, rho, sigma2, lambda, nu, model, type,  m.type)
{

n=length(y)
npar=length(betas)
yeval = CGM(betas, x, type=model, der=0)
D1 = CGM(betas, x, type=model, der=1)
m=d.rho(x, rho, type=m.type, der=0)
Dm=d.rho(x, rho, type=m.type, der=1)
D2m=d.rho(x, rho, type=m.type, der=2)
sig2=sigma2*m

if (type == "N"){
	lambda=delta=delta1=delta2=0
	b=-sqrt(2/pi)
	C = (y-yeval)/sqrt(sig2)
	B = C - b*delta
	A = lambda*B/sqrt(sig2)
	d = B^2
	nu.large=150
	IPhi <- function(w) (2^w)*(nu.large^(nu/2))*gamma((2*w+nu.large)/2)*pt(A*sqrt((2*w+nu.large)/(d+nu.large)),2*w+nu.large)/(gamma(nu.large/2)*(nu.large+d)^((2*w+nu.large)/2))
	Iphi <- function(w) (2^w)*(nu.large^(nu/2))*gamma((2*w+nu.large)/2)/(sqrt(2*pi)*gamma(nu.large/2)*(nu.large+d+A^2)^((2*w+nu.large)/2))
}

if (type == "T"){
	lambda=delta=delta1=delta2=0
	b=-sqrt(2/pi)*sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
	C = (y-yeval)/sqrt(sig2)
	B = C - b*delta
	A = lambda*B/sqrt(sig2)
	d = B^2
	IPhi <- function(w) (2^w)*(nu^(nu/2))*gamma((2*w+nu)/2)*pt(A*sqrt((2*w+nu)/(d+nu)),2*w+nu)/(gamma(nu/2)*(nu+d)^((2*w+nu)/2))
	Iphi <- function(w) (2^w)*(nu^(nu/2))*gamma((2*w+nu)/2)/(sqrt(2*pi)*gamma(nu/2)*(nu+d+A^2)^((2*w+nu)/2))
}

if (type == "SN"){
	delta=lambda/sqrt(1+lambda^2)
	delta1=(1+lambda^2)^(-3/2)
	delta2=-3*lambda*(1+lambda^2)^(-5/2)
	b=-sqrt(2/pi)
	C = (y-yeval)/sqrt(sig2)
	B = C - b*delta
	A = lambda*B/sqrt(sig2)
	d = B^2
	nu.large=150
	IPhi <- function(w) (2^w)*(nu.large^(nu/2))*gamma((2*w+nu.large)/2)*pt(A*sqrt((2*w+nu.large)/(d+nu.large)),2*w+nu.large)/(gamma(nu.large/2)*(nu.large+d)^((2*w+nu.large)/2))
	Iphi <- function(w) (2^w)*(nu.large^(nu/2))*gamma((2*w+nu.large)/2)/(sqrt(2*pi)*gamma(nu.large/2)*(nu.large+d+A^2)^((2*w+nu.large)/2))
}

if (type == "ST"){
	delta=lambda/sqrt(1+lambda^2)
	delta1=(1+lambda^2)^(-3/2)
	delta2=-3*lambda*(1+lambda^2)^(-5/2)
	b = -sqrt(2/pi)*sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2) 
	C = (y-yeval)/sqrt(sig2)
	B = C - b*delta
	A = lambda*B/sqrt(sig2)
	d = B^2
	IPhi <- function(w) (2^w)*(nu^(nu/2))*gamma((2*w+nu)/2)*pt(A*sqrt((2*w+nu)/(d+nu)),2*w+nu)/(gamma(nu/2)*(nu+d)^((2*w+nu)/2))
	Iphi <- function(w) (2^w)*(nu^(nu/2))*gamma((2*w+nu)/2)/(sqrt(2*pi)*gamma(nu/2)*(nu+d+A^2)^((2*w+nu)/2))
}


	
DA <- Dd <- matrix(NA,npar,n)

for(i in 1:npar) {
	DA[i,]=-lambda*D1[,i]/sig2
	Dd[i,]=-2*B*D1[,i]/sqrt(sig2)
}

Iaux <- matrix(0,n,npar)
Ui <- matrix(NA,n,npar)
I1=Iphi(1)
I2=IPhi(3/2)
I6=IPhi(1/2)
for(i in 1:npar) Ui[,i] = Iaux[,i] + (I1*DA[i,] - 0.5*I2*Dd[i,])/I6
return(Ui)
}

	betas = model$betas
	rho = model$rho
	sigma2 = model$sigma2
	n = model$n
	p = length(betas)
	wi = d.rho(x, rho, type=model$m.type, der=0)
	mu.1 = CGM(betas, x, type=model$curve, der=0)
	type=model$distr
	
	if (type == "T" | type == "ST"){
		shape = model$shape
		if(type == "ST") nu = model$nu
		k1 = sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
		k2 = (nu/2)*gamma((nu-2)/2)/gamma(nu/2)
		tau=6
		V <- H.SNI(y, x, betas, rho, sigma2, lambda=shape, nu, model=model$curve, type=type, m.type=model$m.type)
	}	
	if (type == "N" | type == "SN"){
		shape = model$shape
		k1<-1
		k2<-1
		tau=3
		V <- H.SNI(y, x, betas, rho, sigma2, lambda=shape, nu=0, model=model$curve, type=type, m.type=model$m.type)	
	}
	
	lambda=shape
	vari = sigma2*wi*(k2-2/pi*k1^2*lambda^2/(1+lambda^2))
	residuos <-(y-mu.1)/sqrt(vari)
	
	Delta <- V
	mioST <- model$V
   	IfST = abs((Delta)%*%mioST%*%t(Delta))
   	MoST <- diag(IfST)/sum(diag(IfST))

	if(plot.it) {
		#par(mfrow=c(1,2))
		#qqnorm(residuos, main="Normal Q-Q Plot", xlab="Standard Normal Quantiles", ylab="Pearson Residuals", ylim=c(-6,11), type='p', lty=10) 
		#qqline(residuos)
   		plot(MoST,  type='h' , xlim=c(0,n), ylim=c(0,max(MoST)), ylab="M(O)", xlab="Index", main=model$distr)
		level=mean(MoST)+tau*sd(MoST)
		abline(h=level,lty=2)
		abline(h=0,lty=1)
		pos <- rep(0,n)
		for(i in 1:n) {
			if(MoST[i] >= level) {
				pos[i]=i
				text(i, MoST[i], i, cex=0.7)
			}
		}
	}

	outliers=pos[pos!=0]
	return(list(residuals=residuos, outliers=outliers))
}

