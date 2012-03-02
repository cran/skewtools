jn.var <-
function(data, betas, rho, sigma2, shape, nu, weight = FALSE, method="MLE", model = "VB", type = "ST", m.type = "power") 
{
n <- dim(data)[1]
k <- length(betas)
M <- matrix(NA, k, n)
estim <- V <- matrix(0, k, 1)
pos=c(1:n)
options(warn=-1)
	
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


	if(weight) {	
	for(i in 1:n){
		count=100*round(i/n,2)
		print(count)
		n.pos=pos[-i]
		if(method=="MLE") {
			st = c(b1=betas[1], b2=betas[2], b3=betas[3])
			if(model=="Logistic") ml <- nls(formula=y~b1/(1+exp(-b2*(x-b3))),data=data[n.pos,],start=st)
			if(model=="Gompertz") ml <- nls(formula=y~b1*(exp(-exp(-b2*(x-b3)))),data=data[n.pos,],start=st)
			if(model=="VB") ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3))),data=data[n.pos,],start=st)
			if(model=="Richards") {
				st = c(b1=betas[1], b2=betas[2], b3=betas[3], b4=betas[4])
				ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3)))^b4,data=data[n.pos,],start=st)
			}
			M[,i] <- summary(ml)$coefficients[,1]
		}
		if(method=="EM") {
			x=data$x[n.pos]; y=data$y[n.pos]
			parN <- HNL.skew(y, x, betas, rho, sigma2, shape, nu, loglik = TRUE, model, type, m.type, error = 0.00000001)
			M[,i] <- parN$betas
		}
		estim <- estim + M[,i]
	}
	estim <- estim/n
	for(i in 1:n) V <- V + (M[,i]-estim)^2
 	V <- sqrt(V*(n-1)/n)
	}
	
	else {
	eta <- CGM(betas, data$x, type=model, der=1)
	w_n <- det(t(eta)%*%eta)
	w <- rep(NA,n)

	for(i in 1:n){
		count=100*round(i/n,2)
		print(count)
		n.pos=pos[-i]
		if(method=="MLE") {
			st = c(b1=betas[1], b2=betas[2], b3=betas[3])
			if(model=="Logistic") ml <- nls(formula=y~b1/(1+exp(-b2*(x-b3))),data=data[n.pos,],start=st)
			if(model=="Gompertz") ml <- nls(formula=y~b1*(exp(-exp(-b2*(x-b3)))),data=data[n.pos,],start=st)
			if(model=="VB") ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3))),data=data[n.pos,],start=st)
			if(model=="Richards") {
				st = c(b1=betas[1], b2=betas[2], b3=betas[3], b4=betas[4])
				ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3)))^b4,data=data[n.pos,],start=st)
			}
			M[,i] <- summary(ml)$coefficients[,1]
		}
		if(method=="EM") {
			x=data$x[n.pos]; y=data$y[n.pos]
			parN <- HNL.skew(y, x, betas, rho, sigma2, shape, nu, loglik = TRUE, model, type, m.type, error = 0.00000001)
			M[,i] <- parN$betas
		}
		w[i] <- det(t(eta[n.pos,])%*%eta[n.pos,])
		estim <- estim + M[,i]
	}

	estim <- estim/n
	for(i in 1:n) V <- V + w[i]*(M[,i]-estim)^2
 	V <- sqrt(V/w_n)
	}
	
	return(list(V=V))
}

