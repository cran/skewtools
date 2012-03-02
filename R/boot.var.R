boot.var <-
function(data, betas, rho, sigma2, shape, nu, B, size, method = "EM", type = "N", model = "VB", m.type = "power", replace = TRUE, category = TRUE, prop = TRUE) 
{
M <- matrix(NA, length(betas), B)
estim <- V <- matrix(0, length(betas), 1)
n.pos=c(1:length(data$x))
options(warn=-1)

age.category <-
function(x, size, replace=FALSE, prop=TRUE)
{
	# x: vector of ages
	# size: bootstrap sample size	

	cat=as.numeric(names(table(x)))
	k=length(cat)
	Kn=as.numeric(table(x))
	n=length(x)
	pos=1:n	
	
	if(prop==FALSE){
		A=NULL
		if(size <= min(Kn)) B=size
		if(size > min(Kn) & replace==TRUE) B=round(size/k,0)
		if(size > min(Kn) & replace==FALSE) break
		for(i in 1:k) {
			data.new=pos[x==cat[i]]
			pos.new=sample(data.new, B, replace)
			A=c(A,pos.new)
		}
	}
	
	if(prop==TRUE){
		weights=Kn/n
		data.new=pos[x==cat[1]]
		B=round(weights[1]*size,0)
		pos.new=sample(data.new, B, replace)
		A=pos.new

		for(i in 2:k) {
			data.new=pos[x==cat[i]]
			B=round(weights[i]*size,0)
			pos.new=sample(data.new, B, replace)
			A=c(A,pos.new)
		}
	}	
	
	A=A[is.na(A)==FALSE]
	return(A)
}

	for(i in 1:B){
		count=100*i/B
		print(count)
		if(category==FALSE) pos=sample(n.pos, size, replace)
		if(category==TRUE) pos=age.category(data$x, size, replace, prop)
		if(method=="MLE") {
			st = c(b1=betas[1], b2=betas[2], b3=betas[3])
			if(model=="NIST") 	ml <- nls(formula=y~exp(-b1*x)/(b2 + b3*x),data=data[pos,],start=st)
			if(model=="Logistic") 	ml <- nls(formula=y~b1/(1+exp(-b2*(x-b3))),data=data[pos,],start=st)
			if(model=="Gompertz") 	ml <- nls(formula=y~b1*(exp(-exp(-b2*(x-b3)))),data=data[pos,],start=st)
			if(model=="VB") 		ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3))),data=data[pos,],start=st)
			if(model=="Richards") {
				st = c(b1=betas[1], b2=betas[2], b3=betas[3], b4=betas[4])
				ml <- nls(formula=y~b1*(1-exp(-b2*(x-b3)))^b4,data=data[pos,],start=st)
			}
			M[,i] <- summary(ml)$coefficients[,1]
		}
		if(method=="EM") {
			x=data$x
			y=data$y
			parN <- HNL.skew(y=y[pos], x=x[pos], betas, rho, sigma2, shape, nu, loglik = TRUE, model, type, m.type, error = 0.00000001)
			M[,i] <- parN$betas
		}
		estim <- estim + M[,i]
	}
	
	estim <- estim/B
	for(i in 1:B) V <- V + (M[,i]-estim)^2
 	V <- sqrt(V/B)
	return(list(V=V))
}

