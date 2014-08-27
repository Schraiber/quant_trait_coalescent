require(fBasics) #to simulate alpha stable guys, ?dstable
require(ape)
require(gsl) #for better bessel functions
require(stabledist)

jump_dif = function(sigma,lambda,kernel,...) {
	#sigma = rate*t of brownian motion
	#lambda = rate*t of poisson process
	#kernel = function to draw from jump kernel
	#... = parameters for kernel (MUST BE IN CORRECT ORDER!)
	numJumps = rpois(1,lambda)
	#Jbranch <<- append(Jbranch,numJumps)
	curState = 0
	curState = rnorm(1,0,sd=sigma)
	if (numJumps > 0) {
		curState = curState + sum(kernel(numJumps,...))
	}
	return(curState)	
}

sim_stable = function(phy, sigma, alpha, cee, sigma_tip=0) {
	#phy is ape-format tree
	#sigma is rate of brownian motion
	#cee is the scale parameter
	#alpha is the stability parameter
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma*sqrt(curLen)) + rstable(1,alpha=alpha,beta=0,gamma=(curLen)^(1/alpha)*cee)
	}
	return(rnorm(length(phy$tip.label),nodes[1:length(phy$tip.label)],sigma_tip))
	
}

sim_cauchy = function(phy, sigma, cee,sigma_tip = 0) {
	#phy is ape-format tree
	#sigma is rate of brownian motion
	#cee is the scale of the cauchy process
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma*sqrt(curLen)) + rcauchy(1, cee*curLen)
	}
	return(rnorm(length(phy$tip.label),nodes[1:length(phy$tip.label)],sigma_tip))
}

sim_vg = function(phy,sigma_bm,kap_vg,sigma_vg,sigma_tip=0) {
	#phy is ape-format tree
	#sigma_bm is brownian motion rate
	#kap_vg controls kurtosis of VG
	#sigma_vg controls rate of VG
	
	#implements as a time-changed brownian motion. BE SURE THAT THIS IS RIGHT!!!!!
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		curTimeChange = rgamma(1,curLen/kap_vg,scale=kap_vg)
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma_bm*sqrt(curLen)) + rnorm(1, sd = sigma_vg*sqrt(curTimeChange))
	}
	return(rnorm(length(phy$tip.label),nodes[1:length(phy$tip.label)],sigma_tip))
}

#this simulates a compound poisson model with jump kernel given by kernel and parameters given by ...
sim_cpp = function(phy, sigma, lambda, kernel,...,sigma_tip = 0) {
	#phy is an ape-format tree
	#sigma is rate of brownian motion
	#lambda is rate of poisson process
	#... are kernel parameters
	nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
	nodes[phy$edge[1,1]] = 0 #runif(1,-5,5) #sets the root
	for (i in 1:length(phy$edge[,1])) {
		curLen = phy$edge.length[i]
		nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + jump_dif(sigma*sqrt(curLen),lambda*curLen,kernel,...)
		
	}
	nodes <<- nodes
	return(rnorm(length(phy$tip.label),nodes[1:length(phy$tip.label)],sigma_tip))
}

LL_comp_gaus = function(parms,dat,t=2,cutoff=100) {
	#parms[1] = sig_bm, parms[2]=lam_jn, parms[3] = sig_jn
	LL = 0
	sds = sqrt((t*parms[1]^2 + 0:cutoff*parms[3]^2))
	for (i in 1:length(dat)) {
		LL = LL + log(sum(dpois(0:cutoff,parms[2]*t)*dnorm(dat[i],0,sds)))
	}
	return(-LL)
}

dJN = function(x,t,sig_bm,lambda,sig_jn,cutoff=100) {
	sds = sqrt((t*sig_bm^2+0:cutoff*sig_jn^2))
	res = vector();
	for (i in 1:length(x)) {	
		res[i] = sum(dpois(0:cutoff,lambda*t)*dnorm(x[i],0,sds))
	}
	return(res)
}

rVG = function(n,t,kap,sig) {
	gammas = rgamma(n,t/kap,scale=kap)
	to_return = rnorm(n,sd=sig*sqrt(gammas))
	to_return[to_return==0] = .Machine$double.xmin
	return(to_return)
}

dVG = function(x,t,kap,sig,log.p=F) {
	if (!log.p) {
	(2^(3/4 - t/(2*kap))*kap^(-(1/4) - t/(2*kap))*(sig^2/x^2)^((-2*t + kap)/(4*kap))*besselK(sqrt(2)/(sqrt(kap)*sqrt(sig^2/x^2)),abs(1/2 - t/kap)))/(sqrt(pi)*sqrt(sig^2)*gamma(t/kap))
	} else {
	(3/4-t/(2*kap))*log(2)+(-1/4-t/(2*kap))*log(kap) + (-2*t+kap)/(4*kap)*log(sig^2/x^2)+bessel_lnKnu(abs(1/2-t/kap),sqrt(2)/(sqrt(kap)*sqrt(sig^2/x^2)))-(1/2*log(pi*sig^2)+lgamma(t/kap))
	}
}

rVGBM = function(n,t,kap,sig_vg,sig_bm) {
	return(rVG(n,t,kap,sig_vg)+rnorm(n,sd=sig_bm*sqrt(t)))
}

#monte carlo estimation of the density by simulating a shitload of points
#x should be a vector
dVGBM_monte = function(x,t,kap,sig_vg,sig_bm,n=100000) {
	sims = rVGBM(n,t,kap,sig_vg,sig_bm)
	dens = density(sims)
	to_return = vector()
	for (i in 1:length(x)) {
		which_x = which.max(which(dens$x<x[i]))
		if (length(which_x)) {
			to_return[i] = dens$y[which_x]
		} else {
			to_return[i] = .Machine$double.xmin
		}
		
	}
	return(to_return)
}

LL_VG_3_taxon = function(parms,dat,jumps,brlen) {
	#parms[1] = sig_bm, parms[2] = kap_vg, parms[3] = sig_vg
	#dat[1] and dat[2] are sisters, dat[3] is outgroup
	#brlen[1] is 1->(1,2), brlen[2] is 2->(1,2), brlen[3] is 3->(1,2,3) and brlen[4] is (1,2)->(1,2,3)
	#jumps[1] is from (1,2)->1, jumps[2] is from (1,2)->2, jumps[3] is from (1,2,3)->3, and jumps[4] is from (1,2,3)->(1,2)
	
	#compute the jump likelihood
	kj = sum(log(dVG(jumps,brlen,parms[2],parms[3])))
	
	#compute the brownian motion part
	kb = 0
	#node (1,2)
	mu1 = dat[1]-jumps[1]
	mu2 = dat[2]-jumps[2]
	t1 = brlen[1]
	t2 = brlen[2]
	mu12 = (mu1*t2+mu2*t1)/(t1+t2)-jumps[4]
	sig12 = sqrt(parms[1]^2*t1*t2/(t1+t2))
	kb = kb + -(mu1-mu2)^2/(2*parms[1]^2*(t1+t2))-log(sqrt(2*pi*parms[1]^2*(t1+t2))) 
	#node (1,2,3)
	mu3 = dat[3]-jumps[3]
	t3 = brlen[3]
	t12 = brlen[4] + sig12^2/parms[1]^2
	mu123 = (mu12*t3+mu3*t12)/(t12+t3)
	sig123 = sqrt(parms[1]^2*t12*t3/(t12+t3))
	kb = kb + -(mu12-mu3)^2/(2*parms[1]^2*(t12+t3)) - log(sqrt(2*pi*parms[1]^2*(t12+t3)))
	
	#return the likelihood
	return(list(kb=kb,kj=kj,lnL=kb + kj))

		
}

#includes the jumps as part of the likelihood---full data likelihood
LL_VG_timeseries_jumps = function(parms,dat,jumps,t) {
	#parms[1] = sig_bm, parms[2] = kap_vg, parms[3] = sig_vg
	
	#compute jump likelihood
	kj = sum(log(dVG(jumps,t,parms[2],parms[3])))
	
	#compute brownian likelihood
	kb = 0
	for (i in 1:(length(dat)-1)) {
		kb = kb + dnorm(dat[i+1]-jumps[i]-dat[i],mean=0,parms[1]*sqrt(t),log=T)
	}
	c(kj=kj,kb=kb,lnL = kb+kj)
}

#uses an approximation of the integrated likelihood
LL_VG_timeseries = function(parms,dat,t,n=100000) {
	contrast = dat[2:length(dat)]-dat[1:(length(dat)-1)]
	return(sum(log(dVGBM_monte(contrast,t,parms[2],parms[3],parms[1],n=n))))
}

#sample paths of a BM + jump process, broken down in the BM and drift components
#parameters should come in time scaled already, if necessary
sample_path = function(n,len,sigma.bm,jump.process,...) {
	BM = matrix(0,nrow=n,ncol=len)
	jump = matrix(0,nrow=n,ncol=len)
	for (i in 1:n) {
		BM[i,] = cumsum(rnorm(len,sd=sigma.bm))
		jump[i,] = cumsum(jump.process(len,...))	
	}
	combined = BM+jump
	return(list(BM=BM,jump=jump,combined=combined))
}

#samples compound poisson path---sigma.bm should be scaled, the others NOT
compound_poisson = function(n,T,len,sigma.bm,lambda,sigma.jn) {
	BM = matrix(0,nrow=n,ncol=len)
	jump = matrix(0,nrow=n,ncol=len)
	t = seq(0,T,length=len)
	for (i in 1:n) {
		BM[i,] = cumsum(rnorm(len,sd=sigma.bm))
		num.jumps = rpois(1,lambda*T)
		if (num.jumps > 0){ 
			jump.times = runif(num.jumps,0,T)
			jump.sizes = rnorm(num.jumps,0,sd=sigma.jn)
			jump.bins = sapply(jump.times, function(x){max(which(t<x))})
			jump[i,jump.bins] = jump.sizes
			jump[i,] = cumsum(jump[i,])
		}
	}
	combined = BM+jump
	return(list(BM=BM,jump=jump,combined=combined))
}

as_levy = function(x,alpha,c) {
	c^alpha/abs(x)^(1+alpha)
}

vg_levy = function(x,kappa,sigma) {
	1/(kappa*abs(x))*exp(-sqrt(2*x^2/(kappa*sigma^2)))
}

#this simulates the data and generates the figure
#pars should be a list where the elements are jn, vg and as
make_sample_path_figure = function(pars,T=1,len=1000,per.model=3,replicates=2) {
	layout(matrix(c(1,1,2,3,3,4,5,5,6),byrow=T,nrow=3,ncol=3))
	sims = list()
	t = seq(0,T,len=len)
	dt = t[2]-t[1]
	x = seq(-5,5,.001)
	#compound poisson
	#simulate the data
	for (i in 1:per.model) {
		sims[[i]] = compound_poisson(replicates,T,len,pars$jn[1+3*(i-1)],pars$jn[2+3*(i-1)],pars$jn[3+3*(i-1)])
	}
	#find the axis limits
	min.max = unlist(lapply(sims,function(x){apply(x$jump,1,function(y){c(min(y),max(y))})}))
	min.sampled = min(min.max)
	max.sampled = max(min.max)
	#open the plot
	plot("",xlim=c(0,T),ylim=c(min.sampled,max.sampled),xaxt="n",ylab="",xlab="")
	#loop over the samples
	for (i in 1:per.model) {
		matlines(t,t(sims[[i]]$jump),lty=i,col="black")	
	}
	axis(1,labels=F)
	#draw the Levy measures
	levy = matrix(nrow=per.model,ncol=length(x))
	for (i in 1:per.model) {
		levy[i,] = pars$jn[2+3*(i-1)]*dnorm(x,sd=pars$jn[3+3*(i-1)])
	}
	#plot the Levy measures
	plot("",xlim=c(-5,5),ylim=c(0,max(levy)),yaxt="n",xaxt="n",ylab="",xlab="")
	for (i in 1:per.model) {
		lines(x,levy[i,],lty=i)	
	}
	axis(1,labels=F)
	#variance gamma
	for (i in 1:per.model) {
		#fix this shit now
		sims[[i]] = sample_path(replicates,len,pars$vg[1+3*(i-1)],rVG,dt,pars$vg[2+3*(i-1)],pars$vg[3+3*(i-1)])
	}
	#find the axis limits
	min.max = unlist(lapply(sims,function(x){apply(x$jump,1,function(y){c(min(y),max(y))})}))
	min.sampled = min(min.max)
	max.sampled = max(min.max)
	#open the plot
	plot("",xlim=c(0,T),ylim=c(min.sampled,max.sampled),xaxt="n",ylab="",xlab="")
	#loop over the samples
	for (i in 1:per.model) {
		matlines(t,t(sims[[i]]$jump),lty=i,col="black")	
	}
	axis(1,labels=F)
#draw the Levy measures
	levy = matrix(nrow=per.model,ncol=length(x))
	for (i in 1:per.model) {
		levy[i,] = vg_levy(x,pars$vg[2+3*(i-1)],pars$vg[3+3*(i-1)])
	}
	#plot the Levy measures
	plot("",xlim=c(-5,5),ylim=c(0,1),yaxt="n",xaxt="n",ylab="",xlab="")
	for (i in 1:per.model) {
		lines(x,levy[i,],lty=i)	
	}
	axis(1,labels=F)
	#alpha stable
	for (i in 1:per.model) {
		#fix this shit now
		sims[[i]] = sample_path(replicates,len,pars$vg[1+3*(i-1)],rstable,alpha=pars$as[2+3*(i-1)],beta=0,gamma=pars$as[3+3*(i-1)],pm=2)
	}
	#find the axis limits
	min.max = unlist(lapply(sims,function(x){apply(x$jump,1,function(y){c(min(y),max(y))})}))
	min.sampled = min(min.max)
	max.sampled = max(min.max)
	#open the plot
	plot("",xlim=c(0,T),ylim=c(min.sampled,max.sampled),ylab="",xlab="")
	#loop over the samples
	for (i in 1:per.model) {
		matlines(t,t(sims[[i]]$jump),lty=i,col="black")	
	}
#draw the Levy measures
	levy = matrix(nrow=per.model,ncol=length(x))
	for (i in 1:per.model) {
		levy[i,] = as_levy(x,pars$as[2+3*(i-1)],pars$as[3+3*(i-1)])
	}
	#plot the Levy measures
	plot("",xlim=c(-5,5),ylim=c(0,1),yaxt="n",ylab="",xlab="")
	for (i in 1:per.model) {
		lines(x,levy[i,],lty=i)	
	}

}

