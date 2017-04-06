#SDS 383D - Exercise 4
#Function
#Jennifer Starling
#March 2017

#================================================================
# Math Test: Gibbs Sampler ======================================
#================================================================

gibbs.mathtest = function(y,x,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Samler for math test hierarchical normal model.
	#			Unknown parameters: (tau.sq, sig.sq, mu)	
	#-------------------------------------------------------------
	#MODEL:		Hierarchical model:
	#			y_ij ~ N(theta_i, sig.sq)
	#			theta_i ~ N(mu, tau.sq * sig.sq)
	#			mu ~ I(mu) ie indicator for real line (flat prior)
	#			sig.sq ~ 1/sig.sq * I(sig.sq) ie jeffreys prior
	#			tau.sq ~ I(tau.sq) ie indicator for pos real line (flat prior)
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of y_ij math test scores
	#			x = vector of schools
	#			iter = number of posterior samples to generate.
	#			burn = number of burn-in values to discard.
	#			thin = thinning to reduce autocorr of chain.
	#-------------------------------------------------------------
	#OUTPUTS:	theta.post = matrix of posterior theta_i samples (rows = samples)
	#			mu.post = vector of posterior mu samples
	#			sig.sq.post = vector of posterior sig2.sq samples
	#			tau.sq.post = vector of posterior tau2.sq samples
	#-------------------------------------------------------------
	
	n = length(y)							#Number of observations.
	ni = aggregate(y, list(x), length)$x	#Sample sizes for each group.
	p = length(unique(x))					#Number of groups.
	ybar = aggregate(y, list(x), mean)$x	#Sample group means.
	
	#Set up data structures to hold Gibbs samples from posteriors.
	theta 	= matrix(0,iter,p)	#Each row is a p-length Gibbs sample.
	mu 	= rep(0,iter)				#Vector of posterior mu samples
	sig.sq = rep(0,iter)			#Vector of posterior sig2.sq samples
	tau.sq = rep(0,iter)			#Vector of posterior tau2.sq samples

	#Initialize each elemens of chain.
	theta[1,] = rep(0,p)
	mu[1] = mean(y)
	sig.sq[1] = 1
	tau.sq[1] = 1
	
	#Iterate through sampler.
	for (i in 2:iter){
		
		#Update theta_i values (for theta_1...theta_p).
		theta.vars = sig.sq[i-1] * (tau.sq[i-1] / (ni * tau.sq[i-1] + 1))
		theta.means = (tau.sq[i-1] / (ni * tau.sq[i-1] + 1)) * (ni * ybar + (1/tau.sq[i-1]) * mu[i-1])
		theta[i,] = rnorm(p,theta.means,theta.vars)
		
		#Update mu.
		mu.var = sig.sq[i-1] * tau.sq[i-1] / p
		mu.mean = mean(theta[i,])
		mu[i] = rnorm(1,mu.mean,mu.var)
		
		#Update sig.sq.
		alpha = (n+p)/2
		ss1 = sum((y - theta[i,x])^2) 		#sum_{i=1}^{p} sum_{j=1}^{ni} (y_ij - theta_i)^2 
		ss2 = sum((theta[i,]-mu[i])^2)		#sum_{i=1}^{p} (theta_i-mu)^2
		beta = ss1/2 + ss2/(2*tau.sq[i-1]) 
		sig.sq[i] = 1 / rgamma(1, alpha, beta)
		
		#Update tau.sq.
		alpha = p/2-1
		beta = (1 / (2*sig.sq[i])) * sum((theta[i,]-mu[i])^2)
		tau.sq[i] = 1 / rgamma(1, alpha, beta)
	}
	
	#Burn beginning observations.
	if (burn > 0){
		theta = theta[-burn,]
		mu = mu[-burn]
		sig.sq = sig.sq[-burn]
		tau.sq = tau.sq[-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		theta = theta[seq(1,nrow(theta),by=thin),]
		mu = mu[seq(1,length(mu),by=thin)]
		sig.sq = sig.sq[seq(1,length(sig.sq),by=thin)]
		tau.sq = tau.sq[seq(1,length(tau.sq),by=thin)]
	}
	
	#Return results.
	return(list(theta=theta, mu=mu, sig.sq=sig.sq, tau.sq=tau.sq))
}