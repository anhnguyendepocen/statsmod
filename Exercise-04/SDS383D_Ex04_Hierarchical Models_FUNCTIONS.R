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

	#Initialize each element of chain.
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

#================================================================
# Cheese: Hierarchical Model Gibbs Sampler ======================
#================================================================

gibbs.cheese = function(y,X,idx,m,v,C,d,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for cheese demand curve model.
	#			Unknown parameters: (m,v,C,d)
	#-------------------------------------------------------------
	#MODEL:		Hierarchical model:
	#			y_{it} ~ N(X_{it} B_i, sig.sq)
	#				where 
	#				y_{it} = logQ_{it} ie log of demand for store i at time t.
	#				X_{it} = [1	logP_{it} ad_{it}] for store i at time t.
	#					where logP is log of price, and ad is indicator for displayed ad.
	#			
	#			B_i ~ N(mu,Sigma)
	#			(mu,Sigma) ~ Normal-Inv-Wishart(m,v,C,d)
	#			sig.sq ~ 1/sig.sq (Jeffrey's Prior)
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of log-demand for stores i at times t.
	#			X = matrix of responses, with one row for each store/time combo.
	#			m,v,C,d = mean, covariance, scale and df hyperpriors for N-I-W.
	#-------------------------------------------------------------
	#OUTPUTS:	B.post = matrix of posterior B samples (rows = samples)
	#			sig.sq.post = vector of posterior sig2.sq samples
	#			muSig.post = vector of posterior (mu,Sigma) IW draws.
	#-------------------------------------------------------------
	require(MCMCpack)							#Sample from Inverse Wishart.
	require(mvtnorm)							#Sample from Multivariate Normal.
	require(LaplacesDemon)						#Sample from Normal-Inverse-Wishart
	
	n = length(y)									#Number of observations.
	nt = aggregate(y, list(idx$store), length)$x	#Sample sizes for each group (# time points per store).
	s = length(unique(idx$store))					#Number of groups (stores).
	p = 3											#Number of predictors.
	ybar = aggregate(y, list(idx$store), mean)$x	#Sample group means.
	
	#Sample group means for each x covariate.
	Xbar = data.frame(matrix(nrow=s,ncol=p))
	for (k in 1:p){
		Xbar[,k] = aggregate(X[,k],list(idx$store),mean)$x
	}
	Xbar = as.matrix(Xbar)
	X = as.matrix(X)

	#Set up data structures to hold Gibbs samples from posteriors.
	Bi = array(0,dim=c(s,p,iter))		#Each slice is a (sxp) matrix, with each column as a vector of covariates for each store. 
	matrix(0,iter,p)					#Each row is a p-length Gibbs sample of B_(px1) vector.
	sig.sq = rep(0,iter)				#Vector of posterior sig2.sq samples
	mu = matrix(0,iter,p)				#Each row is a p-length Gibbs sample of mu_(px1) vector.
	Sigma = array(0,dim=c(p,p,iter))	#Each Sigma[,,1] is a pxp Gibbs sample of Sigma matrix.

	#Initialize each element of chain.
	Bi[,,1] = rep(0,s*p)
	sig.sq[1] = 1
	mu[1,] = rep(mean(y),p)
	Sigma[,,1] = diag(p)
	
	#Iterate through sampler.
	for (k in 2:iter){
		
		#Update B_i values for each store.
		
		SigInv = solve(Sigma[,,k-1])	#Prechache Sigma inverse.
		
		for (j in 1:s){	#Loop through stores.
			Xi = X[which(j==idx$storenum),]	#Pull only Xi values for current store.
			yi = y[which(j==idx$storenum)]	#Pull out responses for current store.
			Bi.var = solve( solve(SigInv) + (nt[j]/sig.sq[k-1]) * t(Xi) %*% Xi	)
			Bi.mean = Bi.var %*% (SigInv %*% mu[k-1,] + (nt[j]/sig.sq[k-1]) * t(Xi) %*% yi )
			
			Bi[j,,k] = rmvnorm(1,Bi.mean,Bi.var)
		}
		
		#Update sig.sq values.
		Bi.lookup = Bi[idx$store,,k]	#Finds Bi for each yi, based on store number.
		SS = sum((y - X %*% t(Bi.lookup))^2)
		sig.sq[k] = 1 / rgamma(1,n/2,.5 * SS)
		
		#Update mu and Sigma values.
		Bi.bar = colMeans(Bi[,,k])						#Precache.
		S = t(Bi[,,k] - Bi.bar)	%*% (Bi[,,k] - Bi.bar)	#Precache.
		
		mn = (v*m + s*Bi.bar) / (v + s)	#Posterior NIW mean parameter.
		vn = v + s						#Posterior cov parameter.
		dn = d + s						#Posterior df parameter.
		Cn = C + S + (v*s / (v+s)) * (Bi.bar - m) %*% t(Bi.bar - m)		#Posterior scale matrix.
		
		niv.draw = rnorminvwishart(n=1,mn,vn,Cn,dn)	#Draw from Normal-Inverse-Wishart.
		Sigma[,,k] = niv.draw$Sigma
		mu[k,] = niv.draw$mu
	}
	
	#Burn beginning observations.
	if (burn > 0){
		Bi = Bi[,,-burn]
		Sigma = Sigma[,,-burn]
		mu = mu[-burn,]
		sig.sq = sig.sq[-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		Bi = Bi[,,seq(1,dim(Bi)[3],by=thin)]
		Sigma = Sigma[,,seq(1,dim(Sigma)[3],by=thin)]
		mu = mu[seq(1,nrow(mu),by=thin),]
		sig.sq = sig.sq[seq(1,length(sig.sq),by=thin)]
	}
	
	#Calculate posterior means.
	Bi.post.mean = apply(Bi,c(1,2),mean)
	Sigma.post.mean = apply(Sigma,c(1,2),mean)
	mu.post.mean = colMeans(mu)
	sig.sq.post.mean = mean(sig.sq)
	
	#Return results.
	return(list(Bi=Bi,Sigma=Sigma,mu=mu,sig.sq=sig.sq,
		Bi.post.mean=Bi.post.mean, 
		Sigma.post.mean=Sigma.post.mean,
		mu.post.mean=mu.post.mean,
		sig.sq.post.mean=sig.sq.post.mean))
}