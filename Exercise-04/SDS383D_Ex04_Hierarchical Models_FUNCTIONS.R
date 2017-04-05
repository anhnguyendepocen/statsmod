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

#================================================================
# Cheese: Hierarchical Model Gibbs Sampler ======================
#================================================================

gibbs.cheese = function(y,X,idx,muB,V,C,d,iter=11000,burn=1000,thin=2){
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
	require(MASS)								#Generalized inverse.
	
	n = length(y)									#Number of observations.
	ns = aggregate(y, list(idx$store), length)$x	#Sample sizes for each group (# time poinss per store).
	s = length(unique(idx$store))					#Number of groups (stores).
	p = ncol(X)										#Number of predictors.

	#-------------------------------------------------------------
	#Set up yi, Xi, Wi matrices for each store.  Save in a list.
	yi = list()
	Xi = list()
	Wi = list()
	
	#Iterate through each store number.
	for (i in 1:s){	
		yi[[i]] = y[which(idx$store==i)]
		Xi[[i]] = as.matrix(X[which(idx$store==i),])
		Wi[[i]] = Xi[[i]]
	}

	#Set up data structures to hold the Gibbs samples from the posteriors.
	Beta = matrix(0,p,iter)				#Each entry is a (px1) gibbs sample.  (Constant across store.)
	bi = array(0,dim=c(p,s,iter))		#Each gibbs sample slice is a (sxp) matrix, with p covariates for s stores.
	Sigma = array(0,dim=c(p,p,iter))	#Each gibbs sample slice is a (pxp) matrix.
	sig.sq = rep(0,iter)				#Each sig.sq value is a gibbs sample.
		
	#Initialize first element of chain.
	Beta[1] = 0
	bi[,,1] = rep(0,s*p)
	sig.sq[1] = 1
	Sigma[,,1] = diag(p)
	
	#Iterate through sampler.
	for (k in 2:iter){
		
		#print iteration.
		print(k)
		
		#--------------------------------------------------------
		### Update bi for each store.
		
		#Precache values.
		SigInv = solve(Sigma[,,k-1])
		
		SS.Beta.1 = 0	#\sum_{i=1}^{s} Xi^T * Xi for Beta update.
		SS.Beta.2 = 0	#\sum_{i=1}^{s} Xi^T * (yi + Wi*bi) for Beta update.
		SS.sig.sq = 0	#\sum_{i=1}^{s} Xi^T(yi-Wi*bi) for sig.sq update.
		SS.Sigma = 0	#\sum_{i=1}^{s} bi * bi^T for Sigma update.
		
		for (j in 1:s){
			
			#Extract elements for store j.
			X.j = Xi[[j]]
			y.j = yi[[j]]
			W.j = Wi[[j]]
			
			#Precache WTW and y-XB
			WTW = crossprod(W.j)
			YmXB = y.j - X.j %*% Beta[,k-1] 
			
			#Save prior precision and data precision.			
			Var = solve(SigInv + (1/sig.sq[k-1]) * WTW)
			mean = (1/sig.sq[k-1]) * Var %*% crossprod(W.j,YmXB)
			
			bi[,j,k] = rmvnorm(1,mean,Var)
			
			#While in store loop, precache values to be used in Beta, sig.sq, and Sigma updates.
			SS.Beta.1 =+ crossprod(X.j)
			SS.Beta.2 =+ crossprod(X.j, y.j - W.j %*% bi[,j,k])		
			SS.sig.sq =+ crossprod(y.j - X.j %*% Beta[,k-1] - W.j %*% bi[,j,k])
			SS.Sigma =+ tcrossprod(bi[,j,k])
		} #End store loop.
		
		#--------------------------------------------------------
		### Update Beta.
		
		#Precache.
		VInv = solve(V)
		
		Var = solve(VInv + (1/sig.sq[k-1]) * SS.Beta.1)
		mean = Var %*% ( VInv %*% muB + (1/sig.sq[k-1]) * SS.Beta.2)
			
		Beta[,k] = rmvnorm(1,mean,Var)
		
		#--------------------------------------------------------
		### Update sig.sq.
		sig.sq[k] = 1 / rgamma(1,n/2, SS.sig.sq/2)
		
		#--------------------------------------------------------
		### Update Sigma.
		dn = d + s
		Cn = C + SS.Sigma
		
		Sigma[,,k] = riwish(dn,Cn)
	} #End Gibbs sampler loop.
	
	#Burn beginning observations.
	if (burn > 0){
		bi = bi[,,-burn]
		Sigma = Sigma[,,-burn]
		sig.sq = sig.sq[-burn]
	}
	
	#Thin observations.
	if (thin > 0){
		bi = bi[,,seq(1,dim(bi)[3],by=thin)]
		Sigma = Sigma[,,seq(1,dim(Sigma)[3],by=thin)]
		sig.sq = sig.sq[seq(1,length(sig.sq),by=thin)]
	}
	
	#Calculate posterior means.
	bi.pm = apply(bi,c(1,2),mean)
	Beta.pm = rowMeans(Beta)
	sig.sq.pm = mean(sig.sq)
	Sigma.pm = apply(Sigma,c(1,2),mean)
	
	#Return results. 
	return(list(
		bi=bi,
		Beta=Beta,
		sig.sq=sig.sq,
		Sigma=Sigma,
		bi.pm=bi.pm,
		Beta.pm=Beta.pm,
		sig.sq.pm=sig.sq.pm,
		Sigma.pm=Sigma.pm
	))
} #End gibbs.cheese function.

#================================================================
# Polls: Hierarchical Augmented Probit Model Gibbs Sampler ======
#================================================================

gibbs.probit = function(y,X,idx,muB,Sigma,iter=11000,burn=1000,thin=2){
	#-------------------------------------------------------------
	#FUNCTION: 	Gibbs Sampler for augmented probit model.
	#-------------------------------------------------------------
	#MODEL:		Hierarchical model:
	#			y_{it} ~ Bernoulli(p = prob(Phi(mu_i + X_{it}^T %** Beta)) )
	#			y's observed, but we introduce latent z's, where:
	#
	#			z_i ~ N(mu_i + X_i*Beta)
	#			Beta ~ N(muB, V)
	#			mu_i ~ N(0,tau.sq)
	#			tau.sq ~ 1/tau.sq (Jeffreys)
	#-------------------------------------------------------------
	#INPUTS: 	y = vector of log-demand for stores i at times t.
	#			X = matrix of responses, with one row for each store/time combo.
	#			muB,Sigma = hyperparameters for fixed Beta.
	#-------------------------------------------------------------
	#OUTPUTS:	probs.post = vector of posterior p(y=1) samples.
	#			z.post = vector of posterior z samples.
	#			B.post = matrix of posterior B samples (rows = samples)
	#			mui.post = matrix of posterior mui samples (rows = samples)
	#			tau.sq.post = vector of posterior tau.sq samples. 
	#-------------------------------------------------------------
	require(MCMCpack)						#Sample from Inverse Wishart.
	require(mvtnorm)						#Sample from Multivariate Normal.
	
	n = length(y)									#Number of observations.
	ns = aggregate(y, list(idx$state), length)$x	#Sample sizes for each group (# time poinss per store).
	s = length(unique(idx$state))					#Number of groups (stores).
	p = ncol(X)										#Number of predictors.

	#-------------------------------------------------------------
	
	#Set up data structures to hold the Gibbs samples from the posteriors.
	Beta = matrix(0,p,iter)				#Each entry is a (px1) gibbs sample.  (Constant across state.)
	mu.i = matrix(0,s,iter)				#Each slice is a (sxp) matrix, with s state intercept offsets.
	tau.sq = rep(0,iter)				#Each entry is a scalar tau.sq.
	z.i = list()						#Each list element will be a (ni x iter) matrix of z_ij for that state.
	
	#Set up yi, Xi, Wi matrices for each store.  Save in a list.
	yi = list()		#Holds responses for each state.
	Xi = list()		#Design matrix for covariates (fixed effects).
	Wi = list()		#Vector of ones, length ns, for each state.
	z.i = list()	#Matrices of latent variables for each state.
	py.i = list()	#Matrices of P(y_ij=1) for each state.
	
	#Setup: Iterate through each state.
	for (j in 1:s){	
		#Set up y, X and W for each state.
		yi[[j]] = y[which(idx$state==j)]
		Xi[[j]] = as.matrix(X[which(idx$state==j),])
		Wi[[j]] = as.matrix(rep(1,ns[j]),ncol=1,drop=F)
		
		z.i[[j]] = matrix(0,ns[j],iter)
		z.i[[j]][,1] = rnorm(ns[j],0,1)
		
		py.i[[j]] = matrix(0,ns[j],iter)
		py.i[[j]][,1] = rbeta(ns[j],.5,.5)
	}			
	
	#Initialize remaining chains' first elements.
	Beta[,1] = 0
	mu.i[,1] = rep(.01,s)
	tau.sq[1] = 1
	
	#Precache SigInv. (No prior on Sigma.)
	SigInv = solve(Sigma)
	
	#Iterate through sampler.
	for (k in 2:iter){
		
		#print iteration.
		print(k)
		
		#--------------------------------------------------------
		### Update tau.sq.
		rt = sum((mu.i[,k-1])^2)
		tau.sq[k] = 1 / rgamma(1,s/2 + 1,rt/2)
		
		#--------------------------------------------------------
		### Update Beta.
		
		#Precache sums over states for Beta update.
		SS.Beta.1 = 0	#\sum_{j=1}^{s} X.i ^T X.i
		SS.Beta.2 = 0	#\sum_{j=1}^{s} X.j^T * (z.j - mu.j*W.j)
		
		for (j in 1:s){
			#Extract elements for each state.
			z.j = z.i[[j]][,k-1]
			X.j = Xi[[j]]
			W.j = Wi[[j]]
			
			SS.Beta.1 =+ crossprod(X.j)
			SS.Beta.2 =+ crossprod(X.j,z.j - mu.i[j,k-1]*W.j)
		}
		
		#Calculate mean and variance, and draw updated value.
		Var = solve(SigInv + SS.Beta.1)
		mean = Var %*% (SigInv %*% muB + SS.Beta.2)
		Beta[,k] = rmvnorm(1,mean,Var)
		
		#--------------------------------------------------------
		### Iterate through states to update mu.i, z.i, py.i.
		
		#Iterate through states.
		for (j in 1:s){
			
			# Extract elements for each state.
			z.j = z.i[[j]][,k-1]
			X.j = Xi[[j]]
			W.j = Wi[[j]]
			
			#----------------------------------------------------
			# Update mu.i for each state.
			Var = 1 / (tau.sq[k] + ns[j])
			mean = Var * crossprod(W.j, z.j - X.j %*% Beta[,k])
			mu.i[j,k] = rnorm(1,mean,sqrt(Var))
						
			#----------------------------------------------------
			# Update z.j for each state.
			Var = diag(ns[j])
			mean = mu.i[j,k] * W.j + X.j %*% Beta[,k]
			z.i[[j]][,k] = as.numeric(rmvnorm(1,mean,Var))
			
			#----------------------------------------------------
			#Update py.j for each state.
			py.i[[j]][,k] = pnorm(mu.i[j,k]*W.j + X.j %*% Beta[,k])
	
		} #End store loop.
	} #End Gibbs sampler loop.
	
	#Burn beginning observations.
	if (burn > 0){
		mu.i = mu.i[,-burn]
		tau.sq = tau.sq[-burn]	
		z.i = lapply(z.i, function(x) x[,-burn])
		py.i = lapply(py.i, function(x) x[,-burn])
	}
	
	#Thin observations.
	if (thin > 0){
		mu.i = mu.i[,seq(1,ncol(mu.i),by=thin)]
		tau.sq = tau.sq[seq(1,length(tau.sq),by=thin)]
		
		z.i = lapply(z.i, function(x) x[,seq(1,ncol(x),by=thin)])
		py.i = lapply(py.i, function(x) x[,seq(1,ncol(x),by=thin)])
	}
	
	#Calculate posterior means.
	mu.i.pm = rowMeans(mu.i)
	Beta.pm = rowMeans(Beta)
	tau.sq.pm = mean(tau.sq)
	z.i.pm = lapply(z.i, function(x) rowMeans(x))
	py.i.pm = lapply(py.i, function(x) rowMeans(x))

	#Calculate predicted y values based on posterior means of z.i.
	y.pred = unlist(lapply(z.i.pm, function(x) ifelse(x>0,1,0)))
	
	#Return results. 
	return(list(
		mu.i=mu.i,
		Beta=Beta,
		tau.sq=tau.sq,
		z.i = z.i,
		py.i = py.i,
		mu.i.pm=mu.i.pm,
		Beta.pm = Beta.pm,
		tau.sq.pm = tau.sq.pm,
		z.i.pm=z.i.pm,
		py.i.pm=py.i.pm,
		y.pred=y.pred
	))
} #End gibbs.probit function.