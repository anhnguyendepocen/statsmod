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
			SS.Sigma  =+ tcrossprod(bi[,j,k])
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
		
		Sigma[,,k] = rwishart(dn,Cn)
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