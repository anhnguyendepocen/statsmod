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
	#OUTPUTS:	py.i = vector of posterior p(y=1) samples for each state (list).
	#			z.i = vector of posterior z samples for each state (list).
	#			Beta = matrix of posterior Beta samples (rows = samples)
	#			mu.i = matrix of posterior mui samples (rows = samples)
	#			tau.sq = vector of posterior tau.sq samples. 
	#-------------------------------------------------------------
	require(mvtnorm)						#Sample from Multivariate Normal.
	
	n = length(y)									#Number of observations.
	ns = aggregate(y, list(idx$state), length)$x	#Sample sizes for each group (# time poinss per store).
	s = length(unique(idx$state))					#Number of groups (stores).
	p = ncol(X)										#Number of predictors.
	SigInv = solve(Sigma)							#Precache SigInv. (No prior on Sigma.)
	
	#-------------------------------------------------------------
	#Set up data structures to hold the Gibbs samples from the posteriors.
	Beta = matrix(0,p,iter)				#Each entry is a (px1) gibbs sample.  (Constant across state.)
	mu.i = matrix(0,s,iter)				#Each slice is a (sxp) matrix, with s state intercept offsets.
	tau.sq = rep(0,iter)				#Each entry is a scalar tau.sq.
	z.i = list()						#Each list element will be a (ni x iter) matrix of z_ij for that state.
	
	#Set up lists to hold various values by store, since sample sizes vary.
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
	
	#-------------------------------------------------------------
	#Iterate through sampler.
	for (k in 2:iter){
		
		#Progress report.
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

		for (j in 1:s){
			
			# Extract elements for each state.
			X.j = Xi[[j]]
			W.j = Wi[[j]]
			y.j = yi[[j]]
			z.j = z.i[[j]][,k-1]
			
			#----------------------------------------------------
			# Update mu.i for each state.
			Var = 1 / (tau.sq[k] + ns[j])
			mean = Var * crossprod(W.j, z.j - X.j %*% Beta[,k])
			mu.i[j,k] = rnorm(1,mean,sqrt(Var))
						
			#----------------------------------------------------
			# Update z.j for each state.
			mean = mu.i[j,k] * W.j + X.j %*% Beta[,k]
			
			#Set truncation limits based on observed y.i values.
			#If y=1, sample from truncated normal (0,inf).
			#If y=0, sample from truncated normal (-inf,0).
			a = ifelse(y.j==1,0,-Inf)
			b = ifelse(y.j==1,Inf,0)
			
			#Sample z.i as (ns[j] x 1) vector from indep truncated normals.
			z.i[[j]][,k] = rtruncnorm(ns[j],a,b,mean,1)
			
			#----------------------------------------------------
			#Update py.j for each state.
			py.i[[j]][,k] = pnorm(W.j*mu.i[j,k] + X.j %*% Beta[,k])
	
		} #End states loop.
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