#Load Rcpp for Gaussian Processes functions.
library(Rcpp)

#=======================================================================
# Hierarchical Gaussian Processes: Build Full Cov Matrix for Group i ===
#=======================================================================
full.cov = function(Yi,params,log.params=T){
	#-------------------------------------------------------------
	#FUNCTION: 	Calculates the full covariance matrix for a group, Yi.
	#-------------------------------------------------------------
	#INPUTS: 	Yi = data frame for group i, including columns:
	#					group, gene, replicate, time, log2exp
	#			parms = vector of parameters for group.
	#				c(ls.f, var.f, ls.g, var.g, ls.h, var.h, sig.sq)
	#				where ls = lengthscale, var = variance, sig.sq = error variance.
	#			log.params = indicates if input on log-scale.
	#				Use log scale if optimizing by Nelder Mead to ensure non-negs.
	#-------------------------------------------------------------
	#OUTPUTS:	Full covariance matrix for group i.
	#-------------------------------------------------------------
	require(Matrix)		#For bdiag.
	
	#-------------------------------------------------------------
	# Extract parameters and converts parameters from log scale.
	#-------------------------------------------------------------
	if(log.params==T){
		par.f = exp(params[1:2])
		par.g = exp(params[3:4])
		par.h = exp(params[5:6])
		sig.sq = exp(params[7])
	} else{
		par.f  = (params[1:2])
		par.g  = (params[3:4]) 
		par.h  = (params[5:6])
		sig.sq = (params[7])
	}

	#-------------------------------------------------------------
	# Set up sample sizes for each level of hierarchy.
	#-------------------------------------------------------------
	
	Ni = nrow(Yi)						#Number of obs in group.
	reps = unique(Yi$replicate)			#Vector of replicates in group.
	R = length(reps)					#Number of replicates.
	genes = unique(Yi$gene)				#Vector of genes in group.
	G = length(genes)					#Number of genes in group.
	
	#Use first gene in group for Kg block setup.  (All genes have same reps/times.)
	gene.1 = genes[1]					#ID for first gene.
	t.i = Yi$time						#Times for all obs in group.
	t.n = Yi[Yi$gene==gene.1,'time']	#Times for gene 1.
	
	#Use first replicate in first gene for Kf block setup.  (All reps in genes have same times.)
	rep.1 = reps[1]
	t.r = Yi[Yi$gene==gene.1 & Yi$replicate==rep.1,'time']

	#-------------------------------------------------------------
	# Kf mini-blocks assembly to add to the gene block.
	#-------------------------------------------------------------
	Kf.blocks = list()
	
	for (i in 1:R){
		Kf.blocks[[i]] = cov_se(t.r,t.r,c(par.f,0))
	}
		
	Kf = bdiag(Kf.blocks)
	
	#-------------------------------------------------------------
	# Kg blocks assembly to for each gene.
	#-------------------------------------------------------------
	Kg = cov_se(t.n,t.n,c(par.g,0)) 
	
	#-------------------------------------------------------------
	# K.h blocks.
	#-------------------------------------------------------------
	Kh = cov_se(t.i,t.i,c(par.h,0))
	
	#-------------------------------------------------------------
	# Assemble Cov.n matrix for gene
	#-------------------------------------------------------------
	Cov.n = Kg + bdiag(Kf) + diag(sig.sq,length(t.n))
	
	#-------------------------------------------------------------
	# Assemble full cov matrix for group.
	#-------------------------------------------------------------
	Cov = as.matrix(Kh + bdiag(rep(list(Cov.n),G)))
	
	return(Cov)
}

#================================================================
# Marginal Likelihood Function ==================================
#================================================================

marg.neg.ll = function(params,Yi){
	#-------------------------------------------------------------
	#FUNCTION: 	Calculates marginal neg loglhood for group i.
	#-------------------------------------------------------------
	#INPUTS: 	Yi = data frame for group i, including columns:
	#					group, gene, replicate, time, log2exp
	#			parms = vector of parameters for group.
	#				c(ls.f, var.f, ls.g, var.g, ls.h, var.h, sig.sq)
	#				where ls = lengthscale, var = variance, sig.sq = error variance.
	#-------------------------------------------------------------
	#OUTPUTS:	Marginal log-likelihood for group i.
	#-------------------------------------------------------------
		
	#Assemble full cov matrix for group.
	Cov.i = full.cov(Yi, params)
	
	#Marginal neg logl, Yi | Ti ~ N(0,Sigma) where Sigma = full cov matrix.
	m.nll = - dmvnorm(Yi$log2exp, mean=rep(0,nrow(Yi)), sigma=Cov.i,log=T)

	return(m.nll)
}	


#================================================================
# Optimize Parameters over Marginal Likelihood for all groups ===
#================================================================

opt.param = function(Y,init.par){
	#-------------------------------------------------------------
	#FUNCTION: 	Calculates optimal parameters for each group i,
	#			using maximum likelihood estimation via Nedler Mead.
	#-------------------------------------------------------------
	#INPUTS: 	Y = data frame for all groups, including columns:
	#					group, gene, replicate, time, log2exp
	#			init.par = vector of initial parameters for group.
	#				c(ls.f, var.f, ls.g, var.g, ls.h, var.h, sig.sq)
	#				where ls = lengthscale, var = variance, sig.sq = error variance.
	#-------------------------------------------------------------
	#OUTPUTS:	Matrix of optimal parameters; each column is a group.
	#-------------------------------------------------------------

	#Data structure.	
	opt.pars = matrix(NA,7,3)
	colnames(opt.pars) = c('group1','group2','group3')
	rownames(opt.pars) = c('ls.f','var.f','ls.g','var.g','ls.h','var.h','sig.sq')
	
	#Extracts groups in data set.
	G = unique(Y$group)
	
	#Loop over groups and optimize for each group.
	for (g in G){
		
		#Select data for group.
		Yi = Y[Y$group==g,]
		
		#Optimize using Nelder Mead.
		opt.pars[,g] = optim(par=init.par, fn=marg.neg.ll, Yi=Yi)$par
	}

	#For group 1, capping alpha.h value at 10, since ll is very flat.  (Credit: Spencer)
	opt.pars[6,1] = min(opt.pars[6,1],10)	
	
	#Returns log-scale params.
	return(opt.pars)
}

#================================================================
# Gaussian Processes ============================================
#================================================================

gaussian.process = function(x,mu,cov_se,params){
	#-------------------------------------------------------------
	#FUNCTION: 	Generates realizations from the Gaussian Process 
	#			with specified mean and covariance matrix.	
	#-------------------------------------------------------------
	#INPUTS: 	x = vector (x1,...,xn) or matrix of 2-d x values in rows
	#			params = vector(b,tau1.sq,tau2.sq) of 3 hyperparameters,
	#				where:
	#				b = length (of period)
	#				tau1.sq = variance of function
	#				tau2.sq = noise
	#			mu = vector of means, length n.
	#			cov_se = covariance matrix function.
	#-------------------------------------------------------------
	#OUTPUTS:	fx = vector of realizations from gaussian process.
	#-------------------------------------------------------------
	n = length(x)
	
	#Generate covariance matrix.
	cov = cov_se(x,x,kappa=params)
	
	#Generate realizations f(x1)...f(xn).
	#Require the mvtnorm package for random normal generation.
	require(mvtnorm)
	fx = rmvnorm(1,mu,cov,method='chol')
	return(as.numeric(fx))
}

#================================================================
# Squared Exponential Covariance Function =======================
#================================================================

cppFunction(
	#-------------------------------------------------------------
	#FUNCTION: 	Squared Exponential Covariance Function (Euclidean distance)
	#			for 1-d x and y vectors.
	#-------------------------------------------------------------
	#INPUTS:	x,y = two vectors from the same space.
	#			kappa = vector(b,tau1.sq,tau2.sq) of 3 hyperparameters,
	#				b 		= length (of period)
	#				tau1.sq = variance of function
	#				tau2.sq = noise
	#-------------------------------------------------------------
	#OUTPUT:	The squared exponential covariance matrix.
	#-------------------------------------------------------------
	'
	NumericMatrix cov_se(NumericVector x, NumericVector y, NumericVector kappa) {

  	double arg, distance2;

    int n1 = x.size();
    int n2 = y.size();

    NumericMatrix C(n1,n2);

    for(int i = 0; i < n1; i++) {

    	for(int j=0; j < n2; j++) {

    		arg = (x[i] - y[j])/kappa[0];
    		distance2 = arg*arg;
    		C(i,j) = kappa[1]*exp(-0.5*distance2);
    		if(arg == 0.0) C(i,j) += kappa[2];
    	}
    }
    return C;
  }
')

#================================================================
# Predict Genes, Group and Cluster Profiles =====================
#================================================================

hgp.predict = function(Y,params,t.pred){
	#-------------------------------------------------------------
	#FUNCTION: 	Predicts posterior group profiles, gene profiles, and 
	#			gene-replicate profiles.
	#-------------------------------------------------------------
	#INPUTS: 	Y = Data frame for all groups, including columns:
	#					group, gene, replicate, time, log2exp
	#			params = Matrix of parameters for each group.
	#				groups are in columns, where each col consists of:
	#				c(ls.f, var.f, ls.g, var.g, ls.h, var.h, sig.sq)
	#				where ls = lengthscale, var = variance, sig.sq = error variance.
	#			t.pred = Times for predicting new group/gene/rep values.
	#-------------------------------------------------------------
	#OUTPUTS:	A list including:
	#			Y 		= Data set input by user.
	#			params 	= Parameters input by user.
	#			f.mean	= Matrix, cols are posterior means for each gene-rep profile at times t.pred.
	#			f.var	= Matrix, cols are posterior means for each gene-rep profile at times t.pred.
	#			f.lb	= Matrix, cols are posterior CI lower bounds for each gene-rep profile at times t.pred.
	#			f.ub	= Matrix, cols are posterior CI upper bounds for each gene-rep profile at times t.pred.
	#			g.mean	= Matrix, cols are posterior means for each gene profile at times t.pred.
	#			g.var	= Matrix, cols are posterior vars for each gene profile at times t.pred.
	#			g.lb	= Matrix, cols are posterior CI lower bounds for each gene profile at times t.pred.
	#			g.ub	= Matrix, cols are posterior CI upper bounds for each group profile at times t.pred.
	#			h.mean	= Matrix, cols are posterior means for each group profile at times t.pred.
	#			h.var	= Matrix, cols are posterior vars for each group profile at times t.pred.
	#			h.lb	= Matrix, cols are posterior CI lower bounds for each gene profile at times t.pred.
	#			h.ub	= Matrix, cols are posterior CI upper bounds for each ggroup profile at times t.pred.
	#-------------------------------------------------------------

	N = nrow(Y)							#Total number of obs in data set.
	
	groups = unique(Y$group)			#Vector of groups.
	H = length(groups)					#Number of groups.
	
	genes = unique(Y$gene)				#Vector of genes.
	G = length(genes)					#Number of genes.
	
	reps = unique(Y$replicate)			#Vector of reps.
	R = length(reps)					#Number of reps.

	### Set up times.
	T1 = length(Y$time)
	T2 = length(t.pred)
	
	#-------------------------------------------------------------
	# Initialize data structures to hold predictions.
	#-------------------------------------------------------------
	
	#Holds predicted values for each group in cols.
	h.mean = matrix(NA,T2,H)	#Holds posterior means for each gene profile.
	h.lb   = matrix(NA,T2,H)	#Holds 95% CI lower bound.
	h.ub   = matrix(NA,T2,H)	#Holds 95% CI upper bound.
	h.var  = matrix(NA,T2,H)	#Holds covariance matrix diagonals, ie posterior variances for each group profile.
	
	#Holds predicted values for each group in cols.
	g.mean = matrix(NA,T2,G)	#Holds posterior means for each gene profile.
	g.lb   = matrix(NA,T2,G)	#Holds 95% CI lower bound.
	g.ub   = matrix(NA,T2,G)	#Holds 95% CI upper bound.
	g.var  = matrix(NA,T2,G)	#Holds covariance matrix diagonals, ie posterior variances for ea gene profile.
	
	#Holds predicted values for each rep in each group/gene.
	f.mean = matrix(NA,T2,G*R)	#Holds posterior means for each gene/rep profile.
	f.lb   = matrix(NA,T2,G*R)	#Holds 95% CI lower bound.
	f.ub   = matrix(NA,T2,G*R)	#Holds 95% CI upper bound.
	f.var  = matrix(NA,T2,G*R)	#Holds covariance matrix diagonals, ie posterior variances for ea gene/rep profile.
	
	colnames(g.mean) = colnames(g.lb) = colnames(g.ub) = colnames(g.var) = genes
	colnames(f.mean) = colnames(f.lb) = colnames(f.ub) = colnames(f.var) = paste(rep(genes,each=R),rep(reps,G))
	colnames(h.mean) = colnames(h.lb) = colnames(h.ub) = colnames(h.var) = groups
	
	#-------------------------------------------------------------
	# Predict Group Profiles (hi* | y)
	#-------------------------------------------------------------
	
	#Loop over groups.
	for (i in 1:H){
		
		print(paste('Group',i))
		
		#Group name.
		grp = groups[i]
		
		#Data for group.
		Yi = Y[Y$group==grp,]

		#Set up upper-left cov matrix, Cov(Yi).
		Cov.Yi = full.cov(Yi,params[,i],log=T)
		Cov.Yi.Inv = solve(Cov.Yi + diag(1e-6,nrow(Yi)))
		
		#Cov fctn parameters for group (non-log-scale)
		pars.i = c(exp(params[5:7,i]), 0)
		
		#Set up lower-right cov matrix, marginal of hi at t.pred.
		Kh.xx = cov_se(t.pred,t.pred,pars.i)
		
		#Set up off-diagonal blocks for data-hi.star covars.
		t = Y[Y$group==grp,'time']
		Kh.x = cov_se(t.pred,t,pars.i)
		Kh.xT = t(Kh.x)
		
		#Posterior mean for h function for group.
		h.mean[,i] = Kh.x %*% Cov.Yi.Inv %*% Yi$log2exp
		
		#Posterior covar matrix of h for group.
		h.cov = Kh.xx - Kh.x %*% Cov.Yi.Inv %*% Kh.xT
		
		#Posterior variance for groups, and 95% CI bounds.
		h.var[,i] = diag(h.cov)
		h.lb[,i] = h.mean[,i] - 1.96 * sqrt(h.var[,i])
		h.ub[,i] = h.mean[,i] + 1.96 * sqrt(h.var[,i])
		
		#-------------------------------------------------------------
		# For each group: Predict Gene Profiles (gn* | y) 
		#-------------------------------------------------------------
		
		#Genes in group. 
		genes.i = unique(Yi$gene)
		
		#Loop over gene profiles in group.
		for (j in 1:length(genes.i)){
		
			print(paste('Group',i, 'Gene',j))
			
			#Gene name.
			gene = genes.i[j]
			
			#Data for gene in group.
			Yn = Yi[Yi$gene==gene,]
			
			#Cov fctn parameters for gene (non-log-scale)
			pars.n = c(exp(params[3:4,i]),0)

			#Set up lower-right cov matrix, marginal of gn at t.pred.
			Kg.xx = Kh.xx + cov_se(t.pred,t.pred,pars.n) 
		
			#Update off-diagonal blocks for data-gn.star covars.
			tn = Yn$time						#All times for gene.
			Tn = length(tn)						#Length of times vector for gene.
			gn.chg.idx = ((j-1) * Tn + 1):(j*Tn) #Indices in off-diag cov matrix which need upating.
			
			Kg.x.small = cov_se(t.pred,tn,pars.n)
			Kg.x = Kh.x
			Kg.x[ ,gn.chg.idx] = Kg.x[,gn.chg.idx] + Kg.x.small
			Kg.xT = t(Kg.x)

			#Posterior mean for g function for gene.
			g.mean[,gene] = Kg.x %*% Cov.Yi.Inv %*% Yi$log2exp
		
			#Posterior covar matrix of g for gene.
			g.cov = Kg.xx - Kg.x %*% Cov.Yi.Inv %*% Kg.xT
		
			#Posterior variance for genes, and 95% CI bounds.
			g.var[,gene] = diag(g.cov)
			g.lb[,gene] = g.mean[,gene] - 1.96 * sqrt(g.var[,gene])
			g.ub[,gene] = g.mean[,gene] + 1.96 * sqrt(g.var[,gene])	
			
			#-------------------------------------------------------------
			# Predict Rep Profiles For Each Gene (fnr* | y)
			#-------------------------------------------------------------
			for (k in reps){
				
				print(paste('Group',i, 'Gene',j, 'Rep',k))
				
				#Identify rep.
				rep = k
				
				#Gene-rep name for saving into data structures.
				gene.rep = paste(gene,rep)
				
				#Cov fctn parameters for gene (non-log-scale)
				pars.r = c(exp(params[1:2,i]),0)

				#Set up marginal cov matrix for rep.
				Kf.xx = Kg.xx + cov_se(t.pred,t.pred,pars.r) 
				
				#Update off-diagonal blocks for data-fnr.star covars.
				tnr = Yn[Yn$replicate==rep,'time']	#All times for rep in gene.
				Tnr = length(tnr)					#Length of times vector for rep in gene.
				
				rep.chg.idx = tnr + (Tnr * (match(k,reps)-1)) #Indices in off-diag cov matrix which need upating.
				
				Kf.x.small = cov_se(t.pred,tnr,pars.r)
				Kf.x = Kg.x
				Kf.x[, gn.chg.idx[rep.chg.idx]] = Kf.x[, gn.chg.idx[rep.chg.idx]] + Kf.x.small
				Kf.xT = t(Kf.x)
				
				#Posterior mean for f function for rep in gene.
				f.mean[,gene.rep] = Kf.x %*% Cov.Yi.Inv %*% Yi$log2exp
		
				#Posterior covar matrix of f for rep in gene.
				f.cov = Kf.xx - Kf.x %*% Cov.Yi.Inv %*% Kf.xT
		
				#Posterior variance for reps in genes, and 95% CI bounds.
				f.var[,gene.rep] = diag(g.cov)
				f.lb[,gene.rep] = f.mean[,gene.rep] - 1.96 * sqrt(f.var[,gene.rep])
				f.ub[,gene.rep] = f.mean[,gene.rep] + 1.96 * sqrt(f.var[,gene.rep])	
	
			} #End reps in genes in groups loop.
			
		} #End genes in groups loop.
	} #End groups loop.

	#-------------------------------------------------------------
	### Return output.
	return(list(
		Y=Y,
		params=params,
		f.mean=f.mean, f.var=f.var, f.lb=f.lb, f.ub=f.ub,
		g.mean=g.mean, g.var=g.var, g.lb=g.lb, g.ub=g.ub,
		h.mean=h.mean, h.var=h.var, h.lb=h.lb, h.ub=h.ub
	))
	
} #End function.
