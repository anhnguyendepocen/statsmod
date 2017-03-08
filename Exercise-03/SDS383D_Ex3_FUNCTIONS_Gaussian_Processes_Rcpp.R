#SDS 383D - Exercise 3
#Gaussian Process Functions
#Jennifer Starling
#Feb 2017
library(Rcpp)

#================================================================
# Gaussian Processes ============================================
#================================================================

gaussian_process = function(x,mu,cov.fun,params){
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
	#			cov.fun = covariance matrix function.
	#-------------------------------------------------------------
	#OUTPUTS:	fx = vector of realizations from gaussian process.
	#-------------------------------------------------------------
	n = length(x)
	
	#Generate covariance matrix.
	cov = cov.fun(x,x,kappa=params)
	
	#Generate realizations f(x1)...f(xn).
	#Require the mvtnorm package for random normal generation.
	require(mvtnorm)
	fx = rmvnorm(1,mu,cov,method='chol')
	return(as.numeric(fx))
}

gp.predict = function(x,y,x.new,mu,cov.fun,params,sig2=0){
	#-------------------------------------------------------------
	#FUNCTION: 	Generates predictions from the noisy Gaussian Process 
	#			with specified mean and covariance matrix.
	#			Model: y = f(x) + e, with e ~ N(0,sig2*I)
	#-------------------------------------------------------------
	#INPUTS: 	x = vector (x1,...,xn) 
	#			y = observed GP values at each x.  Can be noisy, or not.
	#			params = vector(b,tau1.sq,tau2.sq) of 3 hyperparameters,
	#				where:
	#				b = length (of period)
	#				tau1.sq = variance of function
	#				tau2.sq = noise
	#			mu = vector of means, length n.
	#			cov.fun = covariance matrix function.
	#			sig2 = variance for noise. 0 predicts for a non-noisy GP.
	#-------------------------------------------------------------
	#OUTPUTS:	post.mean = Posterior mean E(f.new|y) 
	#			post.var = Posterior variance Var(f.new|y)
	#-------------------------------------------------------------
	if(class(x)=='matrix'){
		n = nrow(x)
		n.new = nrow(x.new)
	}else{	
		n = length(x)
		n.new = length(x.new)
	}
	
	#Set up partitioned cov matrices.
	C	= cov.fun(x,x,kappa=params)
	Cx	= cov.fun(x,x.new,kappa=params)
	CxT	= t(Cx)
	Cxx	= cov.fun(x.new,x.new,kappa=params)
	
	#Add noise matrix.  (Will be zeros if sig2=0.)
	noise = sig2 * diag(n)
	
	#Calculate posterior means and vars for each predicted value.
	post.mean = CxT %*% solve(C + noise) %*% y
	post.var = diag(  Cxx - CxT %*% solve(C + noise) %*% Cx  )
	
	return(list(post.mean=post.mean,post.var=post.var))	
}

gp.logl.y = function(x,y,mu,cov.fun,params,sig2=0){
	#-------------------------------------------------------------
	#FUNCTION: 	Generates values of the marginal log-likelihood p(y) 
	#			from the noisy Gaussian Process 
	#			with specified mean and covariance matrix.
	#			Model: y = f(x) + e, with e ~ N(0,sig2*I)
	#			Marginal of y: multivariate N(0,sig2I + C)
	#			Logl of Marginal of y: 
	#-------------------------------------------------------------
	#INPUTS: 	x = vector (x1,...,xn) 
	#			y = observed GP values at each x.  Can be noisy, or not.
	#			params = vector(b,tau1.sq,tau2.sq) of 3 hyperparameters,
	#				where:
	#				b = length (of period)
	#				tau1.sq = variance of function
	#				tau2.sq = noise
	#			mu = vector of means, length n.
	#			cov.fun = covariance matrix function.
	#			sig2 = variance for noise. 0 predicts for a non-noisy GP.
	#-------------------------------------------------------------
	#OUTPUTS:	marg.logl = marginal loglikelihood of y (given x, params)
	#-------------------------------------------------------------
	require(mvtnorm)
	n = length(mu)
	
	#Build covariance matrix.
	C	= cov.fun(x,x,kappa=params)
	sig2I = sig2 * diag(n)
	
	#Calculate and return marginal log-likelihood for y.
	marg.logl = dmvnorm(y,mu,sig2I + C,log=T)
	return(marg.logl)
}



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

cppFunction(
	#-------------------------------------------------------------
	#FUNCTION: 	Squared Exponential Covariance Function (non-isotropic distance)
	#			for 2-d x and y matrices.  (Each x, y is a pair of coordinates.)
	#-------------------------------------------------------------
	#INPUTS:	x,y = two matrices containing pairs of coordinates as rows.
	#			kappa = vector(b1,b2,tau1.sq,tau2.sq) of 4 hyperparameters,
	#				b1,b2 		= length (of period) for each dimension
	#				tau1.sq = variance of function
	#				tau2.sq = noise
	#-------------------------------------------------------------
	#OUTPUT:	The squared exponential covariance matrix.
	#-------------------------------------------------------------
	'
	NumericMatrix cov_se_2d(NumericMatrix x, NumericMatrix y, NumericVector kappa) {

  	double arg1, arg2, distance2;

    int n1 = x.nrow();
    int n2 = y.nrow();
    NumericMatrix C(n1,n2);

    for(int i = 0; i < n1; i++) {

    	for(int j=0; j < n2; j++) {
			
    		arg1 = (x(i,0) - y(j,0))/kappa[0];
    		arg2 = (x(i,1) - y(j,1))/kappa[1];
    		distance2 = arg1*arg1 + arg2*arg2;
    		
    			Rcpp::Rcout << "arg1: " << arg1 << std::endl;    			
    			Rcpp::Rcout << "arg2: " << arg2 << std::endl;
    			Rcpp::Rcout << "distance2: " << distance2 << std::endl; 
    			   		
    		C(i,j) = kappa[2]*exp(-0.5*distance2);
    		if(arg1 == 0.0 && arg2 ==0) C(i,j) += kappa[3];
    	}
    }
    return C;
  }
')

cppFunction(
	#-------------------------------------------------------------
	#FUNCTION: 	Matern 5/2 Covariance Function (Euclidean distance)
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
	NumericMatrix cov_m52(NumericVector x, NumericVector y, NumericVector kappa) {

  	double arg, distance2;

    int n1 = x.size();
    int n2 = y.size();

    NumericMatrix C(n1,n2);

    for(int i = 0; i < n1; i++) {

    	for(int j=0; j < n2; j++) {

    		arg = (x[i] - y[j])/kappa[0];
    		distance2 = arg*arg;
    		C(i,j) = kappa[1] * (1 + sqrt(5) * sqrt(distance2) + (5/3) * distance2) * exp(-sqrt(5) * sqrt(distance2));
    		if(arg == 0.0) C(i,j) += kappa[2];
    	}
    }
    return C;
  }
')

cppFunction(
	#-------------------------------------------------------------
	#FUNCTION: 	Matern 5/2 Covariance Function (non-isotropic distance)
	#			for 2-d x and y matrices.  (Each x, y is a pair of coordinates.)
	#-------------------------------------------------------------
	#INPUTS:	x,y = two matrices containing pairs of coordinates as rows.
	#			kappa = vector(b1,b2,tau1.sq,tau2.sq) of 4 hyperparameters,
	#				b1,b2 		= length (of period) for each dimension
	#				tau1.sq = variance of function
	#				tau2.sq = noise
	#-------------------------------------------------------------
	#OUTPUT:	The matern 5/2 covariance matrix.
	#-------------------------------------------------------------
	'
	NumericMatrix cov_se_2d(NumericMatrix x, NumericMatrix y, NumericVector kappa) {

  	double arg1, arg2, distance2;

    int n1 = x.nrow();
    int n2 = y.nrow();
    NumericMatrix C(n1,n2);

    for(int i = 0; i < n1; i++) {

    	for(int j=0; j < n2; j++) {
			
    		arg1 = (x(i,0) - y(j,0))/kappa[0];
    		arg2 = (x(i,1) - y(j,1))/kappa[1];
    		distance2 = arg1*arg1 + arg2*arg2;
    		
    		C(i,j) = kappa[2] * (1 + sqrt(5) * sqrt(distance2) + (5/3) * distance2) * exp(-sqrt(5) * sqrt(distance2));
    		if(arg1 == 0.0 && arg2 ==0) C(i,j) += kappa[3];
    	}
    }
    return C;
  }
')
