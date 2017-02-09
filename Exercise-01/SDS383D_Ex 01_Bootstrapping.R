#SDS 383D. Exercise 01.  January 2016.  Jennifer Starling.

###############################################################
### Quantifying Uncertainty: Some Basic Frequentist Ideas   ###
### Bootstrapping 											###
### Parts A & B												###
###############################################################

### PURPOSE:  Bootstrap covariance matrix and MLE estimates.

#-------------------------------------------------------------------
### Bootstrapping Part A

# Let Sigma_hat = cov matrix of sampling dist of beta_hat.
# Write a function that estimates Sigma_hat via bootstrap for a given
# design matrix X and response vector y.  Use it to compute Sigma_hat
# for the ozone data and compare to the parametric theory estimate.

#---------------------
beta_cov_boot = function(X,y,B){
	#OVERVIEW: This function generates B bootstrap realizations of the beta
	# least-squares coefficients.  It then computes the covariance matrix of beta
	# by computing the variances and covariances of each vector of bootstrapped beta_j's.
	
	#INPUTS:	X = nxp design matrix, with no intercept col of 1's (function will add these)
	#			y = nx1 response vector
	#			B = number of bootstrap samples.
	#OUTPUTS:	cov_hat = Sigma_hat covariance matrix for betas.
	
	n = nrow(X)		#Number of observations in design matrix.
	p = ncol(X)		#Number of predictors.
	
	#Matrix to hold each beta bootstrap sample.  Each row is a sample. Each col is a beta_j.
	betahat_boot = matrix(0,nrow=B,ncol=p)
	
	#Pre-cache (X'X)^-1 X'
	XtX = crossprod(X)
	XtXinv = solve(XtX)
	
	#Fit model and obtain residuals, e.
	beta_hat = XtXinv %*% t(X) %*% y
	yhat = X %*% beta_hat
	e = y - yhat
	
	#NOTE: Bootstrapping the residuals only, as we want to treat X as fixed.
	for (b in 1:B){
		samps = sample(1:n,n,replace=T)	#Select bootstrap indices.
		e_boot = e[samps]				#Sample residuals.
		y_boot = yhat + e_boot			#Bootstrapped y values.
		
		#Calculate bootstrapped beta coefficients.
		betahat_boot[b,] = XtXinv %*% t(X) %*% y_boot
	}
	
	#Estimate cov matrix using var(beta_i,beta_j) for all cols.  
	beta_hat_cov = cov(betahat_boot)

	return(beta_hat_cov)	
} #END FUNCTION
#---------------------

#Test this out with the Ozone covariance matrix.

# Load the library & data
library(mlbench)
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

#Compute cov matrix using lm.  (-1) means don't fit an intercept; we hard-coded it in X.
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))

#Run bootstrap function.
mybetacov = beta_cov_boot(x,y,B=10000)

#Display bootstrap estimate.
round(mybetacov,2)	


#Display the parametric normal theory estimate.
round(betacovlm,2)

#-------------------------------------------------------------------
### Bootstrapping Part B - 1

mvn_simulate = function(n,mu,Sigma){
	#PURPOSE: Simulates mvn random variables given a mean mu and cov Sigma.
	#This function returns a single X ~ MVN(mu,Sigma) realization.

	#ALGORITHM:
	# This algorithm is based on the derivation in Multivariate Normal Part E, which
	# showed that any MVN can be written as a linear combo of standard normals:

	# Simulate Z, a vector of p normal random variables.  (p = desired mvn dimension, p = length(mu))
	# Let mu be the vector of means, Sigma be the specified covariance matrix.
	# Let Sigma = LL^T.  
	# Simulate mvn rv as X = LZ + mu
	
	#INPUTS: 	mu = desired vector of means.  Must be length p.
	#			Sigma = desired covariance matrix.  Must be (pxp), symmetric, pos semidef.
	#			n = number of mvn random variables to generate.
	#OUTPUTS: 	x = a matrix of realizations from MVN(mu,Sigma).  (Each x is a column.)
	
	p = length(mu)		#Set length of mu vector.
	z = matrix(rnorm(n*p,0,1),nrow=p,ncol=n)	#Generate p iid standard normals z_i for each realization.
	
	#Compute L using spectral value decomposition.  V %*% lam %*% solve(V)
	#(See notes below.  Cholesky is 3x faster, spectral is more stable.)
	
	eg = eigen(Sigma)			#Store spectoral value decomposition of Sigma.
	V = eg$vectors				#Extract eigen vectors.
	lam = diag(eg$values)		#Extract diagonal matrix of eigenvalues.
	
	L = V %*% sqrt(lam)			#Assign L so LL^T = Sigma	

	#Compute realizations of x ~ mvn(mu,Sigma)
	#x = L %*% z + mu
	x = apply(z,2,function(a) L %*% a + mu)

	return(x)
}

#Test it out:
mu = c(2,5)
Sigma = matrix(c(10,3,3,2),2,2)
x = mvn_simulate(10000,mu,Sigma)

rowMeans(x)

Sigma_hat = matrix(c(
var(x[1,]),
cov(x[1,],x[2,]),
cov(x[1,],x[2,]),
var(x[2,])),byrow=T,nrow=2)

Sigma_hat

#-------------
#NOTE: Doing LL^T decomposition using eigen() spectral decomposition,
#instead of cholesky.  Cholesky is faster, but eigen is more stable.

#Will show how it is done here both ways.
Sigma <- matrix(c(10,3,3,2),2,2)

#1. Cholesky:
t(chol(Sigma)) %*% chol(Sigma)	#equals Sigma
L = t(chol(Sigma))			#Assign L so LL^T = Sigma
L %*% t(L)					#Verify LL^T = Sigma

#2. Spectral Decomposition: A = V diag(lambda) V^(-1)
eg = eigen(Sigma)			#Store spectoral value decomposition of Sigma.
V = eg$vectors				#Extract eigen vectors.
lam = diag(eg$values)			#Extract diagonal matrix of eigenvalues.
V %*% lam %*% solve(V)		#Check reproducing Sigma.

L = V %*% sqrt(lam)			#Assign L so LL^T = Sigma	
L %*% t(L)					#Verify LL^T = Sigma
#-------------

#-------------------------------------------------------------------
### Bootstrapping Part B - 2

## 2. For a given sample x_1,...,x_n from a mvn, estimate the mean 
## vector and covariance matrix from maximum likelihood.

mvn_mle_est = function(x){
	#PURPOSE: For a sample of x1...xn ~ mvn(mu,Sigma),
	#				estimate mu and Sigma using mle estimates.
	
	#INPUT:		x = matrix of values.  Each row must be a sample.  
	#				Each col must be an x_i ~ mvn, s.t. x is (pxn).
	#OUTPUT:	mu_hat = mle estimate of mvn mean
	#			Sig_hat = mle estimate of mvn cov matrix
	
	#Note: MLE estimates derived and calculated.
	#mu_hat = colMeans(x)
	#sigma_hat = cov(x)
	require(mvnmle)
	
	est = mlest(x)
	return(list(muhat=est$muhat,sigmahat=est$sigmahat,loglik=est$value))	
}

# Generate some simulated data to work with.
library(MASS)
n = 100
mu = c(4,7)
Sigma <- matrix(c(10,3,3,2),2,2)
x = mvrnorm(n=n,mu=mu,Sigma=Sigma)

mvn_mle_est(x)	#Output MLE estimates mu_hat, Sigma_hat.
mu				#Output true mu.
Sigma			#Output true Sigma

#-------------------------------------------------------------------
### Bootstrapping Part B - 3

## 3. Bootstrap a given sample x_1...x_n to estimate the sampling 
## distribution of the MLE.

# Generate some simulated data to work with.
library(MASS)
n = 10
mu = c(4,7)
p = length(mu)
Sigma <- matrix(c(10,3,3,2),2,2)
x = mvrnorm(n=n,mu=mu,Sigma=Sigma)

B = 10000							#Number of bootstrap samples.
mu_boot = array(0,c(p,B))		#px1 array to hold mu vectors.
sig_boot = array(0,c(p,p,B))	#pxpxB array to hold cov matrices.

for (b in 1:B){
	#Rows of x to sample.
	rows_boot = sample(1:n,size=n,replace=T)
	
	#Save bootstrap sample of x.
	xb = x[rows_boot,]	
	
	#Save bootstrap calculations for mu and Sigma.
	mu_boot[,b] = colMeans(xb)
	sig_boot[,,b] = cov(xb)
}

mu_boot_mean = apply(mu_boot,1,mean)
sig_boot_mean = apply(sig_boot,c(1,2),mean)

#Plot outputs of each parameter to observe sampling distribution.
#d=2, so two mu params, and four sigma params.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-01/R Files/boot_mu.pdf')
par(mfrow=c(1,2))
hist(mu_boot[1,])
hist(mu_boot[2,])
dev.off()

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-01/R Files/boot_sigma.pdf')
par(mfrow=c(2,2))
hist(sig_boot[1,1,])
hist(sig_boot[1,2,])
hist(sig_boot[2,1,])
hist(sig_boot[2,2,])
dev.off()