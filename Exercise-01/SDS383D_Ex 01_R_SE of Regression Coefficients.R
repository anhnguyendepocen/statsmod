#SDS 383D. Exercise 01.  January 2016.  Jennifer Starling.

###############################################################
### Quantifying Uncertainty: Some Basic Frequentist Ideas   ###
### Linear Regression										###
### Part B													###
###############################################################

### PURPOSE:  Linear Regression: Estimate Standard Error of Coefficients Beta

# Load the mlbench library containing the ozone data.
library(mlbench)

# Load the ozone data.
ozone = data(Ozone, package='mlbench')

# Look at the help file for details
?Ozone

# Scrub the missing values and extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]				#Extract response.
x = as.matrix(ozone[,2:10])	#Extract design matrix.
x = cbind(1,x)				#Add intercept to design matrix.

#Compute estimated beta_hat value.
betahat = solve(t(x) %*% x) %*% t(x) %*% y

#Compute sigma2_hat estimate.
yhat = x %*% betahat
rss = sum((y-yhat)^2)
sig2_hat = rss / (length(y) - length(betahat))

#Plug in sigma2_hat to obtain beta cov matrix estimate.
covbeta = sig2_hat * solve(t(x) %*% x)

#Standard error estimate for beta_j by sqrt of diagonals.
se_beta = sqrt(diag(covbeta))

#Compare estimates to those produced by lm method.
# Now compare to lm
lm1 = lm(y~x-1)	#The "-1" says not to fit an intercept; we hard-coded it.

summary(lm1)			#Display lm model summary.
betacovlm = vcov(lm1)	#Extract cov matrix from lm model object.
sqrt(diag(betacovlm))	#SE estimates for beta_j from lm object.

#Display comparison of results.
round(se_beta,3)				#my result
round(sqrt(diag(betacovlm)),3)	#lm result




