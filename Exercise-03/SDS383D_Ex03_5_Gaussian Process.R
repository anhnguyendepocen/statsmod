#Stats Modeling 2
#Exercise 3
#Gaussian Processes - Basics - Part A

#================================================================
# Environment Setup & Data Load =================================
#================================================================

#Housekeeping.
rm(list=ls())

#Load functions.
source('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/RCode/SDS383D_Ex3_FUNCTIONS_Gaussian_Processes_Rcpp.R')

#================================================================
# Test Gaussian Process Function ================================
#================================================================

### Generate gaussian process realizations.
n=100
x = rnorm(100,0,1)

b = 1
tau1.sq = .1
tau2.sq = 1e-6
triplet = c(b,tau1.sq,tau2.sq)

x.se = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_se)
x.m52 = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_m52)

#Test plot.
idx = sort(x, index.return = T)$ix
plot(x[idx],x.se[idx],col='black',lwd=1,type='l',ylim=c(-3,3))

#Plot a bunch.
color=rainbow(10)
for (i in 1:10){
	x.se = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_se)
	lines(x[idx],x.se[idx],lwd=1,col=color[i])
}

#================================================================
# Hyperparameter Plotting for Sq Exp Cov Fctn GP ================
#================================================================

#----------------------------------------------------------------
### Test 1: Varying b for Sq Exponential and Matern 5/2 Cov Functions.

tau1.sq = .001
tau2.sq = 1e-6

B=c(.05,.1,.25,1,4)
colors = rainbow(length(B))

#Squared exponential plot.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/GP_Varying_b_SqExp.pdf')
for (i in 1:length(B)){
	b = B[i]
	triplet = c(b,tau1.sq,tau2.sq)

	fx = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_se)
	idx = sort(x, index.return = T)$ix
	
	if (i==1){
		plot(x[idx],fx[idx],lwd=2,type='l',col=colors[i],ylim=c(-.15,.15),main='Squared Exponential Covariance (tau1.sq=.001,tau2.sq=1e-6)')
	} else{
		lines(x[idx],fx[idx],lwd=2,type='l',col=colors[i])
	}
}
legend('topleft',lty=1,lwd=2,legend=paste('b=',B),col=colors)
dev.off()

#Matern 5/2 plot.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/GP_Varying_b_M52.pdf')
for (i in 1:length(B)){
	b = B[i]
	triplet = c(b,tau1.sq,tau2.sq)

	fx = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_m52)
	idx = sort(x, index.return = T)$ix
	
	if (i==1){
		plot(x[idx],fx[idx],lwd=2,type='l',col=colors[i],ylim=c(-.15,.15),main='Matern 5/2 Covariance (tau1.sq=.001,tau2.sq=1e-6)')
	} else{
		lines(x[idx],fx[idx],lwd=2,type='l',col=colors[i])
	}
}
legend('topleft',lty=1,lwd=2,legend=paste('b=',B),col=colors)
dev.off()

#----------------------------------------------------------------
### Test 2: Varying tau1.sq for Sq Exponential and Matern 5/2 Cov Functions.

b = .1
tau1.list = c(.05,.1,.25,1,4)
tau2.sq = 1e-6
colors = rainbow(length(B))

#Squared exponential plot.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/GP_Varying_tau1sq_SqExp.pdf')
for (i in 1:length(tau1.list)){
	tau1.sq = tau1.list[i]
	triplet = c(b,tau1.sq,tau2.sq)

	fx = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_se)
	idx = sort(x, index.return = T)$ix
	
	if (i==1){
		plot(x[idx],fx[idx],lwd=2,type='l',col=colors[i],ylim=c(-5,5),main='Squared Exponential Covariance (b=.1,tau2.sq=1e-6)')
	} else{
		lines(x[idx],fx[idx],lwd=2,type='l',col=colors[i])
	}
}
legend('topleft',lty=1,lwd=2,legend=paste('tau1.sq=',tau1.list),col=colors)
dev.off()

#Matern 5/2 plot.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/GP_Varying_tau1sq_M52.pdf')
for (i in 1:length(B)){
	tau1.sq = tau1.list[i]
	triplet = c(b,tau1.sq,tau2.sq)

	fx = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_m52)
	idx = sort(x, index.return = T)$ix
	
	if (i==1){
		plot(x[idx],fx[idx],lwd=2,type='l',col=colors[i],ylim=c(-5,5),main='Matern 5/2 Covariance (b=.1,tau2.sq=1e-6)')
	} else{
		lines(x[idx],fx[idx],lwd=2,type='l',col=colors[i])
	}
}
legend('topleft',lty=1,lwd=2,legend=paste('tau1.sq=',tau1.list),col=colors)
dev.off()
#----------------------------------------------------------------
### Test 3: Varying tau2.sq for Sq Exponential and Matern 5/2 Cov Functions.

b = .1
tau1.sq = .01
tau2.list = c(.05,.1,.25,1,4)
colors = rainbow(length(B))

#Squared exponential plot.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/GP_Varying_tau2sq_SqExp.pdf')
for (i in 1:length(tau1.list)){
	tau2.sq = tau2.list[i]
	triplet = c(b,tau1.sq,tau2.sq)

	fx = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_se)
	idx = sort(x, index.return = T)$ix
	
	if (i==1){
		plot(x[idx],fx[idx],lwd=2,type='l',col=colors[i],ylim=c(-5,5),main='Squared Exponential Covariance (b=.1,tau1.sq=.01)')
	} else{
		lines(x[idx],fx[idx],lwd=2,type='l',col=colors[i])
	}
}
legend('topleft',lty=1,lwd=2,legend=paste('tau2.sq=',tau2.list),col=colors)
dev.off()

#Matern 5/2 plot.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/GP_Varying_tau2sq_M52.pdf')
for (i in 1:length(B)){
	tau2.sq = tau2.list[i]
	triplet = c(b,tau1.sq,tau2.sq)

	fx = gaussian_process(x,params=triplet,mu=rep(0,length(x)),cov.fun=cov_m52)
	idx = sort(x, index.return = T)$ix
	
	if (i==1){
		plot(x[idx],fx[idx],lwd=2,type='l',col=colors[i],ylim=c(-5,5),main='Matern 5/2 Covariance (b=.1,tau1.sq=.01)')
	} else{
		lines(x[idx],fx[idx],lwd=2,type='l',col=colors[i])
	}
}
legend('topleft',lty=1,lwd=2,legend=paste('tau2.sq=',tau2.list),col=colors)
dev.off()

