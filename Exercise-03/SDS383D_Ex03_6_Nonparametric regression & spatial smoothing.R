#Stats Modeling 2
#Exercise 3
#GPs in Nonparametric regression & spatial smoothing


#================================================================
# Environment Setup & Data Load =================================
#================================================================

#Housekeeping.
rm(list=ls())

#Load functions.
source('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/RCode/SDS383D_Ex3_FUNCTIONS.R')

#Load data.
utilities = read.csv('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/statsmod/Course-Data/utilities.csv',header=T)

#Extract data for model.
x = utilities$temp									#average temp.
y = utilities$gasbill / utilities$billingdays		#avg daily bill
n = length(x)

#================================================================
# C: Pointwise Posterior Mean & 95% CI for all observed xi ======
#================================================================

#Set up hyperparameters.
b = 10
tau1.sq = 5
tau2.sq = 0
triplet = c(b,tau1.sq,tau2.sq)
params=triplet

#Run prediction with sigma2=1 to estimate residuals.
pred = gp.predict(x,y,x.new=x,mu=rep(0,n),cov.fun=cov.se,params=triplet,sig2=1)
sig2 = sum(y-pred$post.mean)^2/(n-1)

#Rerun with estimated sigma2.
pred = gp.predict(x,y,x.new=x,mu=rep(0,n),cov.fun=cov.se,params=triplet,sig2=sig2)

#Vectors to hold posterior mean, var, and CI bounds.
post.mean = pred$post.mean
post.se = sqrt(pred$post.var)
post.ci.lb = post.mean - 1.96*post.se
post.ci.ub = post.mean + 1.96*post.se

#Plotting:
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/6_C.pdf')
x.new=x
idx = sort(x.new, index.return = T)$ix
plot(x,y,col='darkgrey',xlab='x',ylab='y',main='Utilities Data Posterior Mean and 95% CI')
lines(x.new[idx],post.ci.lb[idx],col='red',lty=2)
lines(x.new[idx],post.ci.ub[idx],col='red',lty=2)

#Shade confidence bands.
polygon(c(sort(x.new),rev(sort(x.new))),c(post.ci.ub[idx],rev(post.ci.lb[idx])),col='lightgrey',border=NA)
points(x,y)
lines(x.new[idx],post.mean[idx],col='blue')

dev.off()

#================================================================
# E: Plot log of marginal likelihood for varying tau1.sq and b. =
#	Then use optimal parameters to compute posterior mean       =
#	for f given y.                                              =
#================================================================

#Grid mesh of parameters to test.
tau2.sq = 0
tau1.sq = seq(15,90,by=1)
b = seq(40,90,by=1)
triplets = expand.grid(b,tau1.sq,tau2.sq)

#Empty vector to hold marginal log-likelihoods.
ll = rep(0,nrow(triplets))

#Iterate through triplets.
for (k in 1:length(ll)){
	print(length(ll)-k)
	triplet = unlist(triplets[k,])
	ll[k] = gp.logl.y(x,y,mu=rep(0,length(x)),cov.fun=cov.se,params=triplet,sig2=1)
}

#Save optimal triplet of parameters.
max.idx = which.max(ll)
opt.triplet = unlist(triplets[max.idx,])
ll[max.idx]

#Contour plot of marginal loglikelihood.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/6_MargLogl_Y.pdf')
z = matrix(ll,byrow=F,nrow=length(b))
contour(b,tau1.sq,z,nlevels=20,xlab='b',ylab='tau1.sq',main=paste('Marginal Log-Likelihood of Y | sig2 = ',sig2),
	sub=paste('Optimal params: b = ',opt.triplet[1], " tau1.sq = ",opt.triplet[2]))
points(opt.triplet[1],opt.triplet[2],pch=19,cex=1,col='red')
dev.off()

#---------------------------------------------------------
#Use optimal parameters to calculate posterior mean of f|y,params,sig2.
# (Repeating Part C, but with optimal parameters.)

#Run prediction with sigma2=1 to estimate residuals.
pred = gp.predict(x,y,x.new=x,mu=rep(0,n),cov.fun=cov.se,params=opt.triplet,sig2=1)

#Vectors to hold posterior mean, var, and CI bounds.
post.mean = pred$post.mean
post.se = sqrt(pred$post.var)
post.ci.lb = post.mean - 1.96*post.se
post.ci.ub = post.mean + 1.96*post.se

#Plotting:
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/6_E_opt_params.pdf')
x.new=x
idx = sort(x.new, index.return = T)$ix
plot(x,y,col='darkgrey',xlab='x',ylab='y',main='Utilities Data Posterior Mean and 95% CI with Optimal Parameters',
	sub=paste('Optimal params: b = ',opt.triplet[1], " tau1.sq = ",opt.triplet[2]))
lines(x.new[idx],post.ci.lb[idx],col='red',lty=2)
lines(x.new[idx],post.ci.ub[idx],col='red',lty=2)

#Shade confidence bands.
polygon(c(sort(x.new),rev(sort(x.new))),c(post.ci.ub[idx],rev(post.ci.lb[idx])),col='lightgrey',border=NA)
points(x,y)
lines(x.new[idx],post.mean[idx],col='blue')

dev.off()
