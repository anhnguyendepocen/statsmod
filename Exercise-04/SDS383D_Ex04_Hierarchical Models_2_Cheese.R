#Stats Modeling 2
#Exercise 4
#Cheese - Price Elasticity Modeling (Hierarchical Models)

#================================================================
# Setup =========================================================
#================================================================
rm(list=ls())

library(mosaic)
library(lme4)
library(dplyr)

#Set working directory
setwd('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/')

#Read in data.
data = read.csv('Data/cheese.csv')
attach(data)

#Read in functions.
source('RCode/SDS383D_Ex04_Hierarchical Models_FUNCTIONS.R')

#Add log price and log volume variables.
data$logP = log(data$price)
data$logQ = log(data$vol)

#Add a week variable.
data$week = ave(rep(1:nrow(data)), data$store, FUN=seq_along)

#Add store numbers.
data$storenum =  ave(rep(1:nrow(data)), data$week, FUN=seq_along)

#================================================================
# Exploratory ===================================================
#================================================================

### Inspect grouping structure.

xtabs( ~ factor(week) + factor(store), data=data)
	#Shows that each week has one observation per store.
	#Notice there are 68 weeks, and some stores are missing a few obs.
xtabs( ~ factor(week), data=data)
	#Shows how many observations each week has. 
	#Can see there are some weeks near the end with 87, and
	#a few weeks with only 28.  Not a balanced model.	
	
### Does avg price appear to vary by store?  Choose a sample, since 88 hard to visualize.

vol.avg.by.store = aggregate(vol,list(store),mean)
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/cheese_exp_vol.pdf')
#par(mfrow=c(1,2))
#boxplot(vol ~ store)
plot(vol.avg.by.store$x,pch=19,main='Avg Volume by Store Over All Weeks',xlab='Store',ylab='Vol')
dev.off()

#================================================================
# Naive Approach ================================================
#================================================================

# This approach ignores the pseudoreplication due to weekly 
# measurements of same stores.

# Main effects
lm0 = lm(logQ ~ logP + disp, data = data)
lm1 = lm(logQ ~ logP + disp + store, data=data)
lm2 = lm(logQ ~ logP + store, data=data)

summary(lm0)
summary(lm1)
anova(lm1)

anova(lm1,lm0)	#F-test, concludes store is significant.
anova(lm2,lm0)	#F-test, concludes disp is significant.

#================================================================
# Hierarchical Model: Empirical Bayes ===========================
#================================================================

#----------------------------------------------------------------
#	Model: logQ_{it} = log(alpha_{i}) + b_{i}*log(P_{it}) + gamma_{i} * ad_{it} + v_{i} + e_{it}
#		for i = 1...88 indexing stores, and t indexing time points.
#		logQ_{it}	= Volume at each store at each time.
#
#	This model uses empirical Bayes; errors estimated from data.
#----------------------------------------------------------------

#Hierarchical linear model; allows intercept, logQ and disp to change among stores.
#	Includes ad*price interaction.
hlm = lmer(logQ ~ logP + disp + disp:logP + (1 + logP + disp + disp*logP | store), data=data)
summary(hlm)

#Plot conditional modes (same as means, for linear mixed models) of random effects.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/cheese_hlm_1.pdf')
resid = ranef(hlm, condVar = T)
dotplot(resid,scales=list(cex=c(.5,.5)),layout=c(3,1),main=T, main.title='Random Effects by by Store',)
dev.off()

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/cheese_resid.pdf')
plot(hlm)
dev.off()

#================================================================
# Hierarchical Model: Fully Bayesian w Gibbs Sampler ============
#================================================================

#Data setup for function.
y = data$logQ														#Responses for all i,t
X = as.data.frame(model.matrix(logQ ~ 1 + price + disp,data=data))	#Covariates for all i,t
idx = cbind.data.frame(store=data$store,							#Store/time index info.
	storenum=data$storenum,
	week=data$week)	

#Noninformative hyperpriors for m,v,C,d.
p = 3
m = rep(0,3)
v = 1
C = diag(p)
d = p+1

#Call gibbs sampler function.
output = gibbs.cheese(y,X,idx,m,v,C,d,iter=50,burn=10,thin=2)

#Display results.

output$Bi.post.mean
output$Sigma.post.mean
output$mu.post.mean
output$sig.sq.post.mean

par(mfrow=c(1,3))
plot(output$Bi.post.mean[,1])
plot(output$Bi.post.mean[,2])
plot(output$Bi.post.mean[,3])
