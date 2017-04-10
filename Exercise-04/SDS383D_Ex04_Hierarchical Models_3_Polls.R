#Stats Modeling 2
#Exercise 4
#Polls - Augmented Probit Hierarchical Model (Albert & Chib 1993)

#================================================================
# Setup =========================================================
#================================================================
rm(list=ls())

library(mosaic)
library(lme4)
library(dplyr)
library(plotrix)	#For plotting CIs.

#Set working directory
setwd('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/')

#Read in data.
data = read.csv('Data/polls.csv')
attach(data)

#Read in functions.
source('RCode/SDS383D_Ex04_Hierarchical Models_FUNCTIONS_3_Polls (Augmented Probit).R')

#Remove missing observations (178 rows out of 2193).
#Note: Could use Gibbs to sample and fill in missing y obs, but did not here.
data = data[complete.cases(data),]

#Order data by state.
data = data[with(data,order(data$state)), ]
rownames(data) = 1:nrow(data)

#List and number of unique states.
states.list = sort(unique(data$state))
s = length(states.list)

#Set edu=NoHS and age=18to29 as the reference levels.
data$edu = factor(data$edu, levels = c('NoHS','HS','SomeColl','Bacc'))
data$age = factor(data$age, levels = c('18to29','30to44','45to64','65plus'))

#Set up responses vector and design matrix.
y = data$bush
X = model.matrix(bush ~ edu + age + female + black + weight,data)
idx = cbind.data.frame(state.name = data$state,state=as.numeric(data$state))

#Set up marginal effects coding for ordinal variables.
#HS : NoHS = 1, HS = 1, SomeColl = 0, Bacc = 0
#30-44: 19to29 = 1, 30to44=1, 45to60=0, 60plus=0
#Then SomeColl is marginal effect of some college beyond HS.

for (i in 1:nrow(X)){
	#Adjust edu.
	if (X[i,4] == 1) X[i,3] = X[i,2] = 1
	if (X[i,3] == 1) X[i,2] = 1
	
	#Adjust age.
	if(X[i,7]==1) X[i,6] = X[i,5] = 1
	if(X[i,6]==1) X[i,5] = 1
}

#================================================================
# Hierarchical Model: Empirical Bayes (Using lmer) ==============
#================================================================

# This model uses empirical Bayes; errors estimated from data.
data$wt.sc = scale(data$weight)

#Hierarchical augmented probit model; allows intercept to change between states.
#	Includes ad*price interaction.
hlm = glmer(bush ~ edu + age + female + black  + wt.sc + (1|state), data=data, family='binomial')
summary(hlm)

CI = confint(hlm1)[-c(1:2),]
CI = cbind(CI,'50%'=rowMeans(CI2))

#Plot confidence intervals for fixed coefficients.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/lmer_CI_plot.pdf',height=8,width=16)
plotCI(x=CI[,3],li=CI[,1],ui=CI[,2],xlab='Coefficients',ylab='',main='Coefficient Estimates from LMER',xaxt='n')
axis(1,at=1:length(colnames(X)[-1]),labels=colnames(X)[-1],cex.axis=.75)
abline(h=0,lty=2,col='red')
dev.off()

#Variance for state = 0.1732

#================================================================
# Run Gibbs Sampler =============================================
#================================================================

#Set up hyperparameters.
p = ncol(X)

#Since Beta fixed, no prior; estimating based on OLS coefficients.
ols = lm(bush ~ edu + age + female + black + weight,data)
muB = ols$coef
Sigma = vcov(ols)

#Runs sampler.
output = gibbs.probit(y,X,idx,muB,Sigma,iter=11000,burn=1000,thin=2)

#================================================================
# Trace Plots for MCMC Mixing ===================================
#================================================================
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/mcmc_trace.pdf',width=20,height=10)
par(mfrow=c(3,3))
plot(output$tau.sq,type='l',main='tau.sq')
plot(output$Beta[1,],type='l',main='Beta')
plot(output$Beta[2,],type='l',main='Beta')
plot(output$Beta[3,],type='l',main='Beta')
plot(output$Beta[4,],type='l',main='Beta')
plot(output$mu.i[3,],type='l',main='mu.3')
plot(output$mu.i[10,],type='l',main='mu.10')
plot(output$mu.i[30,],type='l',main='mu.30')
plot(output$mu.i[40,],type='l',main='mu.40')
dev.off()

#Traces of mu.i + Beta[1]
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/mcmc_trace_intercepts.pdf',width=20,height=10)
par(mfrow=c(3,1))
plot(output$Beta[1,] + output$mu.i[1,],type='l',main='Intercept for State 1')
plot(output$Beta[1,] + output$mu.i[2,],type='l',main='Intercept for State 2')
plot(output$Beta[1,] + output$mu.i[3,],type='l',main='Intercept for State 3')
dev.off()

#================================================================
# Check Prediction Error ========================================
#================================================================

#(Since classifier, look at prediction error.)

#Generate Y values using y_{ij} ~ Bernoulli(p_{ij}) where p_{ij} are posterior means.
#Since sorted data by state, can unlist.
y.pred = output$y.pred

#Show table of y actual versus predicted.
table(y,y.pred)

#Error rate.
sum(y!=y.pred)/length(y)

#================================================================
# Plot Demographic Results ======================================
#================================================================

#Goal: How demographics relate to probability of supporting Bush, ie P(Y_ij = 1).

#Bayesian Credible Intervals.
py.state.mean = unlist(lapply(output$py.i,mean))
py.state.sd = unlist(lapply(output$py.i,sd))

lb = py.state.mean - 1.96 * py.state.sd
ub = py.state.mean + 1.96 * py.state.sd
CI = as.matrix(cbind(lb,m=py.state.mean,ub))

pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/State_BCIntervals.pdf',width=20,height=10)
plotCI(x=CI[,2],li=CI[,1],ui=CI[,3],xlab='Coefficients',ylab='',main='Bayesian Credible Intervals by State for P(Supporting Bush)',xaxt='n')
axis(1,at=seq(1:s),labels=states.list,cex.axis=.5)
dev.off()

#PLOT 2: Does education affect probability?
py.all = unlist(py)

pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/edu_boxplot.pdf',width=20,height=10)
boxplot(unlist(py.all)~data$edu,main='Probability of Supporting Bush by Education Level',ylab='P(Supporting Bush)')
dev.off()

#PLOT 3: Does age affect probability?
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/age_boxplot.pdf',width=20,height=10)
boxplot(unlist(py.all)~data$age,main='Probability of Supporting Bush by Age Group',ylab='P(Supporting Bush)')
dev.off()

#PLOT 4: Do race and gender affect probability?
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/race_gender_boxplot.pdf',width=20,height=10)
boxplot(unlist(py.all)~data$female + data$black,main='Probability of Supporting Bush by Race and Gender',
	ylab='P(Supporting Bush)',xaxt='n')
axis(1,at=1:4,labels=c("M-White","F-White","M-Black","F-Black"))
dev.off()

#PLOT 5: Does weight affect probability?
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/weight.pdf',width=20,height=10)
wt.idx = order(data$weight)
plot(data$weight[wt.idx],py.all[wt.idx],pch=1,col='blue',xlab='Weight',ylab='P(Supporting Bush)',	
	main='Probability of Supporting Bush by Weight')
dev.off()

