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
library(plotrix)		#For plotting CIs.
library(truncnorm)		#For sampling truncated normal.
library(matrixStats)	#For colSds.

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
X = model.matrix(bush ~ edu + age + female + black,data)
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

#Hierarchical augmented probit model; allows intercept to change between states.
#	Includes ad*price interaction.
hlm = glmer(bush ~ edu + age + female + black  + (1|state), data=data, family='binomial')
summary(hlm)

CI = confint(hlm)[-c(1:2),]
CI = cbind(CI,'50%'=rowMeans(CI))

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
ols = lm(bush ~ edu + age + female + black,data)
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

#----------------------------------------------------------------
### PLOT 1: Bayesian Credible Intervals by State.

#Calculate means and sds of posterior probabilities of supporting Bush, for each state.
py.state.mean = unlist(lapply(output$py.i,mean))
py.state.sd = unlist(lapply(output$py.i,sd))

#Set up BCI bounds.
lb = py.state.mean - 1.96 * py.state.sd
ub = py.state.mean + 1.96 * py.state.sd
CI = as.matrix(cbind(lb,m=py.state.mean,ub))

#Set up color vector for posterior mean above/below 50%.
state.col = ifelse(py.state.mean>=.5,'red','blue')

#Save plot to pdf file.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/State_BCI.pdf',width=20,height=10)
plotCI(x=CI[,2],li=CI[,1],ui=CI[,3],xlab='Coefficients',ylab='',col=state.col,
	,main='Bayesian Credible Intervals by State for P(Bush)',xaxt='n')
axis(1,at=seq(1:s),labels=states.list,cex.axis=.65)
abline(h=.5,lty=2,col='darkgrey')
dev.off()

#----------------------------------------------------------------
### PLOT 2: Bayesian Credible Intervals by Education.  Does education affect probability?

#Calculate means and sds of posterior probabilities of supporting Bush, for each edu level above HS.
py.all = unlist(output$py.i.pm)
py.edu.mean = colMeans(py.all*X[,2:4])
py.edu.sd = colSds(py.all*X[,2:4])

#Set up BCI bounds.
lb = py.edu.mean - 1.96 * py.edu.sd
ub = py.edu.mean + 1.96 * py.edu.sd
CI = as.matrix(cbind(lb,m=py.edu.mean,ub))

#Set up color vector for posterior mean above/below 50%.
edu.col = ifelse(py.edu.mean>=.5,'red','blue')
edu.list = colnames(X[,2:4])

#Save plot to pdf file.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/edu_BCI.pdf',width=10,height=10)
plotCI(x=CI[,2],li=CI[,1],ui=CI[,3],xlab='Coefficients',ylab='',col=edu.col,
	,main='BC Intervals by Ed Level Above NoHS for P(Bush)',xaxt='n')
axis(1,at=seq(1:3),labels=edu.list,cex.axis=1)
abline(h=.5,lty=2,col='darkgrey')
dev.off()

#----------------------------------------------------------------
### PLOT 3: Does age affect probability?

#Calculate posterior means and sds for P(Bush) for age levels above 18-29.
py.all = unlist(output$py.i.pm)
py.age.mean = colMeans(py.all*X[,5:7])
py.age.sd = colSds(py.all*X[,5:7])

#Set up BCI bounds.
lb = py.age.mean - 1.96 * py.age.sd
ub = py.age.mean + 1.96 * py.age.sd
CI = as.matrix(cbind(lb,m=py.age.mean,ub))

#Set up color vector for posterior mean above/below 50%.
age.col = ifelse(py.age.mean>=.5,'red','blue')
age.list = colnames(X[,5:7])

#Save plot to pdf file.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/age_BCI.pdf',width=10,height=10)
plotCI(x=CI[,2],li=CI[,1],ui=CI[,3],xlab='Coefficients',ylab='',col=age.col,
	,main='BC Intervals by Age Level Above 18-29 for P(Bush)',xaxt='n')
axis(1,at=seq(1:3),labels=age.list,cex.axis=1)
abline(h=.5,lty=2,col='darkgrey')
dev.off()

#----------------------------------------------------------------
### PLOT 4: Do race and gender affect probability?

#Calculate posterior means and sds for each race/gender combo.
py.all = unlist(output$py.i.pm)
py.rg.mean = c(
	mean(py.all[data$female==0 & data$black==0]),	#White Male
	mean(py.all[data$female==1 & data$black==0]),	#White Female
	mean(py.all[data$female==0 & data$black==1]),	#Black Male
	mean(py.all[data$female==1 & data$black==1]))	#Black Female

py.rg.sd = c(
	sd(py.all[data$female==0 & data$black==0]),		#White Male
	sd(py.all[data$female==1 & data$black==0]),		#White Female
	sd(py.all[data$female==0 & data$black==1]),		#Black Male
	sd(py.all[data$female==1 & data$black==1]))		#Black Female

#Set up BCI bounds.
lb = py.rg.mean - 1.96 * py.rg.sd
ub = py.rg.mean + 1.96 * py.rg.sd
CI = as.matrix(cbind(lb,m=py.rg.mean,ub))

#Set up color vector for posterior mean above/below 50%.
rg.col = ifelse(py.rg.mean>=.5,'red','blue')
rg.list = c('Wh-M','Wh-F','B-M','B-F')

#Save plot to pdf file.
pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Probit/race_gender_BCI.pdf',width=10,height=10)
plotCI(x=CI[,2],li=CI[,1],ui=CI[,3],xlab='Coefficients',ylab='',col=rg.col,
	,main='BC Intervals by Race and Gender for P(Bush)',xaxt='n')
axis(1,at=seq(1:4),labels=rg.list,cex.axis=1)
abline(h=.5,lty=2,col='darkgrey')
dev.off()

