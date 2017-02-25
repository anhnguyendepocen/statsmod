#Stats Modeling 2
#Exercise 3
#Cross-Validation - Part B
#Jennifer Starling
#Feb 15, 2017

#------------------------------------------------------------
### Imagine a 2x2 tablefor the unknown, true
###	state of affairs.  Rows are 'wiggly' and 'smooth' functions,
###	and cols are 'highly noisy obs' and 'less noisy obs'.
### Simulate one data set (n=500) for each of four cells of table.
### Split each data set into test and train sets.  Apply method to
### each case.  Apply function from part A to select bandwidth.
#------------------------------------------------------------

#================================================================
# Environment Setup. ============================================
#================================================================

#Housekeeping.
rm(list=ls())

#Load latex table library.
library(xtable)	#For table output to latex.
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

#Load functions.
source('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/RCode/SDS383D_Ex3_FUNCTIONS.R')

#================================================================
# Data and function generation. =================================
#================================================================

#Set up wiggly and less-wiggly functions.
fwiggly = function(x) sin(10*pi*x)
fsmooth = function(x) sin(2*pi*x)

#Set up x values on the unit interval.
n = 500
x = runif(n,0,1)

#Generate noisy and less-noisy data for
#each fwiggly and fsmooth function.
y_wiggly_noisy	= sim_noisy_data(x,fwiggly,sig2=.75)
y_wiggly_clear	= sim_noisy_data(x,fwiggly,sig2=.1)
y_smooth_noisy	= sim_noisy_data(x,fsmooth,sig2=.75)
y_smooth_clear	= sim_noisy_data(x,fsmooth,sig2=.1)

#Combine into single matrix for convenience.
y = matrix(c(y_wiggly_noisy,y_wiggly_clear,y_smooth_noisy,y_smooth_clear),ncol=4,byrow=F)
colnames(y) = c('wiggly_noisy','wiggly_clear','smooth_noisy','smooth_clear')

#Scatter plots to preview the four data sets.
#pdf(file='/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/cv_function_grid.pdf')
	par(mfrow=c(2,2))
	for(j in 1:ncol(y)){
		plot(x,y[,j],main=paste(colnames(y)[j]),xlab='x',ylab='y',ylim=c(-5,5))
	}
#dev.off()

#================================================================
# Cross-Validation to Tune h (One Split) ========================
# (Part B)                               ========================
#================================================================

#------------------------------------------------------------
### Split each data set into train and test sets. (70-30)

test_idx = sample(1:n,size=n*.3,replace=F)
test = cbind(x[test_idx],y[test_idx,])
train = cbind(x[-test_idx],y[-test_idx,])

#------------------------------------------------------------
### Optimize bandwidth h.

# For each case, select a bandwidth parameter. (Used Gaussian kernel.)
h_opt = rep(0,4)	#Vector to hold optimal h values.
names(h_opt) = colnames(test)[-1]
H = seq(.001,1,by=.001)	#Candidate h values.
yhat = list()			#To hold estimated function values.
x_star = test[,1]

#Iterate through each setup.
for (i in 1:4){
	#Extract x and just the y column required.
	tr = train[,c(1,i+1)]
	te = test[,c(1,i+1)]
	
	#Temp vector to hold prediction errors.
	temp_pred_err = rep(0,length(H))
	
	for (j in 1:length(H)){ 
		h = H[j]
		results = tune_h(test=te,train=tr,K=K_gaussian,h=h)	
		temp_pred_err[j] = results$pred_err_test
	}

	h_opt[i] = H[which.min(temp_pred_err)]
	yhat[[i]] = tune_h(te,tr,K=K_gaussian,h_opt[i])$fhat_star
}

#------------------------------------------------------------
### Output results as table and plot.

#Format optimal h results as matrix.
h_opt_mat = matrix(h_opt,nrow=2,byrow=T)
colnames(h_opt_mat) = c('noisy','clear')
rownames(h_opt_mat) = c('wiggly','smooth')
xtable(h_opt_mat,digits=3)	#Output latex table.

#Plot output with fitted data using optimal h values.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/Bandwidth_selection.pdf')
par(mfrow=c(2,2))
for (i in 1:4){
	#Scatterplot of test x/y.
	plot(test[,1],test[,i+1],
		main=paste(names(h_opt)[i],", h=",h_opt[i],sep=''),
		xlab='x_train',ylab='y_train')	
	#Overlay estimated fit.
	idx = sort(x_star, index.return = T)$ix
	lines(sort(x_star),yhat[[i]][idx],col='red') #Fitted line.
}
dev.off()

#================================================================
# LOOCV to Tune h  ==============================================
# (Part C)         ==============================================
#================================================================

# For each case, select a bandwidth parameter. (Used Gaussian kernel.)
h_opt_loocv = rep(0,4)				#Vector to hold optimal h values.
names(h_opt_loocv) = colnames(y)
H = seq(.001,1,by=.001)				#Candidate h values.
yhat_loocv = list()				#To hold estimated function values.

#Iterate through each setup.
for (i in 1:4){
	
	#Extract y column.
	ytemp = y[,i]
	
	#Temp vector to hold prediction errors.
	temp_pred_err = rep(0,length(H))
	
	for (j in 1:length(H)){ 
		h = H[j]
		results = tune_h_loocv(x=x,y=ytemp,K=K_gaussian,h=h)	
		temp_pred_err[j] = results$loocv_err
	}

	h_opt_loocv[i] = H[which.min(temp_pred_err)]
	yhat_loocv[[i]] = tune_h_loocv(x,ytemp,K=K_gaussian,h_opt_loocv[i])$yhat
}

#------------------------------------------------------------
### Output results as table and plot.

#Format optimal h results as matrix.
h_opt_mat_loocv = matrix(h_opt_loocv,nrow=2,byrow=T)
colnames(h_opt_mat_loocv) = c('noisy','clear')
rownames(h_opt_mat_loocv) = c('wiggly','smooth')
xtable(h_opt_mat_loocv,digits=3)	#Output latex table.

#Plot output with fitted data using optimal h values.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/Bandwidth_selection_loocv.pdf')
par(mfrow=c(2,2))
for (i in 1:4){
	#Scatterplot of test x/y.
	plot(x,y[,i],
		main=paste(names(h_opt_loocv)[i],", h=",h_opt_loocv[i],sep=''),
		xlab='x',ylab='y')	
	#Overlay estimated fit.
	idx = sort(x, index.return = T)$ix
	lines(sort(x),yhat_loocv[[i]][idx],col='red') #Fitted line.
}
dev.off()
