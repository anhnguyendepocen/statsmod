#Stats Modeling 2
#Exercise 3
#Local Polynomial Regression - Part E-G
#Jennifer Starling
#Feb 21, 2017

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
x = utilities$temp								#average temp.
y = utilities$gasbill / utilities$billingdays	#avg daily bill

#================================================================
# Analysis of Residuals: ========================================
#================================================================

#Fit D=1 local polynomial model (local linear) using optimal h.
#Inspect residuals.  Does homoscedasticity look reasonable?  If not, propose fix.
h=6.9
yhat = local_linear_smoother(x,y,x,h=h,K_gaussian)$yhat
yhat_log = local_linear_smoother(x,log(y),x,h=h,K_gaussian)$yhat
resids = (y-yhat)
resids_log = (log(y)-yhat_log)

### Plotting
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/Local_Linear_Residuals.pdf',height=6,width=12)
par(mfrow=c(1,2))
plot(x,resids,main=paste('Y: Local Linear Smoothing Residuals, Gaussian Kernel, h = ',sep="",h))
plot(x,resids_log,main=paste('LOG(Y): Local Linear Smoothing Residuals, Gaussian Kernel, h = ',sep="",h))	
dev.off()

#================================================================
# Local Linear Smoother: Optimize bandwidth h using loocv =======
#================================================================

#Candidate h values. 
H = seq(1,8,by=.1)

#Vectors to hold prediction errors and estimated function values for each h.
pred_err = rep(0,length(H))

#Loop through h values.
for (m in 1:length(H)){
	h = H[m]
	results = tune_h_local_poly_loocv(x,y,K=K_gaussian,h=h)
	pred_err[m] = results$loocv_err
}

#Select optimal h and obtain fitted values.
h_opt = H[which.min(pred_err)]
yhat = local_poly_smoother(x,y,x,h=h_opt,K_gaussian,D=1)$yhat


#------------------------------------------------------------
### Plot output with fitted data using optimal h value.

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/Local_Linear_Opt_h.pdf')

#Scatterplot of x/y.
plot(x,y,main=paste('Local Linear Smoothing, Gaussian Kernel, h = ',sep="",h_opt))	

#Overlay estimated fit.
x_star=x
idx = sort(x_star, index.return = T)$ix
	lines(sort(x_star),yhat[idx],col='firebrick3',lwd=2) #Fitted line.

dev.off()


#================================================================
# Confidence Bands: =============================================
#================================================================
#sigma2 for CI band calculation.  Need Var(f*(x)).
RSS = sum(resids^2)
sigma2_hat = RSS / (length(yhat)-1)

#For confidence bands calculation. Other piece of Var(f*(x)).
Hat = local_poly_smoother(x,y,x,h=h_opt,K_gaussian,D=1)$Hat 

#Var(f*(x))
var = rowSums(Hat^2) * sigma2_hat

lb = yhat - 1.96*sqrt(var)
ub = yhat + 1.96*sqrt(var)

#Plot bands.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/LLS_Conf_Bands.pdf')
#Scatterplot of x/y.
plot(x,y,main=paste('Local Lin Smoothing 95% Confidence Bands, h = ',sep="",h_opt))	
	
#Overlay confidence bands.
lines(sort(x_star),lb[idx],col='blue',lwd=2,lty=2) 	
lines(sort(x_star),ub[idx],col='blue',lwd=2,lty=2) 

#Shade confidence bands.
polygon(c(sort(x_star),rev(sort(x_star))),c(ub[idx],rev(lb[idx])),col='thistle',border=NA)
points(x,y,main=paste('Local Lin Smoothing 95% Confidence Bands, h = ',sep="",h_opt))	

#Overlay estimated fit.
x_star=x
idx = sort(x_star, index.return = T)$ix
	lines(sort(x_star),yhat[idx],col='firebrick3',lwd=2) #Fitted line.
	
dev.off()	

#================================================================
# Local Polynomial Smoother Degrees Demo: =======================
#================================================================

x_star = x
h=6.9

#Local linear smoother - James solution.
fx_hat = local_poly_smoother(x,y,x_star,h,K=K_gaussian,D=1)$yhat

#Local linear smoother - via polynomial smoother with D=1.
fx_hat2 = local_poly_smoother(x,y,x_star,h,K=K_gaussian,D=2)$yhat

#Compare to linear smoother from previous ex part.	
fx_hat3 = local_poly_smoother(x,y,x_star,h,K=K_gaussian,D=3)$yhat	

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/Polynomial_Degree_Comparison.pdf')
#Scatterplot of x/y.
plot(x,y,main=paste('Local Polynomial Regression, Gaussian Kernel, h = ',sep="",h))	

#Overlay estimated fits.
idx = sort(x_star, index.return = T)$ix
	lines(sort(x_star),fx_hat[idx],col='firebrick3',lwd=2) #Fitted line.

idx = sort(x_star, index.return = T)$ix
	lines(sort(x_star),fx_hat2[idx],col='forestgreen',lwd=2) #Fitted line.

idx = sort(x_star, index.return = T)$ix
	lines(sort(x_star),fx_hat3[idx],col='blue',lwd=2) #Fitted line.	
	
legend('topright',col=c('firebrick3','forestgreen','blue'),legend=c('D=1','D=2','D=3'),lty=1)	
dev.off()
