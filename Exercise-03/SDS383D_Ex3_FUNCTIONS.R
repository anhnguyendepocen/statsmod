#SDS 383D - Exercise 3
#Functions
#Jennifer Starling
#Feb 2017

#================================================================
# Kernel Functions  =============================================
#================================================================

K_uniform = function(x){
	#-------------------------------------------------------------
	#FUNCTION: 	Uniform kernel.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector of values.
	#OUTPUT:	k = a scalar or vector of smoothed x values.
	#-------------------------------------------------------------
	k = .5 * ifelse(abs(x)<=1,rep(1,length(x)),rep(0,length(x)))
	return(k)
}

K_gaussian = function(x){
	#-------------------------------------------------------------
	#FUNCTION: Gaussian kernel.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector of values.
	#OUTPUT:	k = a scalar or vector of smoothed x values.
	#-------------------------------------------------------------
	k = (1/sqrt(2*pi)) * exp(-x^2/2)
	return(k)
}

#================================================================
# Noisy Data Simulation =========================================
#================================================================

sim_noisy_data = function(x,f,sig2){
	#-------------------------------------------------------------
	#FUNCTION: 	Simulates noisy data from some nonlinear function f(x).
	#-------------------------------------------------------------
	#INPUTS:	x = independent observations.
	#			f = function f(x) to simulate from
	#			sig2 = variance of the e ~ N(0,sig2) noise.
	#-------------------------------------------------------------
	#OUTPUTS:	x = generated x values.
	#			y = generated y = f(x) + e values.
	#-------------------------------------------------------------
	fx = f(x)
	e = rnorm(length(x),0,sqrt(sig2))
	return(y = fx+e)
}

#================================================================
# Smoothing Functions ===========================================
#================================================================

linear_smoother = function(x,y,x_star,h=1,K){
	#-------------------------------------------------------------
	#FUNCTION: 	Linear smoothing function for kernel regression.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector of regression covariates.
	#			x_star = scalar or vector of new x values for prediction.
	#			h = a positive bandwidth.
	#			K = a kernel function.  Default is set to Gaussian kernel.
	#OUTPUT:	yhat = a scalar or vector of smoothed x values.
	#-------------------------------------------------------------
	yhat=0	#Initialize yhat.
	
	for (i in 1:length(x_star)){
		w = (1/h) * K((x-x_star[i])/h) #Calculates weights.
		w = w / sum(w)					#Normalize weights.
		yhat[i] = crossprod(w,y)
	}
	return(yhat) #UPDATE TO INCLUDE WEIGHTS AS PART OF OUTPUT	
}	

local_linear_smoother = function(x,y,x_star,h=1,K){
	#-------------------------------------------------------------
	#FUNCTION: 	Local linear smoothing function for special case D=1 polynomial regression.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector of regression covariates.
	#			x_star = a scalar or vector of target point(s)
	#			h = a positive bandwidth.
	#			K = a kernel function.  Default is set to Gaussian kernel.
	#-------------------------------------------------------------
	#OUTPUT:	yhat = a scalar or vector of smoothed x values.
	
	n = length(x) 					#Number of obs.
	D = 1							#Degrees of local polynomial fctn.
	n_star = length(x_star)			#Number of target points at which to predict.
	yhat = rep(0,n_star)		#Vector to hold predicted function values at target point(s).
	
	#Error check for appropriate h value. 
	#try(if(h < 2/n || h > 1) stop("Enter h between (D+1)/n and 1"))
	 
	#Matrix to hold diags of weight matrices. (Col of weights for each x_star.)
	weights = matrix(0,n,n_star)	
	#NOTE: There is a separate H matrix for each x_star. Diagonals stored as columns.
	
	#Loop through each target point.
	for (b in 1:n_star){
		#Calculate unnormalized weights.
		s = rep(0,2)
	
		for (j in 1:2){
			s[j] = sum( K((x_star[b]-x)/h) * (x-x_star[b])^j )
		}
	
		w = K((x_star[b]-x)/h) * (s[2] - (x-x_star[b])*s[1])
	
		#Normalize weights.
		w = w / sum(w)
		
		#Save normalized weights.
		weights[,b] = w
	
		#Calculate function value for target point b.
		yhat[b] = crossprod(w,y)
		
		#CHECK: Print non-zero weight contributors.
		#print('Non-zero contributors to estimate:')
		#test_df = cbind(x=x,y=round(y,2),w=round(w,2),fhat=as.numeric(fhat_xstar[b]))
		#print(test_df[w>.001,])
	}
	
	#Return function output.
	return(list(yhat=yhat,weights=weights))
}

local_poly_smoother = function(x,y,x_star,h=1,K=K_gaussian,D=1){
	#-------------------------------------------------------------
	#FUNCTION: 	Smoothing function for local polynomial regression.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector of regression covariates.
	#			x_star = scalar new x value for prediction.
	#			h = a positive bandwidth.
	#			K = a kernel function.  Default is set to Gaussian kernel.
	#			D = degree of polynomial.  Default = 1.
	#-------------------------------------------------------------
	#OUTPUT:	yhat = a scalar or vector of smoothed x values.
	#-------------------------------------------------------------
	n = length(x)	#For i indices.
	yhat = rep(0,length(x_star))	#To store function outputs.
	weights.mat = matrix(0,nrow=n,ncol=length(x_star)) #To hold weights for each x-star.
	ahat = matrix(0,ncol=length(x_star),nrow=D+1)
	Hat = matrix(0,length(x_star),n)	#Store projection matrix.  Consists of row 1 of R times each xstar intercept element.
	
	for (b in 1:length(x_star)){	#Loop through x_star points.
		
		xstar.i = x_star[b] #Pick off current x_star.
		
		#Calculate (nxn) weights matrix W.
		W = diag( (1/h) * K((x-xstar.i)/h) )
		
		#Set up R matrix.  {R_ij} = (x_i-x)^j for j in 1...D+1.
		R = matrix(0,nrow=n,ncol=D+1)
	
		for (i in 1:n){
			for (j in 1:(D+1)){
				R[i,j] = (x[i]-xstar.i)^(j-1)
			}
		}
	
		#Precache t(R) %*% W.
		RtW = t(R) %*% W
		
		#Calculate ahat.
		ahat.xstar = solve(RtW %*% R) %*% RtW %*% y
		
		#Calculate hat matrix. First row of ahat.xstar, without y.
		Hat[b,] = (solve(RtW %*% R) %*% RtW)[1,]
		
		#Estimated function value.
		yhat[b] = ahat.xstar[1]
		
		#Save ahat parameters for each x_star.
		ahat[,b] = ahat.xstar	
	}
	return(list(yhat=yhat,weights=weights.mat,ahat=ahat,Hat = Hat))
}

#================================================================
# Cross-Validation / Bandwidth Tuning ===========================
#================================================================

tune_h = function(test,train,K,h){
	#-------------------------------------------------------------
	#FUNCTION: 	Function to tune bandwidth h for 
	#				specified test/train data sets and 
	#				specified kernel K for linear smoother.
	#-------------------------------------------------------------
	#INPUTS:	test = a test data set. Must have two cols, x and y.
	#			train = a training data set. Must have two cols, x and y.
	#			K = the kernel function
	#			h = a scalar or vector of bandwidths to try.
	#OUTPUTS:	pred_err_test = prediction error for testing data for each h.
	#			fhat_test = predicted values for testing data's x vals.
	#-------------------------------------------------------------
	#Extract training and test x and y vectors.
	x = train[,1]
	y = train[,2]
	x_star = test[,1]
	y_star = test[,2]
	
	#Calculate predicted points for test data set, and prediction error.
	yhat = linear_smoother(x,y,x_star,h,K)
	
	#Calculate predicted points for test data set, and prediction error.
	pred_err_test = sum(yhat-y_star)^2 / (length(x))
	
	#Return function outputs:
	return(list(yhat=yhat, pred_err_test = pred_err_test))
}

tune_h_loocv = function(x,y,K,h){
	#-------------------------------------------------------------
	#FUNCTION: 	Function to use leave on out LOOCV to tune  
	#				bandwidth h for specified test/train data sets and 
	#				specified kernel K for linear smoother.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector; the dependent variable.
	#			y = a scalar or vector; the independent var.  Same length as x.
	#			K = the kernel function
	#			h = a scalar or vector of bandwidths to try.
	#OUTPUTS:	pred_err_test = prediction error for testing data for each h.
	#			fhat_test = predicted values for testing data's x vals.
	#-------------------------------------------------------------
	
	#Define ppm H, the 'smoothing matrix'.  {H_ij} = 1/h * K((xi-xj*)/h)
	#For loocv, x=x*, since calculating on single data set instead of test/train.
	
	Hat = matrix(0,nrow=length(x),ncol=length(x)) #Empty matrix.
	
	for (i in 1:length(x)){			#Loop through H rows.
		for (j in 1:length(x)){		#Loop through H cols.
			Hat[i,j] = (1/h) * K((x[j] - x[i])/h)	
		} 	
		#Normalize weights by dividing H by rowsums(H).						
		Hat[i,] = Hat[i,] / sum(Hat[i,])
	}
	
	#Calculate predicted values.
	yhat = Hat %*% y
	
	#Calculate loocv prediction error.
	loocv_err = sum(((y-yhat)/(1-diag(Hat)))^2)
	
	#Return function outputs:	
	return(list(yhat=yhat,loocv_err=loocv_err))	
}

tune_h_local_poly_loocv = function(x,y,K,h){
	#-------------------------------------------------------------
	#FUNCTION: 	Function to use leave on out LOOCV to tune  
	#				bandwidth h for specified test/train data sets and 
	#				specified kernel K for local linear smoother.
	#-------------------------------------------------------------
	#INPUTS:	x = a scalar or vector; the dependent variable.
	#			y = a scalar or vector; the independent var.  Same length as x.
	#			K = the kernel function
	#			h = a scalar or vector of bandwidths to try.
	#-------------------------------------------------------------
	#OUTPUTS:	pred_err_test = prediction error for testing data for each h.
	#			fhat_test = predicted values for testing data's x vals.
	#-------------------------------------------------------------
	#Use existing x obs as target points for LOOCV.
	x_star = x
	
	#Call local_linear_smoother function to obtain yhat for each x point.
	output = local_poly_smoother(x,y,x_star,h=h,K_gaussian,D=1)
	yhat = output$yhat
	Hat = output$Hat
	
	Hii = diag(Hat)
	
	#Calculate loocv prediction error.
	loocv_err = sum(((y-yhat)/(1-Hii))^2)
	
	return(list(yhat=yhat,loocv_err=loocv_err))
}
