#Stats Modeling 2
#Exercise 3
#GPs in Nonparametric regression & spatial smoothing


#================================================================
# Environment Setup & Data Load =================================
#================================================================

#Housekeeping.
rm(list=ls())

#Plotting libraries
library(ggmap)

#Load functions.
source('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/RCode/SDS383D_Ex3_FUNCTIONS_Gaussian_Processes_Rcpp.R')

#Load data.
weather = read.csv('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/statsmod/Course-Data/weather.csv',header=T)

n = nrow(weather)
x = as.matrix(weather[,3:4])

#================================================================
# Fit a GP Model for each of the temp and pressure variables: ===
# Choose Optimal Parameters tau1.sq and b.                    ===
#================================================================

y = weather$temperature

#Run prediction with sigma2=1 to estimate residuals.
pred = gp.predict(x,y,x.new=x,mu=rep(0,nrow(x)),cov_se_2d, params = c(5,1,1,1E-6),sig2=1)$post.mean
sig2 = sum(y-pred)^2/(n-1)

#Choose optimal parameters for model using marginal loglikelihood.
tau2.sq = 1E-6
tau1.sq = .1
b1 = 10
b2= 10
params = c(b1,b2,tau1.sq,tau2.sq)

test.gp = gaussian_process(x,mu=rep(0,nrow(x)),cov.fun=cov_se_2d,params=params)
test.logl =  gp.logl.y(x,y,mu=rep(0,nrow(x)),cov.fun=cov_se_2d,params=params,sig2=1)


#----------------------------------------------------
#Weather: Grid mesh of parameters to test.
tau2.sq = 0
tau1.sq = seq(.1,15,length=20) 
b = b1 = b2 = seq(.1,2,length=20) 
triplets = expand.grid(b1,b2,tau1.sq,tau2.sq)

#Empty vector to hold marginal log-likelihoods.
ll = rep(0,nrow(triplets))

#Iterate through triplets.
for (k in 1:length(ll)){
	print(length(ll)-k)
	triplet = unlist(triplets[k,])
	ll[k] = gp.logl.y(x,y,mu=rep(0,nrow(x)),cov.fun=cov_se_2d,params=triplet,sig2)
}

#Save optimal triplet of parameters.
max.idx = which.max(ll)
opt.params.temp = unlist(triplets[max.idx,])

#----------------------------------------------------
#Pressure: Grid mesh of parameters to test.
y = weather$pressure
tau2.sq = 0
tau1.sq = seq(40000,60000,length=50)
b = b1 = b2 = seq(.01,10,length=50)
triplets = expand.grid(b1,b2,tau1.sq,tau2.sq)

#Empty vector to hold marginal log-likelihoods.
ll = rep(0,nrow(triplets))

#Iterate through triplets.
for (k in 1:length(ll)){
	print(length(ll)-k)
	triplet = unlist(triplets[k,])
	ll[k] = gp.logl.y(x,y,mu=rep(0,nrow(x)),cov.fun=cov_se_2d,params=triplet,sig2)
}

#Save optimal triplet of parameters.
max.idx = which.max(ll)
opt.params.pressure = unlist(triplets[max.idx,])

#================================================================
# Fit a GP Model for each of the temp and pressure variables: ===
# TEMPERATURE & PRESSURE FITTING                              ===
#================================================================

#Set up grid of x values to predict.
lon.range = range(x[,1])
lat.range = range(x[,2])

lon.new = seq(lon.range[1],lon.range[2],length=40)
lat.new = seq(lat.range[1],lat.range[2],length=40)

x.new = as.matrix(expand.grid(lon.new,lat.new))

#Predict temp values.
y = weather$temperature
temp.pred = gp.predict(x,y,x.new,mu=rep(0,nrow(x)),cov_se_2d, params = opt.params.temp,sig2=sig2)$post.mean

#Predict pressure values.
y = weather$pressure
pressure.pred = gp.predict(x,y,x.new,mu=rep(0,nrow(x)),cov_se_2d, params = opt.params.pressure,sig2=sig2)$post.mean

#================================================================
# Fit a GP Model for each of the temp and pressure variables: ===
# TEMPERATURE PLOTTING                                        ===
#================================================================

pred.values= cbind(x.new,temp.pred)
colnames(pred.values) = c('lon','lat','pred')

z = matrix(temp.pred,byrow=F,nrow=60)

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/7_Smoothed_Temp.pdf')

filled.contour(lon.new,lat.new,z,
    key.title = title(main="Temp"),
    color.palette = rainbow,
    plot.axes = {
        axis(1)
        axis(2)
        contour(lon.new,lat.new,z,add=T,lwd=2)             
    }
)
mtext(paste('Optimal Parameters: b1 = ',round(opt.params.temp[1],2),
	', b2 = ',round(opt.params.temp[2],2),
	', tau1.sq = ', round(opt.params.temp[3],2),
	', tau2.sq = 0',sep=""),side=3)
dev.off()
	
#================================================================
# Fit a GP Model for each of the temp and pressure variables: ===
# PRESSURE PLOTTING                                           ===
#================================================================

pred.values= cbind(x.new,pressure.pred)
colnames(pred.values) = c('lon','lat','pred')

z = matrix(pressure.pred,byrow=F,nrow=60)

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-03/Figures/7_Smoothed_Pressure.pdf')

filled.contour(lon.new,lat.new,z,
    key.title = title(main="Pressure"),
    color.palette = rainbow,
    plot.axes = {
        axis(1)
        axis(2)
        contour(lon.new,lat.new,z,add=T,lwd=2)             
    }
)
mtext(paste('Optimal Parameters: b1 = ',round(opt.params.pressure[1],2),
	' b2 = ',round(opt.params.pressure[2],2),
	' tau1.sq = ', round(opt.params.pressure[3],2),
	'tau2.sq = 0'),side=3)	
dev.off()	