#PEER REVIEW FOR SPENCER
#Some useful C++ functions.

#=====================================================================
#=== Data Generation =================================================
#=====================================================================

#Housekeeping
rm(list=ls())

#Load in libraries.
library(Rcpp)
library(RcppArmadillo)

#Load CPP file.
sourceCpp('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Peer Reviews/Useful_Rcpp_Functions.cpp')

#Call Rcpp functions.
rcpp_seq(1,10,2)
mvrnorm(2,c(0,0),diag(2))
rinvwish(1,4,diag(2))

