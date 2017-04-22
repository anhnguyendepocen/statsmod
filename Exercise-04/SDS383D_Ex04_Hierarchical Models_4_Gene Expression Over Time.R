#Stats Modeling 2
#Exercise 4
#Polls - Gene Expression Over Time

#================================================================
# Setup =========================================================
#================================================================
rm(list=ls())

#library(mosaic)
#library(lme4)
#library(dplyr)
#library(plotrix)	#For plotting CIs.
library(lattice)
library(mvtnorm)	#For sampling multivariate normals, and for ll function.

#Set working directory
setwd('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/')

#Read in data.
data = read.csv('Data/droslong.csv')

#Read in functions.
source('RCode/SDS383D_Ex04_Hierarchical Models_FUNCTIONS_4_Genes (Hierarchical GP Regression).R')


#================================================================
# Exploratory ===================================================
#================================================================

#Summarize everything.
xtabs( ~ factor(group) + factor(gene) + factor(replicate), data=data)

#Number of observations in each group.
xtabs( ~ factor(group), data=data)

#Number of each obs of each gene in each group.
xtabs (~ factor(label) + factor(group), data=data)

#How each gene varies over time.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Gene_Expression/Exploratory_GenesOverTime.pdf')
xyplot(log2exp ~ time | gene, data = data,main='Genes Over Time')
dev.off()

#How each group varies over time.
pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Gene_Expression/Exploratory_GroupsOverTime.pdf')
xyplot(log2exp ~ time | group, data = data, main='Groups Over Time')
dev.off()

#================================================================
# Sort and Organize Data ========================================
#================================================================

#Sort data.
data = data[with(data, order(group, gene, replicate, time)),]
rownames(data) = 1:nrow(data)

#Cleaner version of data; remove unneeded cols.
Y = data[,c('group','gene','replicate','time','log2exp')]

#Preview data.
head(Y)

#================================================================
# Optimize Parameters Using Max Lhood (Nelder-Mead) =============
#================================================================

#Optimize a set of parameters for each group.  Same initial values for all groups.  Log scale.
ls.f 	= -2
var.f 	= -2
ls.g 	= 1
var.g 	= 1
ls.h 	= 2
var.h 	= 2
sig.sq 	= -2

init.par = c(ls.f,var.f,ls.g,var.g,ls.h,var.h,sig.sq)

#Optimize parameters for each group.
myopt = opt.param(Y,init.par)

#================================================================
# Predict Group, Gene and Rep Profiles ==========================
#================================================================

t.pred=seq(0.5, 12.5, length.out = 200)
myPred = hgp.predict(Y, params=myopt, t.pred)

#================================================================
# Plot Group Profiles ===========================================
#================================================================

h.mean = myPred$h.mean
h.lb = myPred$h.lb
h.ub = myPred$h.ub

gene.cols = rainbow(14)

pdf('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Gene_Expression/Group_Pred_Profiles.pdf',width=15,height=5)

par(mfrow=c(1,3))

for (i in 1:3){
	#Extract data to plot.
	grp = colnames(h.mean)[i]
	Yi  = Y[Y$group==grp,]
	t = Yi$time
	y = Yi$log2exp
	
	#Scatterplot of genes in group.
	plot(t, y, col = gene.cols[Yi$gene],main = paste('Profile for ',grp), pch=19, 
		ylim=c(min(h.lb[,i])-.01, max(h.ub[,i]+.01)),
		xlab='time',ylab='log-squared expression level')
	
	#Shade CI region.
	polygon(c(t.pred,rev(t.pred)), c(h.ub[,i], rev(h.lb[,i])), col='lightgrey', border=NA)
	
	#Add posterior mean and CI for group.
	lines(t.pred,h.mean[,i],lwd=2,col='black')
	lines(t.pred,h.lb[,i],lwd=2,lty=1,col='lightgrey')
	lines(t.pred,h.ub[,i],lwd=2,lty=1,col='lightgrey')
	
	#Re-add points (shaded CI covers them up).
	points(t, y, col = gene.cols[Yi$gene], pch=19)
}
dev.off()

#================================================================
# Plot Gene Profiles for each group =============================
#================================================================

g.mean = myPred$g.mean
g.lb = myPred$g.lb
g.ub = myPred$g.ub

rep.cols = c('firebrick3','forestgreen','darkblue')

for (i in 1:3){
	
	#Extract data to plot.
	grp = colnames(h.mean)[i]
	Yi  = Y[Y$group==grp,]

	#Genes in group. 
	genes.i = unique(Yi$gene)
	
	pdf(trimws(paste('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Gene_Expression/		Gene_Pred_Profiles_',grp,'.pdf',sep='')), width=10,height=10)
	
	par(mfrow=c(2,3))
	
	#Loop through genes in group.
	for (j in 1:length(genes.i)){
		
		gene = unique(Yi$gene)[j]
		Yn = Yi[Yi$gene==gene,]
		t = Yn$time
		y = Yn$log2exp
		
		#Scatterplot of the gene.
		plot(t, y, col = rep.cols[Yn$replicate],main = paste('Gene Profile: ',gene), pch=19, 
			ylim=c(min(g.lb[,gene])-.01, max(g.ub[,gene]+.01)),
			xlab='time',ylab='log-squared expression level')
		
		#Shade CI region.
		polygon(c(t.pred,rev(t.pred)), c(g.ub[,gene], rev(g.lb[,gene])), col='lightgrey', border=NA)
		
		#Add posterior mean and CI for group.
		lines(t.pred,g.mean[,gene],lwd=2,col='black')
		lines(t.pred,g.lb[,gene],lwd=2,lty=1,col='lightgrey')
		lines(t.pred,g.ub[,gene],lwd=2,lty=1,col='lightgrey')
	
		#Re-add points (shaded CI covers them up).
		points(t, y, col = rep.cols[Yn$replicate], pch=19)
	}
	dev.off()
}


#================================================================
# Plot Gene-Rep Profiles ========================================
#================================================================

f.mean = myPred$f.mean
f.lb = myPred$f.lb
f.ub = myPred$f.ub

rep.cols = c('firebrick3','forestgreen','darkblue')

for (i in 1:3){
	
	#Extract data to plot.
	grp = colnames(h.mean)[i]
	Yi  = Y[Y$group==grp,]

	#Genes in group. 
	genes.i = unique(Yi$gene)
	
	#Reps in each gene.
	reps = unique(Yi$replicate)
	
	#Gene-rep combos in each group.
	gene.rep = paste( rep(genes.i,each=length(reps)), rep(reps,length(genes.i)) )
	
	pdf(trimws(paste('/Users/jennstarling/UTAustin/2017S_Stats Modeling 2/Exercise-04/Figures/Gene_Expression/		Gene_Replicate_Profiles_',grp,'.pdf',sep='')), width=10,height=10)
	
	nrow = length(genes.i)
	ncol = length(reps)
	
	par(mfrow=c(nrow,ncol))
	
	#Loop through genes/replicates in group.
	for (j in 1:length(genes.i)){
		
		for (k in 1:length(reps)){
			
			#Extract data for gene/rep combo.
			gene = genes.i[j]
			rep = reps[k]
			gene.rep = colnames(f.mean)[colnames(f.mean)==paste(gene,rep)] #For f.mean lookup.
			
			Ynr = Yi[Yi$gene==gene & Yi$replicate==rep,]
			t = Ynr$time
			y = Ynr$log2exp
			
			#Scatterplot of the gene.
			plot(t, y, col = rep.cols[Ynr$replicate],main = paste('Gene-Rep Profile: ',gene.rep), pch=19, 
				ylim=c(min(f.lb[,gene.rep])-.01, max(f.ub[,gene.rep]+.01)),
				xlab='time',ylab='log-squared expression level')
		
			#Shade CI region.
			polygon(c(t.pred,rev(t.pred)), c(f.ub[,gene.rep], rev(f.lb[,gene.rep])), col='lightgrey', border=NA)
		
			#Add posterior mean and CI for group.
			lines(t.pred,f.mean[,gene.rep],lwd=2,col='black')
			lines(t.pred,f.lb[,gene.rep],lwd=2,lty=1,col='lightgrey')
			lines(t.pred,f.ub[,gene.rep],lwd=2,lty=1,col='lightgrey')
	
			#Re-add points (shaded CI covers them up).
			points(t, y, col = rep.cols[Ynr$replicate], pch=19)
		
		} #End replicate loop.
	} #End gene loop.
	
	dev.off()
}


