#!/bin/env Rscript


library(dplyr)
library(stringr)
library(rethinking)
library(funk)


## define Bayesian sampling params
iter = 7000
cores = 3
warmup = 1500
chains = 3

############ M5 = biomass change from 1994-2014 ####################
	
biom20<-read.csv(file='data/UVC_biom_change.csv')

## remove species that did not appear in either 1994 or later, per site
biom20<-biom20[!(biom20$biom==0 & biom20$biom94==0),]


## add benthic predictors
load('data/SEY_UVC_benthic.Rdata')
SC<-SC[SC$year=='2014',]
biom20$coralmassives<-SC$coral.massives[match(biom20$Location, SC$location)]
biom20$coralencrusting<-SC$coral.encrusting[match(biom20$Location, SC$location)]
biom20$coralbranching<-SC$coral.branching[match(biom20$Location, SC$location)]
biom20$macroalgae<-SC$macroalgae[match(biom20$Location, SC$location)]
biom20$complexity<-SC$complexity[match(biom20$Location, SC$location)]
biom20$Management<-fish$Management[match(biom20$Location, fish$Location)]

## add dummys and scale predictors
biomS<-scaler(biom20, ID=c('Location', 'ID','Species', 'Family', 'FG', 'FG.fine' ,'biom.change'))

## remove dots for stan
colnames(biomS)<-str_replace_all(colnames(biomS), '\\.', '\\')
## remove spaces for stan
biomS$FG<-str_replace_all(biomS$FG, '\\ ', '')
biomS$FG<-str_replace_all(biomS$FG, '\\-', '')
biomS$FG<-as.factor(biomS$FG)


## Fit bayesian model to biomass change
## m = fixed covariates as in beta models + bodysize + random_ints[location] + random_slopes[FG]
m <- map2stan(
	alist(
		# likelihood for y, normally distributed biomass differences
	    biomchange ~ dnorm( mu , sigma ),

	    # linear model
	    mu <- A +
	    	  BA*RecoveringShifteddummy +
			  BF*FishedProtecteddummy +
			  BG*coralbranching +
			  BH*coralmassives +
			  BI*coralencrusting +
			  BJ*macroalgae +
			  BK*complexity +
			  BL*bodysize,
		# intercept, varying by sites (Location) and functional groups (FG)
	    A <- a + a_loc[Location] + a_fg[FG],
	    # fixed covariate parameters, regime effect varying by FG
		BA <- ba + ba_fg[FG],
		BF <- bf, 
		BG <- bg, 
		BH <- bh, 
		BI <- bi, 
		BJ <- bj, 
		BK <- bk, 
		BL <- bl, 
		
		## adaptive priors
		## Location random intercept
		c(a_loc)[Location] ~ dnorm(0, sigma_loc),
		## Functional group random slopes
		c(a_fg, ba_fg)[FG] ~ dmvnorm2(0, sigma_fg, rho_fg),
		
		## fixed priors	
	    a ~ dnorm(0, 10), ## set to 0 biomass change
	    c(ba, bf, bg, bh, bi, bj, bk, bl) ~ dnorm(0, 10),
	    c(sigma, sigma_loc, sigma_fg) ~ dcauchy(0, 2),
	    rho_fg ~ dlkjcorr(2)
		), 
	data=biomS, iter=iter, warmup = warmup, chains=chains, cores=cores)


save(biomS, biom20, m, file='results/05_biom_change_model.Rdata')


## end of script