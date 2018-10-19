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

############ Beta seq = temporal beta diversity between years ####################

beta<-read.csv(file='data/UVC_beta_seq.csv')

## add predictors
load('data/SEY_UVC_benthic.Rdata')
beta$hardcoral<-SC$hard.coral[match(beta$site.year, SC$site.year)]
beta$macroalgae<-SC$macroalgae[match(beta$site.year, SC$site.year)]
beta$complexity<-SC$complexity[match(beta$site.year, SC$site.year)]
beta$coral.massives<-SC$coral.massives[match(beta$site.year, SC$site.year)]
beta$coral.encrusting<-SC$coral.encrusting[match(beta$site.year, SC$site.year)]
beta$coral.branching<-SC$coral.branching[match(beta$site.year, SC$site.year)]

## add dummys and scale predictors
betaS<-scaler(beta, ID=c('Location', 'site.year', 'beta.bray', 'Location'))
## remove dots for stan
colnames(betaS)<-str_replace_all(colnames(betaS), '\\.', '\\')

## Check non linear year relationships
## With Bayesian - polynomial year effects
### mu is betabray, sigmar is variance
## Trends fitted separately to each UVC site as hierarchical structure

## m1 = Linear year effect
m1 <- map2stan(
	alist(
	    betabray ~ dnorm( mu , sigma ) ,
	    mu <- A + BA*RecoveringShifteddummy +
	    	  BB*Year +
	    	  BC*RecoveringShifteddummy*Year,
	    A <- a + a_loc[Location],
		BA <- ba,
		BB <- bb + bb_loc[Location],
		BC <- bc,
		## adaptive priors
		c(a_loc, bb_loc)[Location] ~ dmvnorm2(0, sigma_loc, rho_loc),
		## fixed priors	
	    a ~ dnorm(0, 10), ## set to null expectation (perfect similarity)
	    c(ba, bb, bc) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)


## m2 = Non-linear year effect, X^2
m2 <- map2stan(
	alist(
	    betabray ~ dnorm( mu , sigma ) ,
	    mu <- A +
	    	  BA*RecoveringShifteddummy +
	    	  BB*Year +
	    	  BC*RecoveringShifteddummy*Year +
	    	  BD*Year^2 +
	    	  BE*RecoveringShifteddummy*Year^2,
	    A <- a + a_loc[Location],
		BA <- ba,
		BB <- bb + bb_loc[Location],
		BC <- bc,
		BD <- bd + bd_loc[Location],
		BE <- be,
		## adaptive priors
		c(a_loc, bb_loc, bd_loc)[Location] ~ dmvnorm2(0, sigma_loc, rho_loc),
		## fixed priors	
	    a ~ dnorm(0, 10), ## set to null expectation (perfect similarity)
	    c(ba, bb, bc, bd, be) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)


## m3 = Non-linear year effect, X^2 + X^3
m3 <- map2stan(
	alist(
	    betabray ~ dnorm( mu , sigma ) ,
	    mu <- A + BA*RecoveringShifteddummy +
	    		  BB*Year +
	    		  BC*RecoveringShifteddummy*Year +
	    		  BD*Year^2 +
	    		  BE*RecoveringShifteddummy*Year^2 +
	    		  BF*Year^3 +
	    		  BG*RecoveringShifteddummy*Year^3,
	    A <- a + a_loc[Location],
		BA <- ba,
		BB <- bb + bb_loc[Location],
		BC <- bc,
		BD <- bd + bd_loc[Location],
		BE <- be,
		BF <- bf + bf_loc[Location],
		BG <- bg,
		## adaptive priors
		c(a_loc, bb_loc, bd_loc, bf_loc)[Location] ~ dmvnorm2(0, sigma_loc, rho_loc),
		## fixed priors	
	    a ~ dnorm(0, 10), ## set to null expectation (perfect similarity)
	    c(ba, bb, bc, bd, be, bf, bg) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)

## compare model fits
compare(m1, m2, m3)
## m1 is top model

## Now add explanatory covariates

## m4 = linear year effect and explanatory covariates
m4 <- map2stan(
	alist(
	    betabray ~ dnorm( mu , sigma ) ,
	    mu <- A +
	    	  BA*RecoveringShifteddummy +
	    	  BB*Year +
	    	  BC*RecoveringShifteddummy*Year +
	    	  # BD*Year^2 +
	    	  # BE*RecoveringShifteddummy*Year^2 +
			  BF*FishedProtecteddummy +
			  BG*coralbranching +
			  BH*coralmassives +
			  BI*coralencrusting +
			  BJ*macroalgae +
			  BK*complexity,
	    A <- a + a_loc[Location],
		BA <- ba,
		BB <- bb + bb_loc[Location],
		BC <- bc,
		# BD <- bd + bd_loc[Location],
		# BE <- be,
		BF <- bf,
		BG <- bg,
		BH <- bh,
		BI <- bi,
		BJ <- bj,
		BK <- bk,
		## adaptive priors
		c(a_loc, bb_loc)[Location] ~ dmvnorm2(0, sigma_loc, rho_loc),
		## fixed priors	
	    a ~ dnorm(0, 10), ## set to null expectation (perfect similarity)
	    c(ba, bb, bc, bf, bg, bh, bi, bj, bk) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)
							
## make posterior predictions for shif and recovering sites
shif.vars<-with(betaS[betaS$RecoveringShifteddummy==1,], 
			data.frame(coralbranching = mean(coralbranching),
						coralencrusting = mean(coralencrusting),  
						coralmassives = mean(coralmassives), 
						macroalgae = mean(macroalgae), 
						complexity = mean(complexity)))

rec.vars<-with(betaS[betaS$RecoveringShifteddummy==0,], 
			data.frame(coralbranching = mean(coralbranching),
						coralencrusting = mean(coralencrusting),  
						coralmassives = mean(coralmassives), 
						macroalgae = mean(macroalgae), 
						complexity = mean(complexity)))

## dataset to simulate over. Location is a placeholder - will be removed with random effects.
pred.dat<-list( Year = rep(seq(min(betaS$Year),max(betaS$Year), length.out=12),times=2),
				RecoveringShifteddummy=rep(c(1,0),each=12),
				Location = rep(factor('Cousin Carbonate'), times=24),
				FishedProtecteddummy=rep(mean(betaS$FishedProtecteddummy),24),
				coralencrusting = rep(c(shif.vars$coralencrusting, rec.vars$coralencrusting),each=12),
				coralmassives = rep(c(shif.vars$coralmassives, rec.vars$coralmassives),each=12),
				coralbranching = rep(c(shif.vars$coralbranching, rec.vars$coralbranching),each=12),
				macroalgae = rep(c(shif.vars$macroalgae, rec.vars$macroalgae),each=12),
				complexity = rep(c(shif.vars$complexity, rec.vars$complexity), each=12))


# replace varying intercept samples with zeros
# 1000 samples by 21 sites
a_location_zeros <- matrix(0,1000,21)
bb_location_zeros <- matrix(0,1000,21)
bd_location_zeros <- matrix(0,1000,21)

means<-link(m4, data=pred.dat, n=1000, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
pred.mean<-apply(means, 2, mean)
pred.PI66 <- apply( means , 2 , HPDI , prob=0.66 )
pred.PI89 <- apply( means , 2 , HPDI , prob=0.89 )


save(betaS, m4, pred.mean, pred.PI66, pred.PI89, file='results/04_beta_seq_model.Rdata')

## end of script