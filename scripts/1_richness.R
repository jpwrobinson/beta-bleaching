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


############ Temporal richness models  ####################

beta<-read.csv('data/UVC_temporal_richness.csv')
beta$site.year<-paste(beta$Location, beta$Year, sep='.')

## add predictors
load('data/SEY_UVC_benthicPV_SC_DEPTH.Rdata')
beta$hardcoral<-SC$hard.coral[match(beta$site.year, SC$site.year)]
beta$macroalgae<-SC$macroalgae[match(beta$site.year, SC$site.year)]
beta$complexity<-SC$complexity[match(beta$site.year, SC$site.year)]
beta$coralmassives<-SC$coral.massives[match(beta$site.year, SC$site.year)]
beta$coralencrusting<-SC$coral.encrusting[match(beta$site.year, SC$site.year)]
beta$coralbranching<-SC$coral.branching[match(beta$site.year, SC$site.year)]
beta <- beta %>% filter(Year != 1994)


## add dummys and scale predictors
betaS<-scaler(beta, ID=c('Location', 'site.year', 'richness'))
## remove dots for stan
colnames(betaS)<-str_replace_all(colnames(betaS), '\\.', '\\')


## Check non linear year relationships
## With Bayesian - polynomial year effects
### mu is richness, sigmar is variance

## m1 = Linear year effect
m1 <- map2stan(
	alist(
	    richness ~ dnorm( mu , sigma ) ,
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
	    a ~ dnorm(53, 10), ## set to mean 1994 richness
	    c(ba, bb, bc) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)


## m2 = Non-linear year effect, X^2
m2 <- map2stan(
	alist(
	    richness ~ dnorm( mu , sigma ) ,
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
	    a ~ dnorm(53, 10), ## set to mean 1994 richness
	    c(ba, bb, bc, bd, be) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)


## m3 = Non-linear year effect, X^2 + X^3
m3 <- map2stan(
	alist(
	    richness ~ dnorm( mu , sigma ) ,
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
	    a ~ dnorm(53, 10), ## set to mean 1994 richness
	    c(ba, bb, bc, bd, be, bf, bg) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup = warmup, chains=chains, cores=cores)


## Now add explanatory covariates

## m4 = Non-linear year effect, X^2 + X^3, and explanatory covariates
m4 <- map2stan(
	alist(
	    richness ~ dnorm( mu , sigma ) ,
	    mu <- A +
	    	  BA*RecoveringShifteddummy +
	    	  BB*Year +
	    	  BC*RecoveringShifteddummy*Year +
	    	  BD*Year^2 +
	    	  BE*RecoveringShifteddummy*Year^2 +
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
		BD <- bd + bd_loc[Location],
		BE <- be,
		BF <- bf,
		BG <- bg,
		BH <- bh,
		BI <- bi,
		BJ <- bj,
		BK <- bk,
		## adaptive priors
		c(a_loc, bb_loc, bd_loc)[Location] ~ dmvnorm2(0, sigma_loc, rho_loc),
		## fixed priors	
	    a ~ dnorm(53, 10), ## set to mean 1994 richness
	    c(ba, bb, bc, bd, be, bf, bg, bh, bi, bj, bk) ~ dnorm(0, 10),
	    c(sigma, sigma_loc) ~ dcauchy(0, 2),
	    rho_loc ~ dlkjcorr(2)
), data=betaS, iter=iter, warmup. = warmup, chains=chains, cores=cores)
							

save(betaS, m1,m2,m3,m4, file='results/richness_mods_remote.Rdata')


sink()
