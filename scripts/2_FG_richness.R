#!/bin/env Rscript

library(here)
library(dplyr)
library(stringr)
library(rethinking)
setwd(here('beta-bleaching'))
source('scaling_function.R')

## save output for later inspection
sink.file = paste0("texts/2A_FG_output.txt")
sink(sink.file)


s1<-read.csv(file='data/UVC_richnessdiff_from94_FG.csv')

## define Bayesian sampling params
iter = 7000
cores = 3
warmup = 1500
chains = 3

## Define model structures

## LINEAR RESPONSES
m1.structure <- alist(
      diff ~ dnorm( mu , sigma ) ,
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
      a ~ dnorm(0, 10), ## set to 0 change in richness
      c(ba, bb, bc) ~ dnorm(0, 10),
      c(sigma, sigma_loc) ~ dcauchy(0, 2),
      rho_loc ~ dlkjcorr(2))

## NON-LINEAR RESPONSES
## m2 = Non-linear year effect, X^2
m2.structure <- alist(
      diff ~ dnorm( mu , sigma ) ,
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
      rho_loc ~ dlkjcorr(2))

## m3 = Non-linear year effect, X^2 + X^3
m3.structure <- alist(
      diff ~ dnorm( mu , sigma ) ,
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
      rho_loc ~ dlkjcorr(2))

## ------- ------- ------- ------- ------- ------- ------- ##
              ### Fit models to each functional group ###
## ------- ------- ------- ------- ------- ------- ------- ##

## ------- ------- ------- ------- ------- ##
              ### Corallivore ###
## ------- ------- ------- ------- ------- ##
focal<-s1 %>% filter(FG == 'Corallivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))

## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

m1 <- map2stan(m1.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m2 <- map2stan(m2.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m3 <- map2stan(m3.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)

#### SAVE TOP MODELS
save(m1,m2,m3, file='results/richness_corallivore_remote.Rdata')

## ------- ------- ------- ------- ------- ##
              ### Herbivore ###
## ------- ------- ------- ------- ------- ##
focal<-s1 %>% filter(FG == 'Herbivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))

## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

m1 <- map2stan(m1.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m2 <- map2stan(m2.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m3 <- map2stan(m3.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)

#### SAVE TOP MODELS
save(m1,m2,m3, file='results/richness_herbivore_remote.Rdata')

## ------- ------- ------- ------- ------- ##
              ### Invertivore ###
## ------- ------- ------- ------- ------- ##
focal<-s1 %>% filter(FG == 'Invertivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))

## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

m1 <- map2stan(m1.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m2 <- map2stan(m2.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m3 <- map2stan(m3.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)

#### SAVE TOP MODELS
save(m1,m2,m3, file='results/richness_invertivore_remote.Rdata')

## ------- ------- ------- ------- ------- ##
              ### Mixed-diet Feeder ###
## ------- ------- ------- ------- ------- ##
focal<-s1 %>% filter(FG == 'Mixed-diet Feeder') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))

## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

m1 <- map2stan(m1.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m2 <- map2stan(m2.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m3 <- map2stan(m3.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)

#### SAVE TOP MODELS
save(m1,m2,m3, file='results/richness_mixed-diet_remote.Rdata')


## ------- ------- ------- ------- ------- ##
              ### Planktivore ###
## ------- ------- ------- ------- ------- ##
focal<-s1 %>% filter(FG == 'Planktivore') %>% mutate(yr = scale(Year)) %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))

## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

m1 <- map2stan(m1.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m2 <- map2stan(m2.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m3 <- map2stan(m3.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)

#### SAVE TOP MODELS
save(m1,m2,m3, file='results/richness_planktivore_remote.Rdata')

## ------- ------- ------- ------- ------- ##
              ### Piscivore ###
## ------- ------- ------- ------- ------- ##
focal<-s1 %>% filter(FG == 'Piscivore') %>% mutate(yr = scale(Year)) %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))

## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

m1 <- map2stan(m1.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m2 <- map2stan(m2.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)
m3 <- map2stan(m3.structure, data=focal.scaled, iter=iter, warmup = warmup, chains=chains, cores=cores)

#### SAVE TOP MODELS
save(m1,m2,m3, file='results/richness_piscivore_remote.Rdata')


sink()