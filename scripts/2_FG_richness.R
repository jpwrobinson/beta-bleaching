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

## read dataset
s1<-read.csv(file='data/UVC_richnessdiff_FG.csv')


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

compare(m1, m2, m3)
### m3 is top model
m<-m3

#### SAVE TOP MODEL
save(m, file='results/02_richness_corallivore.Rdata')

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

compare(m1, m2, m3)
### m1 is top model
m<-m1

#### SAVE TOP MODEL
save(m, file='results/02_richness_herbivore.Rdata')

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

compare(m1, m2, m3)
### m3 is top model
m<-m3

#### SAVE TOP MODEL
save(m, file='results/02_richness_invertivore.Rdata')

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

compare(m1, m2, m3)
### m3 is top model
m<-m3

#### SAVE TOP MODEL
save(m, file='results/02_richness_mixed-diet.Rdata')


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

compare(m1, m2, m3)
### m2 is top model
m<-m2

#### SAVE TOP MODEL
save(m, file='results/02_richness_planktivore.Rdata')

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

compare(m1, m2, m3)
### m1 is top model
m<-m1

#### SAVE TOP MODEL
save(m, file='results/02_richness_piscivore.Rdata')


## end of script