library(here)
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(rethinking)
library(forcats)
library(ggpubr)
library(cowplot)
library(grid)
library(ggridges)
library(funk)


## load models and rename by feeding group
load(file='results/02_richness_planktivore.Rdata'); m <- plank
load(file='results/02_richness_piscivore.Rdata'); m <- cor
load(file='results/02_richness_herbivore.Rdata'); m <- herb
load(file='results/02_richness_invertivore.Rdata'); m <- inv
load(file='results/02_richness_mixed-diet.Rdata'); m <- mix
load(file='results/02_richness_piscivore.Rdata'); m <- pisc

## ---------------------------------------- ##
  			# Read observed data #
## ---------------------------------------- ##

s1<-read.csv(file='data/UVC_richnessdiff_FG.csv')

## get mean and SE
s1.mean<-s1 %>% group_by(Year, FG, state) %>% summarise(mean=mean(diff), se = se(diff))

## assign plotting info
s1.mean$xlab<-focal.scaled$Year[match(s1.mean$Year, focal.scaled$yearraw)]
s1.mean$xlab<-ifelse(s1.mean$state=='Recovering', s1.mean$xlab-0.01, s1.mean$xlab+0.01)
s1.mean$cols2<-ifelse(s1.mean$state=='Shifted', shif.col, r.col)
s1.mean$pch<-ifelse(s1.mean$state=='Shifted', 21, 24)

## ---------------------------------------- ##
  # Make posterior predictions for each FG #
## ---------------------------------------- ##

# replace varying intercept samples with zeros
# 1000 samples by 21 sites
a_location_zeros <- matrix(0,1000,21)
bb_location_zeros <- matrix(0,1000,21)
bd_location_zeros <- matrix(0,1000,21)

#index for separating states
rec.vec<-seq(2,60, by=2)
shif.vec<-seq(1,59, by=2)

## ---------- ##
## Herbivores ##
## ---------- ##
focal<-s1 %>% filter(Year != 1994 & FG == 'Herbivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))
focal.scaled$year.raw<-rep(unique(focal$Year), each=21)[-c(99:101)]
## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

# make predictions
herb.dat<-list(Year = rep(seq(min(focal.scaled$Year), max(focal.scaled$Year), length.out=30),each=2), 
				RecoveringShifteddummy=rep(c(1,0), times=30),
				Location = rep(factor('Cousin Carbonate'), times=60))
means<-link(herb, data=herb.dat, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
herb.mean<-apply(means, 2, meean)
## mean intervals and 89% uncertainty
herb.PI89 <- apply( means , 2 , HPDI, prob=0.95 )

## ---------- ##
## Corallivores ##
## ---------- ##
focal<-s1 %>% filter(Year != 1994 & FG == 'Corallivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))
focal.scaled$year.raw<-rep(unique(focal$Year), each=21)[-c(99:101)]
## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

# make predictions
cor.dat<-list(Year = rep(seq(min(focal.scaled$Year), max(focal.scaled$Year), length.out=30),each=2), 
				RecoveringShifteddummy=rep(c(1,0), times=30),
				Location = rep(factor('Cousin Carbonate'), times=60))
means<-link(cor, data=cor.dat, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
cor.mean<-apply(means, 2, mean)
## mean intervals and 89% uncertainty
cor.PI89 <- apply( means , 2 , HPDI, prob=0.95 )

## ---------- ##
## Invertivores ##
## ---------- ##
focal<-s1 %>% filter(Year != 1994 & FG == 'Invertivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))
focal.scaled$year.raw<-rep(unique(focal$Year), each=21)[-c(99:101)]
## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

# make predictions
inv.dat<-list(Year = rep(seq(min(focal.scaled$Year), max(focal.scaled$Year), length.out=30),each=2), 
				RecoveringShifteddummy=rep(c(1,0), times=30),
				Location = rep(factor('Cousin Carbonate'), times=60))
means<-link(inv, data=inv.dat, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
## mean intervals and 89% uncertainty
inv.mean<-apply(means, 2, mean)
inv.PI89 <- apply( means , 2 , HPDI, prob=0.95 )

## ---------- ##
## Mixed-diet ##
## ---------- ##
focal<-s1 %>% filter(Year != 1994 & FG == 'Mixed-diet Feeder') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))
focal.scaled$year.raw<-rep(unique(focal$Year), each=21)[-c(99:101)]
## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

# make predictions
mix.dat<-list(Year = rep(seq(min(focal.scaled$Year), max(focal.scaled$Year), length.out=30),each=2), 
				RecoveringShifteddummy=rep(c(1,0), times=30),
				Location = rep(factor('Cousin Carbonate'), times=60))
means<-link(mix, data=mix.dat, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
## mean intervals and 89% uncertainty
mix.mean<-apply(means, 2, mean)
mix.PI89 <- apply( means , 2 , HPDI, prob=0.95 )

## ---------- ##
## Piscivores ##
## ---------- ##
focal<-s1 %>% filter(Year != 1994 & FG == 'Piscivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))
focal.scaled$year.raw<-rep(unique(focal$Year), each=21)[-c(99:101)]
## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

# make predictions
pisc.dat<-list(Year = rep(seq(min(focal.scaled$Year), max(focal.scaled$Year), length.out=30),each=2), 
				RecoveringShifteddummy=rep(c(1,0), times=30),
				Location = rep(factor('Cousin Carbonate'), times=60))
means<-link(pisc, data=pisc.dat, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
## mean intervals and 89% uncertainty
pisc.mean<-apply(means, 2, mean)
pisc.PI89 <- apply( means , 2 , HPDI, prob=0.95 )

## ---------- ##
## Planktivores ##
## ---------- ##
focal<-s1 %>% filter(Year != 1994 & FG == 'Planktivore') %>% select(Year, diff, state, Location)
focal.scaled<-scaler(focal, ID=c('diff', 'Location'))
focal.scaled$year.raw<-rep(unique(focal$Year), each=21)[-c(99:101)]
## remove dots for stan
colnames(focal.scaled)<-str_replace_all(colnames(focal.scaled), '\\.', '\\')

# make predictions
plank.dat<-list(Year = rep(seq(min(focal.scaled$Year), max(focal.scaled$Year), length.out=30),each=2), 
				RecoveringShifteddummy=rep(c(1,0), times=30),
				Location = rep(factor('Cousin Carbonate'), times=60))
means<-link(plank, data=plank.dat, replace=list(a_loc=a_location_zeros, bb_loc=bb_location_zeros, bd_loc=bd_location_zeros))
means<-means$mu
## mean intervals and 89% uncertainty
plank.mean<-apply(means, 2, mean)
plank.PI89 <- apply( means , 2 , HPDI, prob=0.95 )

## ---------------------------------------- ##
  			  # Create Figure #
## ---------------------------------------- ##

## plotting format stuff
r.col<-c('#045a8d')	
shif.col<-c('#b30000')

fg.cols<-data.frame(cols<-c('#e41a1c','#4daf4a','#a65628','#ff7f00', '#377eb8', '#984ea3'),
		FG=c('Corallivore', 'Herbivore', 'Invertivore', 'Mixed-diet Feeder',  'Piscivore','Planktivore'))

barW=0.2
mat<-matrix(c(1,1,3,3,5,5,1,1,3,3,5,5,2,2,4,4,6,6,2,2,4,4,6,6), nrow=6, ncol=4, byrow=F)

pdf(file='figures/Figure2.pdf', height=6, width=8)

layout(mat)


par(mar=c(1,4,1,2))
## Panel A: richness turnover - Planktivore
with(s1.mean[s1.mean$FG=='Planktivore',],
plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=21, xlab='', ylab='', axes=F, ylim=c(-3,2), xlim=c(-1.4, 1.5),
	pt.bg='transparent', col='transparent', cex=1.5, sfrac=0, scol='white', bty='n'))
axis(2); mtext(2, text=expression(paste(Delta, ' richness'[1994])),line=2.5, cex=0.8)
axis(1, at=unique(focal.scaled$Year), labels=NA)
par(xpd=F); abline(a = 0, b=0, lty=2)
add_label(0.01, 0.05, label='(a) Planktivore', font=2 ,cex=1)

## post predictions
lines( plank.dat$Year[rec.vec] , plank.mean[rec.vec] , lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Planktivore']))
lines( plank.dat$Year[shif.vec] , plank.mean[shif.vec] , lty=2,lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Planktivore']))
shade( plank.PI89[,rec.vec] , plank.dat$Year[rec.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Planktivore']), 0.2))
shade( plank.PI89[,shif.vec] , plank.dat$Year[shif.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Planktivore']), 0.2))

## add observed points
with(s1.mean[s1.mean$FG=='Planktivore',],
	plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=pch, add=TRUE, pt.bg='black', col='transparent', cex=1.2, sfrac=0, scol='darkgrey'))

legend('topright', legend=c('Recovering', 'Shifted'), pch=c(24,21), pt.bg='black', bty='n', cex=0.8)
legend('topright', legend=c('', ''), lty=c(1,2), inset=c(0.18, 0), col='black', bty='n', cex=0.8)


## Panel B: richness turnover - Corallivore
par(mar=c(1,1,1,3))
with(s1.mean[s1.mean$FG=='Corallivore',],
plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=21, xlab='', ylab='', axes=F, ylim=c(-6,2), xlim=c(-1.4, 1.5),
	pt.bg='transparent', col='transparent', cex=1.5, sfrac=0, scol='white', bty='n'))
axis(2); #mtext(2, text=expression(paste(Delta, ' richness'[1994])),line=2.5, cex=0.8)
axis(1, at=unique(focal.scaled$Year), labels=NA)
par(xpd=F); abline(a = 0, b=0, lty=2)
add_label(0.01, 0.05, label='(b) Corallivore', font=2 ,cex=1)

## post predictions
lines( cor.dat$Year[rec.vec] , cor.mean[rec.vec], lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Corallivore']))
lines( cor.dat$Year[shif.vec] , cor.mean[shif.vec] ,lty=2, lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Corallivore']))
shade( cor.PI89[,rec.vec] , cor.dat$Year[rec.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Corallivore']), 0.2))
shade( cor.PI89[,shif.vec] , cor.dat$Year[shif.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Corallivore']), 0.2))

## add observed points
with(s1.mean[s1.mean$FG=='Corallivore',],
	plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=pch, add=TRUE, pt.bg='black', col='transparent', cex=1.2, sfrac=0, scol='darkgrey'))


## Panel C: richness turnover - Herbivore
par(mar=c(1,4,1,2))
with(s1.mean[s1.mean$FG=='Herbivore',],
plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=21, xlab='', ylab='', axes=F, ylim=c(-6,8), xlim=c(-1.4, 1.5),
	pt.bg='transparent', col='transparent', cex=1.5, sfrac=0, scol='white', bty='n'))
axis(1, at=unique(focal.scaled$Year), labels=NA); 
axis(2); mtext(2, text=expression(paste(Delta, ' richness'[1994])),line=2.5, cex=0.8)
par(xpd=F); abline(a = 0, b=0, lty=2)
add_label(0.01, 0.05, label='(c) Herbivore', font=2 ,cex=1)

## post predictions
lines( herb.dat$Year[rec.vec] , herb.mean[rec.vec] , lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Herbivore']))
lines( herb.dat$Year[shif.vec] , herb.mean[shif.vec] ,lty=2, lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Herbivore']))
shade( herb.PI89[,rec.vec] , cor.dat$Year[rec.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Herbivore']), 0.2))
shade( herb.PI89[,shif.vec] , cor.dat$Year[shif.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Herbivore']), 0.2))

## add observed points
with(s1.mean[s1.mean$FG=='Herbivore',],
	plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=pch, add=TRUE, pt.bg='black', col='transparent', cex=1.2, sfrac=0, scol='darkgrey'))


## Panel D: richness turnover - Invertivore
par(mar=c(1,1,1,3))
with(s1.mean[s1.mean$FG=='Invertivore',],
plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=21, xlab='', ylab='', axes=F, ylim=c(-8,4), xlim=c(-1.4, 1.5),
	pt.bg='transparent', col='transparent', cex=1.5, sfrac=0, scol='white', bty='n'))
axis(1, at=unique(focal.scaled$Year), labels=NA)
axis(2); #mtext(2, text=expression(paste(Delta, ' richness'[1994])),line=2.5, cex=0.8)
par(xpd=F); abline(a = 0, b=0, lty=2)
add_label(0.01, 0.05, label='(d) Invertivore', font=2 ,cex=1)

## post predictions
lines( inv.dat$Year[rec.vec] , inv.mean[rec.vec] ,lty=1, lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Invertivore']))
lines( inv.dat$Year[shif.vec] , inv.mean[shif.vec] ,lty=2, lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Invertivore']))
shade( inv.PI89[,rec.vec] , cor.dat$Year[rec.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Invertivore']), 0.2))
shade( inv.PI89[,shif.vec] , cor.dat$Year[shif.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Invertivore']), 0.2))

## add observed points
with(s1.mean[s1.mean$FG=='Invertivore',],
	plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=pch, add=TRUE, pt.bg='black', col='transparent', cex=1.2, sfrac=0, scol='darkgrey'))



## Panel E: richness turnover - Mixed-diet Feeder
par(mar=c(2,4,0,2))
with(s1.mean[s1.mean$FG=='Mixed-diet Feeder',],
plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=21, xlab='', ylab='', axes=F, ylim=c(-6.5,4.5), xlim=c(-1.4, 1.5),
	pt.bg='transparent', col='transparent', cex=1.5, sfrac=0, scol='white', bty='n'))
axis(2); mtext(2, text=expression(paste(Delta, ' richness'[1994])),line=2.5, cex=0.8)
axis(1, at=unique(focal.scaled$Year), labels=c(2005, 2008, 2011, 2014, 2017))
par(xpd=F); abline(a = 0, b=0, lty=2)
add_label(0.01, 0.05, label='(e) Mixed-diet Feeder', font=2 ,cex=1)

## post predictions
lines( mix.dat$Year[rec.vec] , mix.mean[rec.vec] , lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Mixed-diet Feeder']))
lines( mix.dat$Year[shif.vec] , mix.mean[shif.vec] , lty=2,lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Mixed-diet Feeder']))
shade( mix.PI89[,rec.vec] , mix.dat$Year[rec.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Mixed-diet Feeder']), 0.2))
shade( mix.PI89[,shif.vec] , mix.dat$Year[shif.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Mixed-diet Feeder']), 0.2))

## add observed points
with(s1.mean[s1.mean$FG=='Mixed-diet Feeder',],
	plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=pch, add=TRUE, pt.bg='black', col='transparent', cex=1.2, sfrac=0, scol='darkgrey'))


## Panel F: richness turnover - Piscivore
par(mar=c(2,1,0,3))
with(s1.mean[s1.mean$FG=='Piscivore',],
plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=21, xlab='', ylab='', axes=F, ylim=c(-3,2), xlim=c(-1.4, 1.5),
	pt.bg='transparent', col='transparent', cex=1.5, sfrac=0, scol='white', bty='n'))
axis(2); #mtext(2, text=expression(paste(Delta, ' richness'[1994])),line=2.5, cex=0.8)
axis(1, at=unique(focal.scaled$Year), labels=c(2005, 2008, 2011, 2014, 2017))
par(xpd=F); abline(a = 0, b=0, lty=2)
add_label(0.01, 0.05, label='(f) Piscivore', font=2 ,cex=1)

## post predictions
lines( pisc.dat$Year[rec.vec] , pisc.mean[rec.vec] , lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Piscivore']))
lines( pisc.dat$Year[shif.vec] , pisc.mean[shif.vec] ,lty=2, lwd=2,col=as.character(fg.cols$col[fg.cols$FG=='Piscivore']))
shade( pisc.PI89[,rec.vec] , pisc.dat$Year[rec.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Piscivore']), 0.2))
shade( pisc.PI89[,shif.vec] , pisc.dat$Year[shif.vec] ,col=alpha(as.character(fg.cols$col[fg.cols$FG=='Piscivore']), 0.2))

## add observed points
with(s1.mean[s1.mean$FG=='Piscivore',],
	plotCI(x=xlab, y=mean, ui=mean+se, li=mean-se,pch=pch, add=TRUE, pt.bg='black', col='transparent', cex=1.2, sfrac=0, scol='darkgrey'))



dev.off()


## end of script
