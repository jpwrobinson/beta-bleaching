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


## plotting information
recovering<-c('#045a8d')	
shifted<-c('#b30000')
theme_set(theme_sleek())


pdf(file='figures-pdf/Figure1_ab.pdf', height=4, width=10)

## -------------------------------------------------------------------------- ##
## --------------- predicted richness change from 1994 ---------------------- ##
## --------------------------------------------------------------------------- ##

beta<-read.csv('data/UVC_richness.csv')
r94<-beta %>% filter(Year == 1994) %>% mutate(year.plot=2003.5)
beta<-beta %>% filter(Year != 1994)
beta$pt.col<-ifelse(beta$state == 'Recovering', recovering, shifted)

## load posteriors
load(file='results/01_richness_model.Rdata')


pred.clean<-data.frame(richness = c(pred.mean[c(1:12)],pred.mean[c(13:24)]),
						state = rep(c('Shifted', 'Recovering'), each=12),
						Year = rep(seq(min(beta$Year),max(beta$Year), length.out=12),times=2))

pred.clean$PI95.lower<-c(pred.PI95[1,c(1:12)], pred.PI95[1,c(13:24)])
pred.clean$PI95.upper<-c(pred.PI95[2,c(1:12)], pred.PI95[2,c(13:24)])

bs0<-ggplot(pred.clean, aes(Year, richness, fill=state, group=state, col=state)) + 
		geom_ribbon(alpha=0.1, show.legend=F,linetype=0,aes(ymin = PI95.lower, ymax=PI95.upper)) +
		geom_line(lwd=1, aes(linetype=state)) + geom_vline(xintercept=2004.25, linetype='dashed', col='grey')+
		scale_fill_manual(values=c(recovering,shifted)) +
		scale_linetype_manual(values=c(1,5)) +
		scale_colour_manual(values=c(recovering,shifted)) +
		geom_boxplot(data=r94, lwd=0.25,outlier.size=0.5, show.legend=F, position=position_dodge(0.8),
			aes(year.plot, richness, group=state, col=factor(state),fill=factor(state))) +
		stat_summary(data=r94, geom = "crossbar", width=0.7, show.legend=F, fatten=0, position=position_dodge(0.8),color="white",
			aes(year.plot, richness, fill=state, group=state),  fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
		theme(axis.text.x=element_text(size=10, colour='black'),
			axis.text.y=element_text(size=10, colour='black'),
			axis.title.y=element_text(size=12, colour='black'),
				panel.border = element_blank(),
				 panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
		 		axis.line = element_line(colour = "black"),
				legend.position=c(0.6, 0.99), 
				legend.title=element_blank(), 
				legend.text=element_text(size=9),
				legend.spacing.x = unit(0.25, 'cm'),
				legend.key.width = unit(1.1,"cm")) +
		guides(fill=guide_legend(nrow=1)) +
		scale_y_continuous(limits=c(20, 75)) +
		scale_x_continuous(breaks=c(2003.5,2005, 2008, 2011, 2014, 2017), labels=c(1994, 2005, 2008, 2011, 2014, 2017), limits=c(2003, 2017.5)) +
		labs(x = '', y ='Species richness') + 
		geom_jitter(data = beta, aes(Year, richness), width=0.25, size=1, alpha=0.3)


## now parameter estimates
ylabs<-c('Shifted regime', 'Year','Regime shift * Year', expression('Year'^2), expression('Shifted regime * Year'^2),
				 'Protection', 'Branching coral', 'Massive coral', 'Encrusting coral', 'Macroalgae', 'Complexity')
post<-as.data.frame(extract.samples(m4)) %>% gather(param, dist) %>% 
		filter(param %in% c('ba', 'bb', 'bc', 'bd', 'be', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk'))

params0<-ggplot(post, aes(x = dist, y = param)) + 
		geom_halfeyeh(size=0.5, .width=0.95, fill=NA, density.color=NA) + 
		geom_halfeyeh(size=5, .width=0.50, fill=NA, density.color=NA) + 
		geom_vline(xintercept=0, linetype='dashed') +
		labs(y = '', x ='') +
		scale_x_continuous(breaks=seq(-10, 6, 2), lim = c(-10, 6)) +
		scale_y_discrete( labels=ylabs, position='right') + 
		theme(axis.text.y=element_text(size=10, colour='black'),
				axis.text.x=element_text(size=10, colour='black'),
				panel.border = element_blank(), 
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
 				axis.line = element_line(colour = "black"))


## plot panels together
plot_grid(bs0,params0, labels=c('(a)', '(b)'), align='h',nrow=1, rel_widths=c(1.5,1), label_size=12)


dev.off()

## end of script

