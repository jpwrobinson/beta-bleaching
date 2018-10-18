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

## ------------------------------------------------------------
##--------------- predicted beta change year on year ---------------------------------------------
# ---------------------------------------------------------------------------

## load posteriors
load(file='results/03_beta_94_model.Rdata')

pred.clean<-data.frame(beta.bray = c(pred.mean[c(1:12)],pred.mean[c(13:24)]),
						state = rep(c('Shifted', 'Recovering'), each=12),
						Year = rep(seq(min(betaS$Year),max(betaS$Year), length.out=12),times=2))

pred.clean$PI89.lower<-c(pred.PI89[1,c(1:12)], pred.PI89[1,c(13:24)])
pred.clean$PI89.upper<-c(pred.PI89[2,c(1:12)], pred.PI89[2,c(13:24)])

betaS$pt.col<-ifelse(betaS$state == 'Recovering', recovering, shifted)

bs1<-ggplot(pred.clean, aes(Year, beta.bray, fill=state, col=state)) + 
		geom_ribbon(alpha=0.1, show.legend=F,linetype=0,aes( ymin = PI89.lower, ymax=PI89.upper)) +
		geom_line(lwd=1, aes(linetype=state)) + 
		scale_fill_manual(values=c(recovering,shifted)) +
		scale_linetype_manual(values=c(1,5)) +
		scale_colour_manual(values=c(recovering,shifted)) +
		theme(axis.text.x=element_text(size=10, colour='black'),
			axis.text.y=element_text(size=10, colour='black'),
			axis.title.y=element_text(size=15, colour='black'),
				panel.border = element_blank(),
				 panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
		 		axis.line = element_line(colour = "black"),
				legend.position='none', 
				legend.title=element_blank(), 
				legend.text=element_text(size=8)) +
		guides(fill=guide_legend(nrow=1)) +
		scale_y_continuous(limits=c(0.25, 0.75)) +
		scale_x_continuous(breaks=c(2005, 2008, 2011, 2014, 2017), labels=c(2005, 2008, 2011, 2014, 2017), limits=c(2004.5, 2017.5)) +
		labs(x = '', y =expression(paste(beta['seq']))) +
		geom_jitter(data = betaS, aes(Year, beta.bray), width=0.25, size=1, alpha=0.3)


## now parameter estimates
ylabs<-c('Shifted regime', 'Year','Regime shift * Year', 'Protection', 
		'Branching coral', 'Massive coral', 'Encrusting coral', 'Macroalgae', 'Complexity')
post<-as.data.frame(extract.samples(m4)) %>% gather(param, dist) %>% 
		filter(param %in% c('ba', 'bb', 'bc', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk'))

params1<-ggplot(post, aes(x = dist, y = param)) + 
		geom_halfeyeh(size=0.5, .prob=0.95, fill=NA, density.color=NA) + 
		geom_halfeyeh(size=5, .prob=0.50, fill=NA, density.color=NA) + 
		geom_vline(xintercept=0, linetype='dashed') +
		labs(y = '', x ='Standardized effect size') +
		scale_x_continuous(breaks=seq(-0.10, 0.10, 0.05), lim = c(-0.10, 0.10)) +
		scale_y_discrete( labels=ylabs, position='right') + 
		theme(axis.text.y=element_text(size=10, colour='black'),
				axis.text.x=element_text(size=10, colour='black'),
				axis.title.x=element_text(size=12, colour='black'),
				panel.border = element_blank(), 
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
 				axis.line = element_line(colour = "black"))



## plot panels together
plot_grid(bs1,params1, labels=c('(c)', '(d)'), align='h',nrow=1, rel_widths=c(1.5,1), label_size=12)


dev.off()

## end of script
