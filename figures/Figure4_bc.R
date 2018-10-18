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

pdf(file='figures/Figure_4bc.pdf', height=4.5, width=15)
par(mfrow=c(1,2))

############################################
####### BIOMASS BAYESIAN POSTERIORS ###########
############################################

## prepare colours 
fg.cols<-data.frame(col=c('#984ea3','#e41a1c','#4daf4a','#a65628','#ff7f00','#377eb8'),
		FG=c('Planktivore','Corallivore', 'Herbivore',  'Invertivore','Mixed-diet Feeder','Piscivore'))
cols<-fg.cols$col; names(cols)<-fg.cols$FG

load('results/07_biom_change_model.Rdata')
biom.full<-biom20


# names for ranef levels
ranef1<-c("Planktivore", 'Corallivore', 'Herbivore', 'Invertivore', 'MixeddietFeeder', 'Piscivore')
ranef2<-rep(ranef1, each=2)
ranef2[seq(1, 11, 2)]<-paste0('Rec.', ranef2[seq(1,11,2)])
ranef2[seq(2, 12, 2)]<-paste0('Shif.', ranef2[seq(2,12,2)])
fg<-levels(biomS$FG)

## now parameter estimates

post<-as.data.frame(extract.samples(m)) %>% gather(param, dist) %>% 
			filter(param %in% c(
					"a_fg.6","a_fg.1","a_fg.2","a_fg.3","a_fg.4","a_fg.5", ## FG intercepts = recovering biomass
					"ba_fg.6","ba_fg.1","ba_fg.2","ba_fg.3","ba_fg.4","ba_fg.5", ## FG slopes = shifted biomass
					'bl', 'ba', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk', 'a')) ## fixed parameters


cols.reduced<-c(rep(c('#984ea3','#e41a1c','#4daf4a','#a65628','#ff7f00','#377eb8'), each=2), rep('#969696',8))
names(cols.reduced)<-c(ranef2,'bl', 'ba', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk', 'a')

ylabs<-c("Planktivore",'Corallivore', 'Herbivore', 'Invertivore', 'Mixed-diet\n Feeder', 'Piscivore', 'Body size',
	'Regime shift','Protection', 'Branching coral','Massive coral', 'Encrusting coral', 'Macroalgae', 'Complexity')

post$param<-str_replace_all(post$param, 'ba_fg', 'Shif.')
post$param<-str_replace_all(post$param, 'a_fg', 'Rec.')
post$param<-str_replace_all(post$param, '.1', fg[1])
post$param<-str_replace_all(post$param, '.2', fg[2])
post$param<-str_replace_all(post$param, '.3', fg[3])
post$param<-str_replace_all(post$param, '.4', fg[4])
post$param<-str_replace_all(post$param, '.5', fg[5])
post$param<-str_replace_all(post$param, '.6', fg[6])

post$state<-str_split_fixed(post$param, '\\.', 2)[,1]


post$groups<-post$param
post$groups<-str_replace_all(post$groups, 'Shif.', '')
post$groups<-str_replace_all(post$groups, 'Rec.', '')

post$x.num<-as.numeric(factor(post$groups, levels=c(ranef1,'bl', 'ba', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk')))
post$x.num[grepl('Rec.', post$param)]<-post$x.num[grepl('Rec.', post$param)]-0.15
post$x.num[grepl('Shif.', post$param)]<-post$x.num[grepl('Shif.', post$param)]+0.15
post$alpha<-ifelse(grepl('Shif.',post$param), 0.5, 1)
post$alpha[post$param=='bA']<-1

cols.reduced2<-c(c('#984ea3','#e41a1c','#4daf4a','#a65628','#ff7f00','#377eb8'),rep('black',8))
names(cols.reduced2)<-c(ranef1,'bl', 'ba', 'bf', 'bg', 'bh', 'bi', 'bj', 'bk')

## divide into ranef and fixef
post.fg<-post[post$x.num < 7,]
post.fix<-post[which(post$x.num >= 7),]


params.fg<-ggplot(post.fg, aes(x = dist, y = x.num, alpha=alpha, col=state)) + 
		geom_vline(xintercept=0, linetype='dashed', size=0.5, col='black') +
		geom_halfeyeh(size=5, trim=TRUE,  fill=NA,alpha=0.8, .prob=0.5) + 
		# geom_density_ridges(size=5, trim=TRUE, fill='grey90', alpha=0.8, .prob=0.5) + 
		geom_halfeyeh(size=1, trim=TRUE, fill=NA, alpha=1, .prob=0.95) + 
		labs(y = '', x ='Standardized effect size') +
		scale_x_continuous(breaks=seq(-15, 15, 5), lim = c(-15, 15)) +
		scale_y_continuous(lim = c(0.5, 6.5), breaks=1:6, labels=ylabs[1:6], position='right') + 
		guides(fill=FALSE, alpha=FALSE, col=FALSE) +
		scale_alpha(range = c(1)) +
		# scale_fill_manual(name='state', values = c('azure4', 'grey90')) +
		scale_color_manual(name='state', values = tot.cols) +
		theme(axis.text.y=element_text(size=12, colour='black'),
				axis.text.x=element_text(size=12, colour='black'),
				axis.title.x=element_text(size=13, colour='black'),
				axis.title.y=element_text(size=13, colour='black'),
				panel.border = element_blank(), 
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
 				axis.line = element_line(colour = "black"))

params.fix<-ggplot(post.fix, aes(x = dist, y = x.num, alpha=alpha)) + 
		geom_vline(xintercept=0, linetype='dashed', size=0.5, col='black') +
		geom_halfeyeh(size=1, fill=NA, alpha=0, .prob=0.95) + 
		geom_halfeyeh(size=5, fill=NA, alpha=0, .prob=0.5) + 
		labs(y = '', x ='Standardized effect size') +
		scale_x_continuous(breaks=seq(-7, 7, 2), lim = c(-7, 7)) +
		scale_y_continuous(lim = c(7, 14), breaks=7:14, labels=ylabs[7:14], position='right') + 
		guides(fill=FALSE, alpha=FALSE) +
		scale_alpha(range = c(1)) +
		# scale_fill_manual(name='param', values = cols.reduced) +
		# scale_color_manual(name='state', values = tot.cols) +
		theme(axis.text.y=element_text(size=12, colour='black'),
				axis.text.x=element_text(size=12, colour='black'),
				axis.title.x=element_text(size=13, colour='black'),
				axis.title.y=element_text(size=13, colour='black'),
				panel.border = element_blank(), 
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
 				axis.line = element_line(colour = "black"))


## print ggs with viewport
vp <- viewport(height = unit(1,"npc"), width=unit(0.33, "npc"), 
              just = c("centre"),
              y = 0.5, x = 0.52)
print(params.fg, vp = vp)

vp <- viewport(height = unit(1,"npc"), width=unit(0.33, "npc"), 
              just = c("right"),
              y = 0.5, x = 1)
print(params.fix, vp = vp)

## add labels
par(xpd=NA)
add_label(1.15, 0.01, label='(b)', font=2 ,cex=1.6)
add_label(2.15, 0.01, label='(c)', font=2 ,cex=1.6)


dev.off()

## end of script

