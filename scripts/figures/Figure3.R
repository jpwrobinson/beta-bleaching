
print('Creating Figure 3')

library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)


############ Figure 3 - spatial beta patterns ####################

## plot details
recovering<-c('#045a8d')	
shifted<-c('#b30000')
theme_set(theme_bw())

## read data
beta<-read.csv(file='data/UVC_beta_spatial.csv')

## baseline: mean 1994 estimates
beta.base<-beta %>% 
            filter(Year==1994) %>% 
            group_by(Year, type) %>% 
            summarise(beta.bray = mean(beta.bray))

## 2005-2017 beta, all estimates for boxplots
beta.plot<-beta %>% 
            filter(Year!=1994) %>% 
            mutate(state=factor(type))


## labels
grob1 <- grobTree(textGrob("Recovering baseline", x=0.01,  y=0.28, hjust=0,
  gp=gpar(col="#636363", fontsize=9)))
grob2 <- grobTree(textGrob("Shifted baseline", x=0.01,  y=0.3025, hjust=0,
  gp=gpar(col="#636363", fontsize=9)))

## create plots

## recovering
bs0<-ggplot(beta.plot[beta.plot$state=='Recovering',], aes(factor(Year), beta.bray, fill=state, col=state)) + 
		geom_hline(yintercept=beta.base$beta.bray[beta.base$type=='Recovering'], 
			linetype='longdash', col='grey') +
		scale_fill_manual(values=c(recovering)) +
		scale_colour_manual(values=c(recovering)) +
		geom_boxplot(show.legend=F, aes(factor(Year), beta.bray, fill=state)) +
		stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 1, col='white', show.legend=F) +
		theme(axis.text.x=element_text(size=12, colour='black'),
			axis.text.y=element_text(size=12, colour='black'),
			axis.title.y=element_text(size=16, colour='black'),
			axis.title.x=element_text(size=0, colour='black'),
				panel.border = element_blank(),
				 panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
		 		axis.line = element_line(colour = "black"),
				legend.text=element_text(size=8)) +
		scale_y_continuous(limits=c(0.44,0.71)) +
		labs(x = '', y =expression(beta['spatial'])) +
		annotation_custom(grob1)

## shifted
bs1<-ggplot(beta.plot[beta.plot$state=='Shifted',], aes(factor(Year), beta.bray, fill=state, col=state)) + 
		geom_hline(yintercept=beta.base$beta.bray[beta.base$type=='Shifted'], 
			linetype='longdash', col='grey') +
		scale_fill_manual(values=c(shifted)) +
		scale_colour_manual(values=c(shifted)) +
		geom_boxplot(show.legend=F, aes(factor(Year), beta.bray, fill=state)) +
		stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 1, col='white', show.legend=F) +
		theme(axis.text.x=element_text(size=12, colour='black'),
			axis.text.y=element_text(size=12, colour='black'),
			axis.title.y=element_text(size=0, colour='black'),
			axis.title.x=element_text(size=0, colour='black'),
				panel.border = element_blank(),
				 panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
		 		axis.line = element_line(colour = "black"),
				legend.text=element_text(size=8)) +
		scale_y_continuous(limits=c(0.44,0.71)) +
		labs(x = '', y ='') +
		annotation_custom(grob2) 

pdf(file='figures-pdf/Figure3.pdf', height=7, width=9)
plot_grid(bs0, bs1, nrow=1, labels=c('(a)', '(b)'), label_size=14)
dev.off()
