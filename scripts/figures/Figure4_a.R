
############################################
####### SIMPER ANALYSIS ####################
############################################

pdf(file='figures-pdf/Figure4_a.pdf', height=6, width=9)


## data load
load(file='results/06_simper.Rdata')

## set plot details
par(mfrow=c(1,1))
r.col<-c('#045a8d')	
shif.col<-c('#b30000')
tot.cols<-c(r.col, shif.col)
barW=0.19 ## set bar width on plot

s3<-simp.FG.biom %>% select(-score) %>% 
				group_by(comparison) %>% 
				complete(FG,change, fill=list(contribution.change=0))

## add levels to FG
s3$FG<-factor(s3$FG, levels=c('Planktivore','Corallivore', 'Herbivore',  'Invertivore','Mixed-diet Feeder','Piscivore'))

## assign plotting info
s3$cols<-fg.cols$cols[match(s3$FG, fg.cols$FG)]
s3$cols2<-ifelse(s3$comparison=='Pre.Shifted.Post.Shifted', shif.col, r.col)
s3$xlab<-ifelse(s3$comparison=='Pre.Shifted.Post.Shifted', as.numeric(s3$FG)+0.2, as.numeric(s3$FG)-0.2)

## get mean changes
means<-s3 %>% group_by(comparison, FG) %>% summarise(mean = mean(contribution.change))


par(mar=c(3.4,4,1,0))
plot(0, 0, xlim=c(0.7, 6.2), col='transparent', ylab='', xlab='', axes=F,ylim=c(-15,20))

with(s3[s3$change=='Gain',],
rect(xleft=xlab-barW, ybottom=0, xright=xlab+barW, ytop=contribution.change*100, 
	col=as.character(cols2), border=NA))
with(s3[s3$change=='Loss',],
rect(xleft=xlab-barW, ybottom=0, xright=xlab+barW, ytop=contribution.change*100, 
	col=alpha(as.character(cols2), 0.5), border=NA))
axis(2, cex.axis = 1.4); mtext(2, text='Relative contribution (%)',line=2.5, cex=1.4)
par(xpd=F); abline(a = 0, b=0, lty=2)
par(xpd=T); add_label(0.01, 0.01, label='(a)', font=2 ,cex=1.6)
legend('topright', legend=c('Recovering', 'Shifted'), inset=c(0, 0.1), 
			fill=c(r.col, shif.col), cex=1.6, bty='n',border=NA, pt.cex=2)
axis(1, at=seq(1, 6, by=1), padj=0.5, cex.axis=1.1, 
	labels=c('Planktivore','Corallivore', 'Herbivore', 'Invertivore','Mixed-diet\n feeder','Piscivore'),las=1)

## add mean changes as white lines
a<-100*means$mean[means$FG=='Planktivore' & means$comparison=='Pre.Recovering.Post.Recovering']
segments(0.59, a, 1, a, col='white')
a<-100*means$mean[means$FG=='Planktivore' & means$comparison=='Pre.Shifted.Post.Shifted']
segments(1, a, 1.41, a, col='white')
a<-100*means$mean[means$FG=='Corallivore' & means$comparison=='Pre.Recovering.Post.Recovering']
segments(1.59, a, 2, a, col='white')
a<-100*means$mean[means$FG=='Corallivore' & means$comparison=='Pre.Shifted.Post.Shifted']
segments(2, a, 2.41, a, col='white')
a<-100*means$mean[means$FG=='Herbivore' & means$comparison=='Pre.Recovering.Post.Recovering']
segments(2.59, a, 3, a, col='white')
a<-100*means$mean[means$FG=='Herbivore' & means$comparison=='Pre.Shifted.Post.Shifted']
segments(3, a, 3.41, a, col='white')
a<-100*means$mean[means$FG=='Invertivore' & means$comparison=='Pre.Recovering.Post.Recovering']
segments(3.59, a, 4, a, col='white')
a<-100*means$mean[means$FG=='Invertivore' & means$comparison=='Pre.Shifted.Post.Shifted']
segments(4, a, 4.41, a, col='white')
a<-100*means$mean[means$FG=='Mixed-diet Feeder' & means$comparison=='Pre.Recovering.Post.Recovering']
segments(4.59, a, 5, a, col='white')
a<-100*means$mean[means$FG=='Mixed-diet Feeder' & means$comparison=='Pre.Shifted.Post.Shifted']
segments(5, a, 5.41, a, col='white')
a<-100*means$mean[means$FG=='Piscivore' & means$comparison=='Pre.Recovering.Post.Recovering']
segments(5.59, a, 6, a, col='white')
a<-100*means$mean[means$FG=='Piscivore' & means$comparison=='Pre.Shifted.Post.Shifted']
segments(6, a, 6.41, a, col='white')

dev.off()
