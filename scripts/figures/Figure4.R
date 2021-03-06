
print('Creating Figure 5')

library(ggplot2) 
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(funk)

## data load
biom20<-read.csv(file='data/UVC_biom_change.csv')
biom20$FG <- factor(biom20$FG, levels = levels(biom20$FG)[c(6,1,2,3,4,5)])

## Get mean and SE
biom20.mean<-aggregate(biom ~ FG + Species + state, biom20, mean)
biom20.mean$biom.SE<-aggregate(biom ~ FG + Species + state, biom20, se)[,4]

biom20.mean$biom94<-aggregate(biom94 ~ FG + Species + state, biom20, mean)[,4]
biom20.mean$biom94.SE<-aggregate(biom94 ~ FG + Species + state, biom20, se)[,4]

biom20.mean$biom.change.raw<-with(biom20.mean, biom94 - biom)

## change to log biom
biom20.mean$biom<-log10(biom20.mean$biom +1)
biom20.mean$biom94<-log10(biom20.mean$biom94 +1)


## simplify naming
biom20<-biom20.mean

## add some colours
fg.cols<-data.frame(col=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628'),
    col2=c('#fc8d59', '#a6bddb', '#a1d99b', '#bcbddc', '#fed976', '#d8b365'),
    FG=c('Corallivore', 'Piscivore', 'Herbivore','Planktivore', 'Mixed-diet Feeder', 'Invertivore'))
biom20$col<-as.character(fg.cols$col[match(biom20$FG, fg.cols$FG)])
biom20$col.border<-'white'


## add absolute biomass difference, logged
biom20$absdiff<-log10(abs(biom20$biom.change.raw) + 1)
biom20$absdiff<-ifelse(biom20$biom.change.raw > 0 , -biom20$absdiff, biom20$absdiff)



## separate into recovering and shifted
t.rec<-biom20 %>% filter(state == 'Recovering')
t.rec<-t.rec %>% group_by(Species) %>% arrange(FG, desc(biom.change.raw))

t.shif<-biom20 %>% filter(state == 'Shifted')
t.shif<-t.shif[match(t.rec$Species,t.shif$Species),]



pdf(file='figures-pdf/Figure5.pdf', height=7, width=8.5)

## plot details
cex.ax = 1
cx.name = 0.45
big.lab.cx = 1



split.screen(rbind(c(0,0.47,0, 1), c(0.53, 1, 0, 1)))
screen(1)
par(mar=c(0,0,0,0), xpd=T, oma=c(3.5,1,1.5,1),tcl = -0.25,mgp = c(2, 0.6, 0))

## recovering reefs
segs<-barplot(t.rec$biom, horiz=T, yaxs='i',space=0,xlim=c(-2, 2),axisnames=F, axes=F,las=2, cex.names=0.5, plot=F)
barplot(t.rec$absdiff, col=as.character(t.rec$col), border=t.rec$col.border, horiz=T, yaxs='i',space=0, names.arg=t.rec$Species,
  xlim=c(-1.7, 2.30103),axisnames=F, axes=F,tck=-0.01,las=2, cex.names=cx.name)

### add log axes labels
xBig = log10(c(seq(1, 10, 1), seq(20, 100, 10),  200))
axis(1, at=c(0, 1, log10(20), log10(50), log10(100),  log10(200)), labels=c(1, 10, 20, 50, 100,  200), cex.axis=cex.ax, tcl=0)
axis(1, xBig , labels=rep("", length(xBig)), tcl=-0.4)
axis(1, -xBig , labels=rep("", length(xBig)), tcl=-0.4)
axis(1, at=-c(1, log10(20), log10(50), log10(100),  log10(200)), labels=c(-10, -20, -50, -100, -200), cex.axis=cex.ax, tcl=0)

par(xpd=F)
add_label(0, 0.018, '(a)', cex =1.4, font = 2, pos =4)


par(xpd=NA)
add_label(0.5, -0.05, 'Recovering', cex =big.lab.cx, font = 2, pos =1)

### legend names
add_label(1.05, 0.05, fg.cols$FG[2], col=as.character(fg.cols$col[2]), cex=1, font=2, pos =3)
add_label(1.05, 0.225, 'Mixed-diet feeder', col=as.character(fg.cols$col[5]), cex=1, font=2, pos =3)
add_label(1.05, 0.44, fg.cols$FG[6], col=as.character(fg.cols$col[6]), cex=1, font=2, pos =3)
add_label(1.05, 0.7, fg.cols$FG[3], col=as.character(fg.cols$col[3]), cex=1, font=2, pos =3)
add_label(1.05, 0.94, fg.cols$FG[1], col=as.character(fg.cols$col[1]), cex=1, font=2, pos =3)
add_label(1.05, 0.99, fg.cols$FG[4], col=as.character(fg.cols$col[4]), cex=1, font=2, pos =3)


screen(2)
par(mar=c(0,0,0,0), xpd=T, oma=c(3.5,1,1.5,1),tcl = -0.25,mgp = c(2, 0.6, 0))

## shifted reefs
segs<-barplot(t.shif$biom, horiz=T, yaxs='i',space=0,xlim=c(-2, 2),axisnames=F, axes=F,las=2, cex.names=0.5, plot=F)
barplot(t.shif$absdiff, col=as.character(t.shif$col), border=t.shif$col.border, horiz=T, yaxs='i',space=0, names.arg=t.shif$Species,
  xlim=c(-1.7, 2.30103),axisnames=F, axes=F,tck=-0.01,las=2, cex.names=cx.name)


### add log axes labels
xBig = log10(c(seq(1, 10, 1), seq(20, 100, 10),  200))
axis(1, at=c(0, 1, log10(20), log10(50), log10(100),  log10(200)), labels=c(1, 10, 20, 50, 100,  200), cex.axis=cex.ax, tcl=0)
axis(1, xBig , labels=rep("", length(xBig)), tcl=-0.4)
axis(1, -xBig[-20] , labels=rep("", length(xBig)-1), tcl=-0.4)
axis(1, at=-c(1, log10(20), log10(50), log10(100)), labels=c(-10, -20, -50, -100), cex.axis=cex.ax, tcl=0)


par(xpd=NA)
abline(h=6.1, lty=2, lwd=0.5)
abline(h=14.1, lty=2, lwd=0.5)
abline(h=51, lty=2, lwd=0.5)
abline(h=87.1, lty=2, lwd=0.5)
abline(h=118, lty=2, lwd=0.5)
add_label(0.94, 0.018, '(b)', cex =1.4, font = 2, pos =4)

par(xpd=NA)
add_label(0.5, -0.05, 'Regime shifted', cex =big.lab.cx, font = 2, pos =1)



mtext(1, text=expression(paste('biomass + 1 (kg ha'^'-1', ')', sep='')), line=2, outer=T, cex=big.lab.cx)

dev.off()

