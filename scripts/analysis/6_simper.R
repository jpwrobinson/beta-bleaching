#!/bin/env Rscript

library(vegan)
library(dplyr)

## load dataframe fish.mat.full (site-level community biomass matrices)
load('data/SEY_biomass_matrix.Rdata') 

# create labels by year + state
lab<-data.frame(Location=rownames(fish.mat.full))
lab$year<-c(rep(c(1994, 2014), each=21))
lab$state<-fish$state[match(lab$Location, fish$Location)]

## assign simper groups for comparisons
lab$simper.group<-ifelse(lab$label=='1994.Shifted', 'Pre.Shifted', 'Post.Shifted')
lab$simper.group<-ifelse(lab$label=='1994.Recovering', 'Pre.Recovering', lab$simper.group)
lab$simper.group<-ifelse(lab$state=='Recovering' & lab$year != 1994, 'Post.Recovering', lab$simper.group)

## ----------------------------------------- ###
### GET SIMPER SCORES FROM BIOMASS MATRICES ###
## ----------------------------------------- ###

## Simper analysis
simp<-simper(fish.mat.full, group=lab$simper.group)

## extract relevant state comparisons
simp2<-data.frame(species=simp[[2]]$species,contribution=simp[[2]]$average, sd = simp[[2]]$sd, ava=simp[[2]]$ava, avb=simp[[2]]$avb,comparison='Pre.Shifted.Post.Shifted')
simp5<-data.frame(species=simp[[5]]$species,contribution=simp[[5]]$average, sd = simp[[5]]$sd, ava=simp[[5]]$ava, avb=simp[[5]]$avb,comparison='Pre.Recovering.Post.Recovering')
simp.scores<-rbind(simp2, simp5)

## add functional group to simp.scores
simp.scores$FG<-fish$FG.coarse[match(simp.scores$species, fish$Species)]
## add biomass change as plus or minus
simp.scores$abund.diff<-with(simp.scores, avb-ava)
simp.scores$change <- ifelse(simp.scores$abund.diff>0, 'Gain', 'Loss')

## sum by functional feeding group (FG)
simp.FG.biom<-simp.scores %>% group_by(comparison, change,FG) %>% summarise(score=sum(contribution))
simp.FG.biom$contribution.change<-ifelse(simp.FG.biom$change=='Loss', -simp.FG.biom$score, simp.FG.biom$score)



## end save scores and lab plotting info 
simp.scores.biom<-simp.scores
save(simp.FG.biom, simp.scores.biom, file = 'results/06_simper.Rdata')