
library(R.utils)
library(here)

## navigate to beta-bleaching
setwd(here('beta-bleaching'))

## recreate figures
dir.create('figures-pdf')
sourceDirectory("scripts/figures/", print.eval = TRUE)

## rerun statistical analyses
## WARNING! running all bayesian models can take several hours/days on a desktop computer
sourceDirectory("scripts/analysis/")