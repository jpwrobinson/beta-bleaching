# beta-bleaching
R code and data accompanying **Climate change induces long-term persistence of novel coral reef fish assemblages**.  *In review.*

The following R packages were used to analyse data and create figures.

```
install.packages(c("tidyverse", "rethinking", "devtools", "betapart", "here", "vegan")
library(devtools)
install_github('jpwrobinson/funk')
```

****

#### R scripts for statistical analyses 

Scripts contains analyses underlying main results, including spatial and null diversity functions, Bayesian model structures, and SIMPER analysis:

- Community richness model in [1_richness.R](/scripts/1_richness.R)
- Functional group richness change models in [2_FG_richness.R](scripts/2_FG_richness.R)
- beta_94 models in [3_beta_94.R](scripts/3_beta_94.R)
- beta_seq models in [4_beta_seq.R](scripts/4_beta_seq.R)
- Function for pairwise beta spatial estimates in [5_beta_spatial.R](scripts/5_beta_spatial.R)
- SIMPER dissimilarity analysis in [6_simper.R](scripts/6_simper_spatial.R)
- Model of species-level biomass change from 1994 to 2014 in [7_biom_diff.R](scripts/7_biom_diff.R)
- Null model for reshuffling species among sites in [8_null_model.R](scripts/8_null_model.R)

****

#### Figures

Scripts load Bayesian models, generate predictions from posterior distributions, then create figures.

- Community richness in [Figure1a,b](figures/Figure1_ab.R)
- beta_94 in [Figure1c,d](figures/Figure1_cd.R)
- beta_seq in [Figure1e,f](figures/Figure1_ef.R)
- Functional group richness in [Figure 2](figures/Figure2.R)
- Spatial beta diversity in [Figure 3][scripts/Figure3.R]
- SIMPER analysis in [Figure 4a](figures/Figure4_a.R)
- Modelled biomass change in [Figure4b,c]('figures/Figure4_bc.R')
- Observed biomass change in [Figure 5]('figures/Figure5.R')

****
#### Datasets and model predictions

*Datasets*

* Community biomass matrices for 1994 and 2014 in [SEY_biomass_matrix.Rdata](data/SEY_biomass_matrix.Rdata)
* Beta_1994 estimates in [UVC_beta_1994.csv](data/UVC_beta_1994.csv)
* Beta_seq estimates in [UVC_beta_seq.csv](data/UVC_beta_seq.csv)
* Spatial beta diversity estimates in [UVC_beta_spatial.csv](data/UVC_beta_spatial.csv)
* Species richness estimates in [UVC_richness.csv](data/UVC_richness.csv)
* Functional group richness estimates in [UVC_richnessdiff_FG.csv](data/UVC_richnessdiff_FG.csv)
* Species-level biomass change from 1994 to 2014 in [UVC_biom_change.csv](data/UVC_biom_change.csv)

*Model predictions*

* Species richness model model (Figure 1a,b): [01_richness_model.Rdata](results/01_richness_model.Rdata)
* Corallivore richness model (Figure 2b): [02_richness_corallivore.Rdata](results/02_richness_corallivore.Rdata)
* Herbivore richness model (Figure 2c): [02_richness_herbivore.Rdata](results/02_richness_herbivore.Rdata)
* Invertivore richness model (Figure 2d): [02_richness_invertivore.Rdata](results/02_richness_invertivore.Rdata)
* Mixed-diet richness model (Figure 2e)[02_richness_mixed-diet.Rresults/data](02_richness_mixed-diet.Rdata)
* Planktivore richness model (Figure 2a): [02_richness_piscivore.Rdata](results/02_richness_piscivore.Rdata)
* Piscivore richness model (Figure 2f): [02_richness_planktivore.Rdata](results/02_richness_planktivore.Rdata)
* Beta_1994 estimates (Figure 1c,d): [03_beta_94_model.Rdata](results/03_beta_94_model.Rdata)
* Beta_seq estimates (Figure 1e,f): [04_beta_seq_model.Rdata](results/04_beta_seq_model.Rdata)
* SIMPER analysis (Figure 4a): [06_simper.Rdata](results/06_simper.Rdata)
* Biomass change from 1994 to 2014 (Figure 4b,c): [07_biom_change_model.Rdata](results/07_biom_change_model.Rdata)

****
Data and scripts are licensed under Creative Commons Attribution 4.0 International.
