# beta-bleaching
R code and data accompanying **Climate change induces long-term persistence of novel coral reef fish assemblages**.  *In review.*

The following R packages were used to analyse data and create figures.

```
## Data analysis
install.packages(c("tidyverse", "rethinking", "devtools", "betapart", "here")
library(devtools)
install_github('jpwrobinson/funk')

## Figures
install.packages(c("tidyverse", "rethinking", "XXXYYY")
```

#### R scripts for statistical models 

Scripts load diversity estimates and fit Bayesian models underlying main results.



- Community richness model in [1_richness.R](/scripts/1_richness.R)
- Functional group richness change models in [2_FG_richness.R](scripts/2_FG_richness.R)
- beta_94 models in [3_beta_94.R](scripts/3_beta_94.R)
- beta_seq models in [4_beta_seq.R](scripts/4_beta_seq.R)
- Model of species-level biomass change from 1994 to 2014 in [5_biom_diff.R](scripts/5_biom_diff.R)
- Null model for reshuffling species among sites in [6_null_model.R](scripts/6_null_model.R)

#### **Model predictions and Figures**

Scripts load Bayesian models, generate predictions from posterior distributions, then create figures.



- Community richness in [Figure1a,b]('figures/richness_fig.R')

- beta_94 in [Figure1c,d]('figures/beta_94_fig.R')

- beta_seq in [Figure1e,f]('figures/beta_seq_fig.R')

- Biomass change in [Figure4b,c]('figures/biomass_fig.R')
