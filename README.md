# Script for JAGS analysis of Defossez et al. accepted Oikos

This repository list the files used in the analysis of ground vegetation and
adult canopy direct and indirect effect on tree seedlings survival and growth
along climatic gradient for five species in the French Alps. Variation of micro-environmental variables with ground vegetation and canopy cover is analysed. The
R scripts use [JAGS](http://mcmc-jags.sourceforge.net/) and
[R2jags](http://cran.r-project.org/web/packages/R2jags/index.html) and
also reauire [ggplot2](http://ggplot2.org/) for the figures.

Data are available on Dryad DOI: doi:10.5061/dryad.2j5s7 (provisional)

**TODO INCLUDE DATA DOWLOAD FROM DRYAD IN R**

## Scripts
- seedling survival analysis: survival.analysis.2014.R
- seedling growth analysis: growth.analysis.2014.R
- PCA of climate from one weather station per site: climate.pca.R
- Canopy effect on micro-environmental variables: Canop_clim.R
- Herb effect on micro-environmental variables: Herb_clim.R


## Folders

- **data** : raw data used in the growth and survival analysis and
  associated metadata
- **R** : R functions used in the analysis
- **output** : output of R (dir created by the R script)
- **figs** : pdf figure (dir created by R scripts)
- **jags.output** : jags files (dir created by R scripts)

