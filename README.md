# Script for JAGS analysis of Defossez et al. accepted Oikos

This repository list the files used in the analysis of ground vegetation and
adult canopy direct and indirect effect on tree seedlings survival and growth
along climatic gradient for five species in the French Alps. Variation of micro-environmental variables with ground vegetation and canopy cover is analysed. The
R scripts use [JAGS](http://mcmc-jags.sourceforge.net/) and
[R2jags](http://cran.r-project.org/web/packages/R2jags/index.html) and
also require [ggplot2](http://ggplot2.org/) for the figures. The R
script are run with the package
[remake](https://github.com/richfitz/remake) that facilitate
reproducible workflow. 

Data are available on
[Dryad package](http://datadryad.org/resource/doi:10.5061/dryad.2j5s7)

The code run with the package [remake]
**TODO INCLUDE DATA DOWLOAD FROM DRYAD IN R**

## Run analysis
The R code will download the data from dryad, run the jags analysis
and create the output tables and figures.

```r
remake::make('all')
```



## Folders
- **R** : R functions used in the analysis

- **download** : created folder with downloaded data from dryad.
- **output** : created folder of output of R 
- **figures** : created folder for pdf figures 
- **jags.model** : created folder with jags models 

