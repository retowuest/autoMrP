# autoMrP
autoMrP improves the prediction performance of multilevel regression with post-stratification (MrP) by combining a number of machine learning methods through ensemble bayesian model averaging (EBMA).

## Installation
Step 1) autoMrP depends on EBMAforecast which is currently only available on the cran archieve. To install EBMAforecast run the following code.

```R
packageurl <- "https://cran.r-project.org/src/contrib/Archive/EBMAforecast/EBMAforecast_0.52.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

Note that EBMAforecast depends on the stats, graphics, separationplot, plyr, methods, Hmisc and abind packages.

Step 2) To install autoMrP from GitHub run

```R
devtools::install_github("retowuest/autoMrP")
```