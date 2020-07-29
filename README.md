# autoMrP
autoMrP improves the prediction performance of multilevel regression with post-stratification (MrP) by combining a number of machine learning methods through ensemble Bayesian model averaging (EBMA). The full discussion can be found here: Broniecki, Leemann, and Wüest. 2020. "Improving Multilevel Regression with Post-Stratification Through Machine Learning (autoMrP)." Journal of Politics, conditionally accepted.

## Installation
To install autoMrP from GitHub run:

```R
devtools::install_github("retowuest/autoMrP")
```

## Examples

autoMrP can be used to improve forecasts by combining predictions from various classifiers and it can be used to estimate the standard MrP model. We run through both applications.

For this example, we use the data that is included in the package (see `?autoMrP::survey` and `?autoMrP::census`).

```R
census <- autoMrP::census
survey <- autoMrP::survey_item
```

### The Standard MrP Model

We include the three individual level variables L1x1, L1x2 and L1x3 as well as the context-level variables L2.x1 and L2.x2. In addition, we include context-level random effects for states nested in regions.

Note, that we set the argument `L2.x.scale = FALSE` which stops `autoMrP()` from normalizing the context-level variables. In our example data, the context-level variables are already normalized. Context-level variables must be normalized because some of our classifiers are not scale invariant.

```R
mrp_model <- autoMrP::auto_MrP(
  y = "YES",
  L1.x = c("L1x1", "L1x2", "L1x3"),
  L2.x = c("L2.x1", "L2.x2"),
  L2.unit = "state",
  L2.reg = "region",
  L2.x.scale = FALSE,
  survey = survey,
  census = census,
  bin.proportion = "proportion",
  best.subset = FALSE,
  lasso = FALSE,
  pca = FALSE,
  gb = FALSE,
  svm = FALSE,
  mrp = TRUE
)
```


### Better Predictions Through Machine Learning

In this example, we accept the default settings for the tuning parameters and combine predictions from four classifiers (best subset, lasso, pca, gb, svm) via EBMA to one overall prediction. The classifiers make use of all six context-level variables in the data. In addition, we include the standard MrP model with the context-level variables that we used in the previous example in EBMA. To exclude the standard MrP model from EBMA, we would need to set the argument `mrp = FALSE`. Running this example will take some time (~20 minutes).

Note, that because we want the standard MrP model to use only two context-level variables while we want our classifiers to use all six context-level variables, we specify the context-level variables to be used for the standard MrP model in the argument `mrp.L2.x = c("L2.x1", "L2.x2")`.

EBMA tuning is quite exhaustive with the default settings and takes quite long. We reduce the search grid to speed up the proccess by setting the argument `ebma.tol` and `ebma.n.draws`.

Finally, we set the argument `gb.L2.unit = FALSE` because, in our applications, the gradient boosting classifier tended to perform worse with many state context-level indicators (there are 48 states in our example data).

```R
out <- autoMrP::auto_MrP(
  y = "YES",
  L1.x = c("L1x1", "L1x2", "L1x3"),
  L2.x = c("L2.x1", "L2.x2", "L2.x3", "L2.x4", "L2.x5", "L2.x6"),
  L2.unit = "state",
  L2.reg = "region",
  L2.x.scale = FALSE,
  survey = survey,
  census = census,
  bin.proportion = "proportion",
  best.subset = TRUE,
  lasso = TRUE,
  pca = TRUE,
  gb = TRUE,
  svm = TRUE,
  mrp = TRUE,
  mrp.L2.x = c("L2.x1", "L2.x2"),
  ebma.tol = c(0.001, 0.0005),
  ebma.n.draws = 10,
  gb.L2.unit = FALSE
)
```

During parameter tuning, models may not converge for certain parameter values and cause warning messages if `verbose = TRUE`.
