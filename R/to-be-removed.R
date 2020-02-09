# This file contains code used for testing
# File should be removed after finishing package

# Load packages
library(magrittr)

# Clean working environment
rm(list = ls())

# Load survey data
load(here::here("data", "survey_sample.RData"))

# Load census data
load(here::here("data", "census_data.RData"))

# Assign data sets to objects with "correct" names
survey <- survey_sample
census <- census_data

# Define individual-level covariates
L1.x <- c("age", "educ", "gXr")

# Define context-level covariates
L2.x <- c("pvote", "religcon", "urban", "unemp", "hispanics", "white")

# Define context-level unit
L2.unit <- "stateid"

# Define region in which context-level units are nested
L2.reg <- "region"

# Define outcome variable
y <- "y"

# Define column in census containing bin size
bin.size <- "n"

# Define fraction of EBMA size (NULL defaults to 1/3)
ebma.size <- NULL

# Define number of folds
k.folds <- 5

# Unit to be used in sampling to create CV folds
cv.sampling <- "L2 units"

# Set seed (NULL defaults to 12345)
seed <- NULL

# Set verbose option
verbose <- TRUE

# Define unit for loss function
loss.unit <- "individual"

# Define measure for loss function
loss.measure <- "mse"

# Define lambdas as vector for Lasso
lasso.lambda.set <- c(1, 2, 3)

# Define lambdas as data frame for Lasso
lasso.lambda.set <- data.frame(step_size = c(0.1, 1),
                               threshold = c(1, 5))

# Define maximum number of iterations w/o improvement for Lasso
lasso.iterations.max <- NULL

# Define set of interaction depths for GB
gb.interaction.set <- c(1, 2, 3)

# Define set of learning rates for GB
gb.shrinkage.set <- c(0.04, 0.01)

# Define initial number of total trees
gb.tree.start <- 50

# Define increase in number of trees as scalar for GB
gb.tree.increase.set <- 50

# Define increase in number of trees as vector for GB
gb.tree.increase.set <- c(50, 100)

# Define maximum number of trees as scalar for GB
gb.trees.max.set <- 200

# Define maximum number of trees as vector for GB
gb.trees.max.set <- c(200, 400)

# Define maximum number of iterations w/o improvement for GB
gb.iterations.max <- NULL

# Define minimum number of observations in terminal nodes for GB
gb.n.minobsinnode <- 5

# Define whether L2.unit should be inlcuded in GB
gb.L2.unit.include <- FALSE

# Define whether L2.reg should be included in GB
gb.L2.reg.include <- FALSE

# Define kernel for SVM
svm.kernel <- "radial"

# Define error function for SVM
svm.error.fun <- "MSE"   # // might need to be changed to NULL

# Define gamma parameters for SVM
svm.gamma.set <- c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4)

# Define cost parameters for SVM
svm.cost.set <- c(1, 10)

# Define number of draws from EBMA sample with equal obs/state
Ndraws <- 100

# Define tolerance values for EBMA tuning
tol.values <- c(0.01, 0.005, 0.001)

#_______________________________________________________
# load functions

# aux functions
source("./utils.R")

# classifiers
source("./best_subset.R")
source("./best_subset_classifier.R")
source("./lasso.R")
source("./lasso_classifier.R")
source("./pca.R")
source("./gb.R")
source("./gb_classifier.R")
source("./post_stratification.R")
source("./ebma.R")
source("./svm_classifier.R")
source("./svm.R")


