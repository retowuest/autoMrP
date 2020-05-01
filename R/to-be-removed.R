# This file contains code used for testing
# File should be removed after finishing package

# Load packages
library(magrittr)

# Clean working environment
rm(list = ls())

# Load survey data
load(here::here("data", "survey_paper_i11.RData"))

# Load census data
load(here::here("data", "census_paper_i11.RData"))

# Assign data sets to objects with "correct" names
#survey <- survey_sample
#census <- census_data

# Set seed (NULL defaults to 12345)
seed <- NULL

# Define outcome variable
#y <- "y"
y <- "YES"

# Define individual-level covariates
#L1.x <- c("age", "educ", "gXr")
L1.x <- c("L1x1", "L1x2", "L1x3")

# Define context-level covariates
#L2.x <- c("pvote", "religcon", "urban", "unemp", "hispanics", "white")
L2.x <- c("L2.x1", "L2.x2", "L2.x3", "L2.x4", "L2.x5", "L2.x6")

# Define context-level unit
#L2.unit <- "stateid"
L2.unit <- "L2.unit"

# Define region in which context-level units are nested
L2.reg <- "region"

# Whether to scale context level variables (defaults to TRUE)
L2.x.scale <- FALSE

# provide user controlled principal components
pcs <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

# provide user controlled fold assignment
folds <- "fold"

# Define column in census containing proportion
bin.proportion <- "proportion"

# Define column in census containing bin size
#bin.size <- "n"
bin.size <- NULL

# Define fraction of EBMA size (defaults to 1/3)
ebma.size <- 1/3

# Define number of folds
k.folds <- 5

# Unit to be used in sampling to create CV folds
cv.sampling <- "L2 units"

# Define unit for loss function
loss.unit <- "individual"

# Define measure for loss function
loss.fun <- "MSE"

# Whether to include best subset classifier
best.subset <- TRUE

# Whether to include lasso classifier
lasso <- TRUE

# Whether to include pca classifier
pca <- TRUE

# Whether to include gb classifier
gb <- TRUE

# Whether to include svm classifier
svm <- TRUE

# Whether to include mrp classifier
mrp <- FALSE

# Whether to include forward selection
forward.select <- FALSE

# Define best subset context-level covariates
best.subset.L2.x <- NULL

# Define lasso context-level covariates
lasso.L2.x <- NULL

# Define gb context-level covariates
gb.L2.x <- NULL

# Define svm context-level covariates
svm.L2.x <- NULL

# Define mrp context-level covariates
mrp.L2.x <- NULL

# Define whether L2.unit should be inlcuded in GB
gb.L2.unit <- FALSE

# Define whether L2.reg should be included in GB
gb.L2.reg <- TRUE

# Define lambdas as vector for Lasso
lasso.lambda.set <- c(1, 2)

# Define lambdas as data frame for Lasso
#lasso.lambda.set <- data.frame(step_size = c(2, 5),
#                               threshold = c(1, 5))

# Define maximum number of iterations w/o improvement for Lasso
lasso.iterations.max <- NULL

# Define set of interaction depths for GB
gb.interaction.set <- c(1, 2, 3)

# Define set of learning rates for GB
gb.shrinkage.set <- c(0.04, 0.01, 0.008, 0.005, 0.001)

# Define initial number of total trees
gb.tree.start <- 50

# Define increase in number of trees as scalar for GB
gb.tree.increase.set <- 50

# Define increase in number of trees as vector for GB
#gb.tree.increase.set <- c(50, 100)

# Define maximum number of trees as scalar for GB
gb.trees.max.set <- 1000

# Define maximum number of trees as vector for GB
#gb.trees.max.set <- c(200, 400)

# Define maximum number of iterations w/o improvement for GB
gb.iterations.max <- 70

# Define minimum number of observations in terminal nodes for GB
gb.n.minobsinnode <- 5

# Define kernel for SVM
svm.kernel <- "radial"

# Define error function for SVM
svm.error.fun <- "MSE"   # // might need to be changed to NULL

# Define gamma parameters for SVM
svm.gamma.set <- c(0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1, 2, 3, 4)

# Define cost parameters for SVM
svm.cost.set <- c(1, 10)

# Define context-level covariates for standard MRP
#mrp.L2.x <- c("pvote", "religcon", "urban", "unemp", "white")
mrp.L2.x <- c("L2.x1", "L2.x2", "L2.x3", "L2.x4", "L2.x6")

# EBMA draws
ebma.n.draws <- 100

# tolerance values for ebma tuning
ebma.tol.values <- c(0.01, 0.001)

# Whether to obtain bootstrapped prediction uncertainty
uncertainty <- FALSE

# Set verbose option
verbose <- TRUE


#_______________________________________________________
# load functions

# aux functions
source(here::here("R", "utils.R"))

# classifiers
source(here::here("R", "run_best_subset.R"))
source(here::here("R", "best_subset_classifier.R"))
source(here::here("R", "run_lasso.R"))
source(here::here("R", "lasso_classifier.R"))
source(here::here("R", "run_pca.R"))
source(here::here("R", "run_gb.R"))
source(here::here("R", "gb_classifier.R"))
source(here::here("R", "svm_classifier.R"))
source(here::here("R", "run_svm.R"))
source(here::here("R", "post_stratification.R"))
source(here::here("R", "ebma.R"))


