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

# Define region variable which groups context-level units
L2.re <- "region"

# Define outcome variable
y <- "y"

# Define column in census containing bin size
bin.size <- "n"

# Define fraction of EBMA size (NULL defaults to 1/4)
ebma.size <- NULL

# Define number of folds
k.folds <- 5

# Unit to be used in sampling to create CV folds
cv.sampling <- "units"

# Set seed (NULL defaults to 12345)
seed <- NULL
