library(magrittr)

# Set seed
set.seed(1234)

# Load Stata data files
survey_data <- foreign::read.dta("./data/.tmp/survey_data.dta")
census_data <- foreign::read.dta("./data/.tmp/census_data.dta")

# Select relevant variables
survey_data <- survey_data %>%
  dplyr::select(y, age, educ, gXr, state, stateid, region,
                pvote, religcon, urban, unemp, hispanics, white)

census_data <- census_data %>%
  dplyr::select(age, educ, gXr, stateid, region,
                pvote, religcon, urban, unemp, hispanics, white,
                n)

# Draw a sample of size n_sample from full survey data set
# (We make sure to have at least n_min respondents per geographic unit)
n_sample <- 1500
n_min <- 6
select_rows <- sample(1:nrow(survey_data), size = n_sample, replace = FALSE)
survey_sample <- survey_data[select_rows, ]

while(min(table(survey_sample$stateid)) < n_min) {
  select_rows <- sample(1:nrow(survey_data), size = n_sample, replace = FALSE)
  survey_sample <- survey_data[select_rows, ]
}

# Save data files
save(survey_data, file = "./data/survey_data.RData")
save(census_data, file = "./data/census_data.RData")
save(survey_sample, file = "./data/survey_sample.RData")

# Remove data files from workspace
rm(survey_data, census_data, survey_sample,
   n_min, n_sample, select_rows)
