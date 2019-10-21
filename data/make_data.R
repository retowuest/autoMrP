library(magrittr)

# Load dta data files
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

# Save data files
save(survey_data, file = "./data/survey_data.RData")
save(census_data, file = "./data/census_data.RData")

# Remove data files from workspace
rm(survey_data, census_data)
