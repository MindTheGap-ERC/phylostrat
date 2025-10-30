#### Import constants & libraries ####
source("code/constants.R")

path = "sim_data/fbd_gap_prior/" # path to store outputs
# create output directory if it does not exist
if (!dir.exists(path)){
  dir.create(path, recursive = TRUE)
}

set.seed(44)

#### Aux functions ####
source("code/utils.R")

#### Stratigraphic context ####
source("code/strat_scenarios.R")

source("code/case_sinusoid.R")
