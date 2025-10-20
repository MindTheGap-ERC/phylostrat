#### Import constants & libraries ####
source("code/constants.R")

path = "data/sim/fbd_strat/" # path to store outputs
# create output directory if it does not exist
if (!dir.exists(path)){
  dir.create(path, recursive = TRUE)
}

set.seed(43)

#### Aux functions ####
source("code/utils.R")

source("code/case_miller.R")
source("code/case_sinusoid.R")