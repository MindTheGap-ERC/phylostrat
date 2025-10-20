#### Import constants & libraries ####
source("code/constants.R")

path = "data/sim/fbd_gap_est/" # path to store outputs
# create output directory if it does not exist
if (!dir.exists(path)){
  dir.create(path, recursive = TRUE)
}

set.seed(45)

#### Aux functions ####
source("code/utils.R")

#### Stratigraphic context ####
sinusoid_data = read.csv("data/strat/sinusoid_adm.csv")
adm_sinusoid = tp_to_adm(t = sinusoid_data$time..Myr., h = sinusoid_data$adm_1..m.)

source("code/case_sinusoid.R")
