# load requried packages
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(phyclust)
library(admtools)
library(StratPal)

#### Constants ####
n_chars = c(30, 100, 300, 1000) # number of characters sampled per fossil
t_max = 5 # duration of simulation in Myr
lambda = 0.8 # origination rate
mu = 0.6 # extinction rate
clock_rate_morph = 0.05 # morph. clock rate (strict clock)
clock_rate_mol = 0.005 # mol clock rate (strict clock)
sampling_rate = 500 # total sampling rate for fossils
n_fossils = 30 # no of preserved fossils
n_rep = 50 # number of replicates
length_alignment = 2000 # lenhth of molecular alignment
mod = "-mHKY -t5 -a0.25 -g5" # HKY model, 5 categories for the gamma distribution. alpha shape parameter set to 0.25, transition to transversion ratio is 5
ids = seq_len(n_rep) # ids for the replicates
opts = paste0(mod, " -l", length_alignment)