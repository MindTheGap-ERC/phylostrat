#### Load packages and utility functions ####
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(admtools)
library(StratPal)
library(phyclust)

source("code/utils.r")

#### Fix seed ####
set.seed(42)

#### Load age-depth models from Hohmann et al (2024) ####
da = load("data/osf/R_outputs/ageDepthModelsScenariosAandB.Rdata")
adm_A_2km = admtools::tp_to_adm(t = ageDepthModels$A$`2 km`$time,
                                h = ageDepthModels$A$`2 km`$height,
                                T_unit = "Myr",
                                L_unit = "m")
# shorten ADM from scenario B to have the same duration
t = ageDepthModels$B$`6 km`$time  -  max(ageDepthModels$B$`6 km`$time) + max_time(adm_A_2km)
sel = t>=0
h = ageDepthModels$B$`6 km`$height[sel]
adm_B_6km = admtools::tp_to_adm(t = t[sel],
                                h = h - min(h),
                                T_unit = "Myr",
                                L_unit = "m")
#### Constants ####

t_max = max_time(adm_A_2km) # total duration of simulation
n_sim = 1 # number of trees simulated
lambda = 1 # origination rate
mu = 0.3 # extinction rate
rho = 0  # sampling fraction in the present (= top of section)
rate_sampling_true = 100 # rate of fossil recovery per lineage
rate_bin_evol = 1 # rate of evolution for the binary characters
n_char = 1000 # number of characters sampled per fossil
hiat_min = 0.5 # minimum duration of hiatus (Myr) to be considered in the skyline model
n_pres_fossils = 30 # number of preserved fossils after stratigraphic effects

#### simulate complete tree ####

# make sure the tree is not trivial (e.g., one tip or goes fully extinct before t_max)
repeat{
  tree_complete = TreeSim::sim.bd.age(age = t_max,
                                      numbsim = n_sim,
                                      lambda = lambda,
                                      mu = mu,
                                      mrca =FALSE,
                                      complete = TRUE)[[1]]
  if (class(tree_complete) == "phylo"){
    break
  }
}

ape::write.nexus(tree_complete,
                 file = "data/tree_complete.nex")
plot.phylo(tree_complete)
axis(1)

#### Simulate fossil record ####
# simulate complete fossil record along the tree
fossils_full = FossilSim::sim.fossils.poisson(rate = rate_sampling_true,
                                              tree = tree_complete)

#### Remove fossils coinciding with gaps ####
plot(x = adm_A_2km$t, 
     y = adm_A_2km$t |> get_collection_prob(adm_A_2km)(),
     xlab = "Time",
     ylab = "Collection probability")
# generate sampled ancestor tree

tree_sa_complete = FossilSim::SAtree.from.fossils(tree = tree_complete, 
                                                  fossils = fossils_full)
fossils_full = tree_sa_complete$fossils # get fossils with tip names
# apply preservation of fossils
pres = identity
ctc = get_collection_prob(adm_A_2km)
fossils_inc_A = tree_sa_complete$fossils |>
  admtools::rev_dir(ref = t_max) |> 
  StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
  admtools::rev_dir(ref = t_max) |>
  subsample_fossils(n = n_pres_fossils)
ctc = get_collection_prob(adm_B_6km)
fossils_inc_B = tree_sa_complete$fossils |>
  admtools::rev_dir(ref = t_max) |> 
  StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
  admtools::rev_dir(ref = t_max) |>
  subsample_fossils(n = n_pres_fossils)

fossils_cont = tree_sa_complete$fossils |> subsample_fossils(n = n_pres_fossils)

plot(fossils_full, tree = tree_complete, rho = rho)
plot(fossils_cont, tree = tree_complete, rho = rho)
plot(fossils_inc_A, tree = tree_complete, rho = rho)
plot(fossils_inc_B, tree = tree_complete, rho = rho)

#### simulate character matrix ####
# tree with only fossils
tree_sf_complete = FossilSim::sampled.tree.from.combined(tree = tree_sa_complete$tree,
                                                         rho = rho)
# simulate full character matrix
char_mat_full = sim_bin_char(tree = tree_sf_complete,
                        rate = rate_bin_evol, 
                        nchar = n_char)
## subset character matrix
# select specimens preserved from the adm
char_mat_inc_A = char_mat_full[fossils_inc_A$tip.label,]
char_mat_inc_B = char_mat_full[fossils_inc_B$tip.label,]
# select equal no. of specimens, continuous sampling
char_mat_cont = char_mat_full[fossils_cont$tip.label,]
# write character matrices
ape::write.nexus.data(x = char_mat_full, 
                      file = "data/char_mat_full.nex",
                      format = "standard")
ape::write.nexus.data(x = char_mat_inc, 
                      file = "data/char_mat_inc.nex",
                      format = "standard")
ape::write.nexus.data(x = char_mat_cont, 
                      file = "data/char_mat_cont.nex",
                      format = "standard")

#### Export specimen ages ####
df_full = fossils_full |>
  get_fossil_ages()

df_inc = fossils_inc |>
  get_fossil_ages()

df_cont = fossils_cont |>
  get_fossil_ages()

# export 
write.table(df_full, 
            file = "data/taxa_full.tsv",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)
write.table(df_inc, 
            file = "data/taxa_inc.tsv",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)
write.table(df_cont, 
            file = "data/taxa_cont.tsv",
            sep = "\t",
            quote = FALSE, 
            row.names = FALSE)

#### Determine time points where no preservation is happening  for skyline model ####
get_skyline_breakpoints = function(adm, hiat_min){
  ind = get_hiat_list(adm) |> 
    sapply( \(x){x["end"] - x["start"]})|>
    unname() |> 
    (\(x){which(x > hiat_min)})()
  skyline = c(min_time(adm), max_time(adm))
  for (i in ind){
    skyline = c(skyline, get_hiat_list(adm)[[i]][c("start", "end")])
  }
  return(skyline |> unname() |> sort())
}
skyline = get_skyline_breakpoints(adm, hiat_min = hiat_min)

skyline_bp_age = t_max - skyline |> rev()
