#### Load packages and utility functions ####
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(admtools)
library(StratPal)

source("code/utils.r")

#### Fix seed ####
set.seed(42)

#### Constants ####
dist = "2km" # distance from shore of adm
adm = tp_to_adm(t = scenarioA$t_myr,      # age-depth model
                h = scenarioA$h_m[,dist],
                L_unit = "m",
                T_unit = "Myr")
t_max = max_time(adm) # total duration of simulation
n_sim = 1 # number of trees simulated
lambda = 1 # origination rate
mu = 0.3 # extinction rate
rho = 0  # sampling fraction in the present (= top of section)
rate_sampling_true = 10 # rate of fossil recovery per lineage
rate_bin_evol = 1 # rate of evolution for the binary characters
n_char = 1000 # number of characters sampled per fossil
hiat_min = 0.5 # minimum duration of hiatus (Myr) to be considered in the skyline model

#### simulate complete tree ####
tree_complete = TreeSim::sim.bd.age(age = t_max,
                                    numbsim = n_sim,
                                    lambda = lambda,
                                    mu = mu,
                                    mrca =TRUE,
                                    complete = TRUE)[[1]]
ape::write.nexus(tree_complete,
                 file = "data/tree_complete.nex")
plot.phylo(tree_complete)
axis(1)

#### Simulate fossil record ####
# simulate complete fossil record along the tree
fossils_full = FossilSim::sim.fossils.poisson(rate = rate_sampling_true,
                                              tree = tree_complete)

#### Remove fossils coinciding with gaps ####
plot(x = scenarioA$t_myr, 
     y = scenarioA$t_myr |> get_collection_prob(adm)(),
     xlab = "Time",
     ylab = "Collection probability")
# generate sampled ancestor tree
tree_sa_complete = FossilSim::SAtree.from.fossils(tree = tree_complete, 
                                                  fossils = fossils_full)
# apply preservation of fossils
pres = identity
ctc = get_collection_prob(adm)
fossils_inc = tree_sa_complete$fossils |>
  admtools::rev_dir(ref = t_max) |> 
  StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
  admtools::rev_dir(ref = t_max)

n_fossils = FossilSim::count.fossils(fossils_inc)
fossils_cont = tree_sa_complete$fossils |> subsample_fossils(n = n_fossils)

plot(fossils_full, tree = tree_complete, rho = rho)
plot(fossils_cont, tree = tree_complete, rho = rho)
plot(fossils_inc, tree = tree_complete, rho = rho)


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
char_mat_inc = char_mat_full[fossils_inc$tip.label,]
# select equal no. of specimens, continuous sampling
char_mat_cont = char_mat_full[fossils_cont$tip.label,]
# write character matrices
ape::write.nexus.data(x = char_mat, 
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
  FossilSim::SAtree.from.fossils(tree_complete, fossils = _) |>
  (\(x){x$tree})() |>
  FossilSim::sampled.tree.from.combined(rho = rho) |>
  get_fossil_ages(t_max = t_max)

df_inc = fossils_inc |>
  FossilSim::SAtree.from.fossils(tree_complete, fossils = _) |>
  (\(x){x$tree})() |>
  FossilSim::sampled.tree.from.combined(rho = rho) |>
  get_fossil_ages(t_max = t_max)

df_cont = fossils_cont |>
  FossilSim::SAtree.from.fossils(tree_complete, fossils = _) |>
  (\(x){x$tree})() |>
  FossilSim::sampled.tree.from.combined(rho = rho) |>
  get_fossil_ages(t_max = t_max)

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
get_skyline_breakpoints(adm, hiat_min = hiat_min)

