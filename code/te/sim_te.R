#### Load packages and utility functions ####
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(admtools)
library(StratPal)
library(phyclust)

#### Create directories ####
path = "data/sim/te/" # path to sim outputs
path_osf = "data/osf" # data to osf 
path_mcmc = "output"
if (!dir.exists(path)){
  dir.create(path)
}
if(!dir.exists(path_osf)){
  dir.create(path_osf)
}
if(!dir.exists(path_mcmc)){
  dir.create(path_mcmc)
}

#### Load age-depth models from Hohmann et al (2024) ####
if (!file.exists("data/osf/R_outputs/ageDepthModelsScenariosAandB.Rdata")){
  # download raw data if needed
  files = osfr::osf_retrieve_node("https://osf.io/zbpwa/") |>
    osfr::osf_ls_files() |>
    osfr::osf_download(path = "data/osf", conflicts = "overwrite")
}
da = load("data/osf/R_outputs/ageDepthModelsScenariosAandB.Rdata")
adm_A_2km = admtools::tp_to_adm(t = ageDepthModels$A$`2 km`$time,
                                h = ageDepthModels$A$`2 km`$height,
                                T_unit = "Myr",
                                L_unit = "m")
# shorten ADM from scenario B to have the same duration
t = ageDepthModels$B$`6 km`$time  -  max(ageDepthModels$B$`6 km`$time) + admtools::max_time(adm_A_2km)
sel = t>=0
h = ageDepthModels$B$`6 km`$height[sel]
adm_B_6km = admtools::tp_to_adm(t = t[sel],
                                h = h - min(h),
                                T_unit = "Myr",
                                L_unit = "m")
#### Constants ####
n_chars =   c(100, 300, 1000) # number of characters sampled per fossil
runs =  50 # no of runs 

t_max = admtools::max_time(adm_A_2km) # total duration of simulation
lambda = 1 # origination rate
mu = 0.3 # extinction rate
rate_sampling_true = 100 # rate of fossil recovery per lineage
hiat_min = 0.5 # minimum duration of hiatus (Myr) to be considered in the skyline model
n_pres_fossils = 30 # number of preserved fossils after stratigraphic effects
length_alignment = 2000
mod = "-mHKY -t5 -a0.25 -g5" # HKY model, 5 categories for the gamma distribution. alpha shape parameter set to 0.25, transition to transversion ratio is 5
morph_rate = 0.1 # clock rate of strict morphological clock
cases = c("inc_A", "inc_B", "cont")

# parameters for the strict molecular clock
#sd_mol = 0.1 # standard deviation (log scale) for the strict molecular clock
#obs_mean_mol = 0.1 # observed mean of the lognormal distribution for the strict molecular clock
clock_rate_mol = 0.1

#### set seed ####
seed = 4
set.seed(seed)

#### aux functions ####
get_fossil_ages = function(fossils, extant_taxa = NULL){
  #' @title get ages of samples from fossils object
  #' 
  #' @param fossils a fossil object
  #' @param recent_taxa list of taxon names to be added as recent (i.e., with age 0)
  if (is.null(extant_taxa)){
    df = data.frame(taxon = fossils$tip.label,
                    min_age = fossils$hmin,
                    max_age = fossils$hmax)
  }
  if (!is.null(extant_taxa)){
    n_recent = length(extant_taxa)
    df = data.frame(taxon = c(fossils$tip.label, extant_taxa),
                    min_age = c(fossils$hmin, rep(0, n_recent)),
                    max_age = c(fossils$hmax, rep(0, n_recent)))
  }
  
  return(df)
}

get_collection_prob = function(adm){
  #' @title determine collection probability
  #' 
  #' @param adm an age-depth modek
  #' 
  #' @description
  #' given an age-depth model, returns a function that returns 0 if the time coincides with a hiatus
  #' and 1 if not. Outside of the domain of the age-depth model, it returns 1 
  #' 
  collection_prob = function(x){
    admtools::is_destructive(adm, x, out_dom_mode = "conservative") |> (\(x){!x})() |> as.numeric()
  }
  return(collection_prob)
}

subsample_fossils = function(fossils, n){
  #' randomly select a number of fossils
  ind = sample(seq_len(length(fossils$sp)), size = n)
  return(fossils[ind,])
}

select_fossils = function(fossils, case, n){
  #' select fossils according to one of the sampling cases
  #'     pres = identity
  pres = identity
  if (case == "inc_A"){
    ctc = get_collection_prob(adm_A_2km)
    fossils_i= fossils |>
      admtools::rev_dir(ref = t_max) |> 
      StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
      admtools::rev_dir(ref = t_max) |>
      subsample_fossils(n = n_pres_fossils)
    return(fossils_i)
  }
  if (case == "inc_B"){
    ctc = get_collection_prob(adm_A_2km)
    fossils_i = fossils |>
      admtools::rev_dir(ref = t_max) |> 
      StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
      admtools::rev_dir(ref = t_max) |>
      subsample_fossils(n = n_pres_fossils)
    return(fossils_i)
  }
  if (case =="cont"){
    fossils_i = fossils |> subsample_fossils(n = n_pres_fossils)
    return(fossils_i)
  }
  stop("unknown sampling case")
}

sim_bin_char = function(tree, nchar){
  #' @title simulate discrete traits along tree
  #' 
  #' @param tree tree
  #' @param nchar number of characters to simulate
  #' 
  par = matrix(data = c(-1, 1, 1, -1), 
               nrow = 2, 
               ncol = 2)
  char_mat = geiger::sim.char(phy = tree,
                              par = par,
                              model = "discrete",
                              nsim = nchar)
  char_mat = apply(char_mat, 1, function(x) x -1) |> t()
  return(char_mat)
}

for (id in seq_len(runs)){
  cat(paste0(id, "\n"))
  for (n_char in n_chars){
    for (case in cases){
    #### simulate tree ####
    # make sure the tree is not trivial (e.g., one tip or goes fully extinct before t_max)
    repeat{
      tree_complete = TreeSim::sim.bd.age(age = t_max,
                                          numbsim = 1,
                                          lambda = lambda,
                                          mu = mu,
                                          mrca =FALSE,
                                          complete = TRUE)[[1]]
      if (class(tree_complete) == "phylo"){
        break
      }
    }
    
    fossils_full = FossilSim::sim.fossils.poisson(rate = rate_sampling_true,
                                                  tree = tree_complete)
    
    fossils = select_fossils(fossils_full, case = case, n = n_pres_fossils)
    
    tree_sa = FossilSim::SAtree.from.fossils(tree = tree_complete, 
                                             fossils = fossils)
    
    tree_rho1 = FossilSim::sampled.tree.from.combined(tree = tree_sa$tree,
                                                      rho = 1)
    tree_rho0 = FossilSim::sampled.tree.from.combined(tree = tree_sa$tree,
                                                      rho = 0)
    extant_tips = setdiff(tree_rho1$tip.label, tree_rho0$tip.label )
    
    fossil_file_name = paste0(path, "fossils_", id, "_nchar", n_char, "_", case, "_rho1.csv")
    tree_sa$fossils |>
      get_fossil_ages(extant_taxa = extant_tips)|>
      write.table(file = fossil_file_name,
                  quote = FALSE,
                  row.names = FALSE)
    
    tree_file_name = paste0(path, "tree_", id, "_nchar", n_char, "_", case, "_rho1.nex")
    ape::write.nexus(tree_rho1,
                     file = tree_file_name)
    
    tree_w_rate_morph = tree_rho1
    tree_w_rate_morph$edge.length = tree_w_rate_morph$edge.length * morph_rate
    char_mat = sim_bin_char(tree = tree_w_rate_morph,
                            nchar = n_char)
    char_mat_name = paste0(path, "char_mat_", id, "_nchar", n_char, "_", case, "_rho1.nex")
    ape::write.nexus.data(char_mat,
                          file = char_mat_name,
                          format = "standard")
    
    tree_w_rate_mol = tree_rho1
    tree_w_rate_mol$edge.length = tree_w_rate_mol$edge.length * clock_rate_mol
    opts = paste0(mod, " -l", length_alignment)
    a = phyclust::seqgen(opts = opts, rooted.tree = tree_w_rate_mol) |>
      strsplit(" +")
    l = list()
    for (i in 1:length(a)){
      if (a[[i]][1] %in% extant_tips){
        l[[ a[[i]][1] ]] = unlist(strsplit(a[[i]][2], split = ""))
      }
    }
    mol_dat_name = paste0(path, "mol_dat_", id, "_nchar", n_char, "_", case, "_rho1.nex")
    ape::write.nexus.data(l,
                          file = mol_dat_name)
    }
  }
}

#### Determine time points where no preservation is happening  for skyline model ####
# get_skyline_breakpoints = function(adm, hiat_min){
#   ind = get_hiat_list(adm) |> 
#     sapply( \(x){x["end"] - x["start"]})|>
#     unname() |> 
#     (\(x){which(x > hiat_min)})()
#   skyline = c(min_time(adm), max_time(adm))
#   for (i in ind){
#     skyline = c(skyline, get_hiat_list(adm)[[i]][c("start", "end")])
#   }
#   return(skyline |> unname() |> sort())
# }
# skyline = get_skyline_breakpoints(adm_A_2km, hiat_min = hiat_min)
# 
# skyline_bp_age = t_max - skyline |> rev()


cat("Done! \n")
