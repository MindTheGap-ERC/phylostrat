# simulation for a baseline - no model violation

library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(admtools)
library(StratPal)
library(phyclust)

#### Constants ####
n_chars = c(30, 300, 1000) # number of characters sampled per fossil
runs = 10 # no of runs 
path = "data/baseline/" # path to store outputs
t_max = 2 #admtools::max_time(adm_A_2km) # total duration of simulation
n_sim = 1 # number of trees simulated
lambda = 1 # origination rate
mu = 0.3 # extinction rate
rate_sampling_true = 3 # rate of fossil recovery per lineage
rate_bin_evol = 1 # rate of evolution for the binary characters
hiat_min = 0.5 # minimum duration of hiatus (Myr) to be considered in the skyline model
n_pres_fossils = 30 # number of preserved fossils after stratigraphic effects
length_alignment = 2000
mod = "-mHKY -t5 -a0.25 -g5" # HKY model, 5 categories for the gamma distribution. alpha shape parameter set to 0.25, transition to transversion ratio is 5
# parameters for the strict morphological clock
sd_morph = 0.1 # standard deviation (log scale) for strict morph. clock
obs_mean_morph = 0.1 # observed mean of the lognormal distribution for strict morph. clock

# parameters for the strict molecular clock
sd_mol = 0.1 # standard deviation (log scale) for the strict molecular clock
obs_mean_mol = 0.1 # observed mean of the lognormal distribution for the strict molecular clock 

#### set seed ####
seed = 4
set.seed(seed)

if (!dir.exists("data/baseline")){
  dir.create("data/baseline")
}




for (id in seq_len(runs)){
  
  repeat{
    tree_complete = TreeSim::sim.bd.age(age = t_max,
                                        numbsim = n_sim,
                                        lambda = lambda,
                                        mu = mu,
                                        mrca =FALSE,
                                        complete = TRUE)[[1]]
    if (class(tree_complete) == "phylo"){
      h = tree_complete$root.edge + node.depth.edgelength(tree_complete)
      n_recent = tree_complete$tip.label[h > t_max - 100 * .Machine$double.eps] |> length()
      break
    }
  }
  
  
  ape::write.nexus(tree_complete,
                   file = paste0(path, "tree_complete_",id,".nex"))
  plot.phylo(tree_complete, main = as.character(id))
  axis(1)
  
  fossils_full = FossilSim::sim.fossils.poisson(rate = rate_sampling_true,
                                                tree = tree_complete)
  
  plot(fossils_full, tree = tree_complete, rho = 0)
  
  
  # generate sampled ancestor tree
  
  tree_sa_complete = FossilSim::SAtree.from.fossils(tree = tree_complete, 
                                                    fossils = fossils_full)
  plot(tree_sa_complete$tree)
  axis(1)
  ## tree with fossils also sampled at the present
  tree_fp_complete = FossilSim::sampled.tree.from.combined(tree = tree_sa_complete$tree,
                                                           rho = 1)
  fossils_full = tree_sa_complete$fossils # get fossils with tip names
  
  h = tree_fp_complete$root.edge + ape::node.depth.edgelength(tree_fp_complete)
  recent_tips = tree_fp_complete$tip.label[h > t_max - 100 * .Machine$double.eps]
  
  # save fossils
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
  
  
  # save all fossils
  fossils_full |>
    get_fossil_ages()|>
    write.table(file = paste0(path,"fossils_full_rho0_",id,".tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  
  fossils_full |>
    get_fossil_ages(extant_taxa = recent_tips)|>
    write.table(file = paste0(path,"fossils_full_rho1_",id,".tsv"),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  
  
  # simulate character data
  #### Simulate character matrices ####
  # strict morphological clock with clock rate following a lognormal distribution with mean 1
  mu_clock_morph = log(obs_mean_morph) - sd_morph^2/2
  clock_rate_morph = rlnorm(1, meanlog = mu_clock_morph, sdlog = sd_morph)
  
  tree_w_rate_morph = tree_fp_complete
  tree_w_rate_morph$edge.length = tree_w_rate_morph$edge.length * clock_rate_morph
  sim_bin_char = function(tree, rate, nchar){
    #' @title simulate discrete traits along tree
    #' 
    #' @param tree tree
    #' @param rate rate for Q matrix
    #' @param nchar number of characters to simulate
    #' 
    par = matrix(data = c(-rate, rate, rate, -rate), 
                 nrow = 2, 
                 ncol = 2)
    char_mat = geiger::sim.char(phy = tree,
                                par = par,
                                model = "discrete",
                                nsim = nchar)
    char_mat = apply(char_mat, 1, function(x) x -1) |> t()
    return(char_mat)
  }
  
  # simulate full morph character matrix
  char_mat_full = sim_bin_char(tree = tree_w_rate_morph,
                               rate = rate_bin_evol, 
                               nchar = max(n_chars))
  
  # export character data
  for (n_char in n_chars){
    ## subset character matrix
    # select specimens preserved from the adm and correct no of characters
    char_mat_rho0 = char_mat_full[fossils_full$tip.label,seq_len(n_char)] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_rho0_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    char_mat_rho1 = char_mat_full[c(fossils_full$tip.label, recent_tips),seq_len(n_char)] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_rho1_nchar", n_char, "_",id,".nex"),
                            format = "standard")
  }
  
  
  # simulate mol data
  #### simulate molecular data ####
  mu_clock_mol = log(obs_mean_mol) - sd_mol^2/2
  clock_rate_mol = rlnorm(1, meanlog = mu_clock_mol, sdlog = sd_mol)
  
  tree_w_rate_mol = tree_fp_complete
  tree_w_rate_mol$edge.length = tree_w_rate_mol$edge.length * clock_rate_mol
  
  opts = paste0(mod, " -l", length_alignment)
  a = phyclust::seqgen(opts = opts, rooted.tree = tree_w_rate_mol) |>
    strsplit(" +")
  l = list()
  for (i in 1:length(a)){
    if (a[[i]][1] %in% recent_tips){
      l[[ a[[i]][1] ]] = unlist(strsplit(a[[i]][2], split = ""))
    }
  }
  write.nexus.data(l,
                   file = paste0(path, "mol_data_", id, ".nex"))

}

