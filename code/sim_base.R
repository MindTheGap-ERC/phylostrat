# load requried packages
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(phyclust)

# set seed
set.seed(42)

#### Constants ####
n_chars = c(30, 100, 300, 1000) # number of characters sampled per fossil
t_max = 5 # duration of simulation in Myr
lambda = 0.8 # origination rate
mu = 0.6 # extinction rate
clock_rate_morph = 0.05 # morph. clock rate (strict clock)
clock_rate_mol = 0.005 # mol clock rate (strict clock)
sampling_rate = 50 # total sampling rate for fossils
n_fossils = 30 # no of preserved fossils
path = "data/sim/fbd_base/" # path to store outputs
n_rep = 50 # number of replicates
length_alignment = 2000 # lenhth of molecular alignment
mod = "-mHKY -t5 -a0.25 -g5" # HKY model, 5 categories for the gamma distribution. alpha shape parameter set to 0.25, transition to transversion ratio is 5
ids = seq_len(n_rep) # ids for the replicates

# create output directory if it does not exist
if (!dir.exists(path)){
  dir.create(path, recursive = TRUE)
}

#### Auxiliary functions ####
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

subsample_fossils = function(fossils, n){
  #' randomly select a number of fossils
  ind = sample(seq_len(length(fossils$sp)), size = n)
  return(fossils[ind,])
}

sim_bin_char = function(tree, nchar){
  #' @title simulate discrete traits along tree
  #' 
  #' @param tree tree
  #' @param rate rate for Q matrix
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

sim_mol_char = function(tree, extant_tips, opts){
  #' @title simulate molecular partition for extant tips
  #' 
  a = phyclust::seqgen(opts = opts, rooted.tree = tree_w_rate_mol) |>
    strsplit(" +")
  l = list()
  for (i in 1:length(a)){
    if (a[[i]][1] %in% extant_tips){
      l[[ a[[i]][1] ]] = unlist(strsplit(a[[i]][2], split = ""))
    }
  }
  return(l)
}


n_succ = 0
n_fail = 0
for (id in ids){
  cat(paste0("Running simulation no ", id, "\n"))
  for (n_char in n_chars){
    # simulate tree, making sure it makes it to the present
    repeat{
      tree_complete = TreeSim::sim.bd.age(age = t_max,
                                          numbsim = 1,
                                          lambda = lambda,
                                          mu = mu,
                                          mrca =FALSE,
                                          complete = TRUE)[[1]]
      if (class(tree_complete) == "phylo"){ # check if multiple tips are in the present, see ?sim.bd.age
        n_succ = n_succ +1
        break 
      }
      if (class(tree_complete) != "phylo"){
        n_fail = n_fail + 1
      }
    }
    
    # full fossil record
    fossils_full_base = FossilSim::sim.fossils.poisson(rate = sampling_rate,
                                                       tree = tree_complete) |> 
      subsample_fossils(n = n_fossils)
    # turn into sampled ancestor tree
    tree_sa_complete = FossilSim::SAtree.from.fossils(tree = tree_complete,
                                                      fossils = fossils_full_base)
    # fossils with tip label names
    fossils_full = tree_sa_complete$fossils
    
    tree_fp_rho1 = FossilSim::sampled.tree.from.combined(tree = tree_sa_complete$tree,
                                                         rho = 1)
    tree_fp_rho0 = FossilSim::sampled.tree.from.combined(tree = tree_sa_complete$tree,
                                                         rho = 0)
    extant_tips = setdiff(tree_fp_rho1$tip.label, tree_fp_rho0$tip.label )
    
    if(length(extant_tips) == 0){
      stop("No extant tips in tree")
    }
    
    # save full fossil record
    fossil_name = paste0(path,"fossils_", id , "_nchar", n_char, ".csv")
    fossils_full |>
      get_fossil_ages(extant_taxa = extant_tips)|>
      write.table(file = fossil_name,
                  quote = FALSE,
                  row.names = FALSE)
    
    tree_name = paste0(path,"tree_", id , "_nchar", n_char, ".nex")
    ape::write.nexus(tree_fp_rho1,
                     file= tree_name)
    
    tree_w_rate_morph = tree_fp_rho1
    tree_w_rate_morph$edge.length = tree_w_rate_morph$edge.length * clock_rate_morph
    char_mat = sim_bin_char(tree = tree_w_rate_morph,
                            nchar = n_char)
    char_mat_name = paste0(path,"char_mat_", id , "_nchar", n_char, ".nex")
    ape::write.nexus.data(char_mat,
                          file = char_mat_name,
                          format = "standard")
    #plot(tree_fp_rho1)
    
    ## molecular alignment
    opts = paste0(mod, " -l", length_alignment)
    
    tree_w_rate_mol = tree_fp_rho1
    tree_w_rate_mol$edge.lenth = tree_w_rate_mol$edge_length * clock_rate_mol
    mol = sim_mol_char(tree = tree_w_rate_mol, extant_tips = extant_tips, opts = opts)
    ape::write.nexus.data(mol,
                          file = paste0(path, "mol_dat_", id, "_nchar", n_char, ".nex"))
  }
}

cat(paste0(n_succ, " successful trees\n"))
cat(paste0(n_fail, " failed trees\n"))

cat("Done \n")