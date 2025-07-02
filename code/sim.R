#### Load packages and utility functions ####
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
library(admtools)
library(StratPal)
library(phyclust)

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
n_chars = c(30, 300, 1000) # number of characters sampled per fossil
runs = 1 # no of runs 
path = "data/sim/" # path to store outputs
t_max = admtools::max_time(adm_A_2km) # total duration of simulation
n_sim = 50 # number of trees simulated
lambda = 1 # origination rate
mu = 0.3 # extinction rate
rate_sampling_true = 100 # rate of fossil recovery per lineage
rate_bin_evol = 1 # rate of evolution for the binary characters
hiat_min = 0.5 # minimum duration of hiatus (Myr) to be considered in the skyline model
n_pres_fossils = 30 # number of preserved fossils after stratigraphic effects
length_alignment = 2000
mod = "-mHKY -t5 -a0.25 -g5" # HKY model, 5 categories for the gamma distribution. alpha shape parameter set to 0.25, transition to transversion ratio is 5

#### set seed ####
seed = 1234
set.seed(seed)

for (id in seq_len(runs)){
    #### simulate tree ####
    # make sure the tree is not trivial (e.g., one tip or goes fully extinct before t_max)
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
    plot.phylo(tree_complete)
    axis(1)
    
    #### Simulate fossil record ####
    # simulate complete fossil record along the tree
    fossils_full = FossilSim::sim.fossils.poisson(rate = rate_sampling_true,
                                                  tree = tree_complete)
    
    #### Remove fossils coinciding with gaps ####
    
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
    plot(x = adm_A_2km$t, 
         y = adm_A_2km$t |> get_collection_prob(adm_A_2km)(),
         xlab = "Time",
         ylab = "Collection probability")
    # generate sampled ancestor tree
    
    tree_sa_complete = FossilSim::SAtree.from.fossils(tree = tree_complete, 
                                                      fossils = fossils_full)
    
    ## tree with fossils also sampled at the present
    tree_fp_complete = FossilSim::sampled.tree.from.combined(tree = tree_sa_complete$tree,
                                                             rho = 1)
    # save tree with all fossils as sanity check
    ape::write.nexus(tree_fp_complete,
                     file = paste0(path, "tree_all_fossils_",id,".nex"))
    # identify tips that are recent
    h = tree_fp_complete$root.edge + ape::node.depth.edgelength(tree_fp_complete)
    recent_tips = tree_fp_complete$tip.label[h > t_max - 100 * .Machine$double.eps]
    # check if no of recent tips is correct - to catch numeric inconsistencies ins node.depth.edgelength
    if (length(recent_tips) != n_recent){
      stop("error in no of recent tips")
    }
    
    fossils_full = tree_sa_complete$fossils # get fossils with tip names
    
    
    subsample_fossils = function(fossils, n){
      #' randomly select a number of fossils
      ind = sample(seq_len(length(fossils$sp)), size = n)
      return(fossils[ind,])
    }
    # sampled according to scenario A
    pres = identity
    ctc = get_collection_prob(adm_A_2km)
    fossils_inc_A = tree_sa_complete$fossils |>
      admtools::rev_dir(ref = t_max) |> 
      StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
      admtools::rev_dir(ref = t_max) |>
      subsample_fossils(n = n_pres_fossils)
    # sampled according to scenario B
    ctc = get_collection_prob(adm_B_6km)
    fossils_inc_B = tree_sa_complete$fossils |>
      admtools::rev_dir(ref = t_max) |> 
      StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
      admtools::rev_dir(ref = t_max) |>
      subsample_fossils(n = n_pres_fossils)
    # continuoiusly sampled fossils
    fossils_cont = tree_sa_complete$fossils |> subsample_fossils(n = n_pres_fossils)
    
    plot(fossils_full, tree = tree_complete, rho = 0)
    plot(fossils_cont, tree = tree_complete, rho = 0)
    plot(fossils_inc_A, tree = tree_complete, rho = 0)
    plot(fossils_inc_B, tree = tree_complete, rho = 0)
    
    # export all fossil ages
    
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
    
    fossils_full |>
      get_fossil_ages()|>
      write.table(file = paste0(path,"fossils_full_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    # continuously sampled fossils, with rho  = 0
    fossils_cont |>
      get_fossil_ages()|>
      write.table(file = paste0(path,"fossils_cont_rho0_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    # continuously sampled fossils with rho = 1
    fossils_cont |>
      get_fossil_ages(extant_taxa = recent_tips)|>
      write.table(file = paste0(path,"fossils_cont_rho1_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    # fossil ages according to scenario A with rho = 0
    fossils_inc_A |>
      get_fossil_ages()|>
      write.table(file = paste0(path,"fossils_inc_A_rho0_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    # fossil ages according to scenario A with rho = 1
    fossils_inc_A |>
      get_fossil_ages(extant_taxa = recent_tips)|>
      write.table(file = paste0(path,"fossils_inc_A_rho1_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    # fossil ages according to scenario B with rho = 0
    fossils_inc_B |>
      get_fossil_ages()|>
      write.table(file = paste0(path,"fossils_inc_B_rho0_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    # fossil ages according to scenarion B with rho = 1
    fossils_inc_B |>
      get_fossil_ages(extant_taxa = recent_tips)|>
      write.table(file = paste0(path,"fossils_inc_B_rho1_",id,".tsv"),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE)
    
    
    #### simulate character matrix ####
    for (n_char in n_chars){
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
                                  nsim = 1000)
      char_mat = apply(char_mat, 1, function(x) x -1) |> t()
      return(char_mat)
    }
    
    # simulate full morph character matrix
    char_mat_full = sim_bin_char(tree = tree_fp_complete,
                                 rate = rate_bin_evol, 
                                 nchar = n_char)
    ## subset character matrix
    # select specimens preserved from the adm
    char_mat_inc_A_rho0 = char_mat_full[fossils_inc_A$tip.label,] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_inc_A_rho0_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    char_mat_inc_A_rho1 = char_mat_full[c(fossils_inc_A$tip.label, recent_tips),] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_inc_A_rho1_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    char_mat_inc_B_rho0 = char_mat_full[fossils_inc_B$tip.label,] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_inc_B_rho0_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    char_mat_inc_B_rho1 = char_mat_full[c(fossils_inc_B$tip.label, recent_tips),] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_inc_B_rho1_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    # select equal no. of specimens, continuous sampling
    char_mat_cont_rho0 = char_mat_full[fossils_cont$tip.label,] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_cont_rho0_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    char_mat_cont_rho1 = char_mat_full[c(fossils_cont$tip.label, recent_tips),] |>
      ape::write.nexus.data(file = paste0(path,"char_mat_cont_rho1_nchar", n_char, "_",id,".nex"),
                            format = "standard")
    }
    
    #### simulate molecular data ####
    opts = paste0(mod, " -l", length_alignment)
    a = phyclust::seqgen(opts = opts, rooted.tree = tree_fp_complete) |>
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
skyline = get_skyline_breakpoints(adm_A_2km, hiat_min = hiat_min)

skyline_bp_age = t_max - skyline |> rev()



