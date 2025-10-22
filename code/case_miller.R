# simulation for the miller case
case = "miller"

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
    
    # save original sim tree without modifications bc I don't trust myself
    tree_original_name = paste0(path,"original_tree_", id , "_nchar", n_char, "_", case,  ".nex")
    ape::write.nexus(tree_complete,
                     file= tree_original_name)
    
    if(case == "sinusoid"){adm = adm_sinusoid}
    if(case == "miller"){adm = adm_miller}
    pres = identity
    ctc = get_collection_prob(adm)
    
    # full fossil record
    fossils_full_base = FossilSim::sim.fossils.poisson(rate = sampling_rate,
                                                       tree = tree_complete) |>
      admtools::rev_dir(ref = t_max) |> 
      StratPal::apply_taphonomy(pres_potential = pres, ctc =  ctc) |>
      admtools::rev_dir(ref = t_max) |>
      subsample_fossils(n = n_fossils)
    
    fossils_orig_filename = paste0(path,"original_fossils_", id , "_nchar", n_char, "_", case,  ".csv")
    write.csv(fossils_full_base,
              file = fossils_orig_filename)
    
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
    fossil_name = paste0(path,"fossils_", id , "_nchar", n_char, "_", case,  ".csv")
    fossils_full |>
      get_fossil_ages(extant_taxa = extant_tips)|>
      write.table(file = fossil_name,
                  quote = FALSE,
                  row.names = FALSE)
    
    tree_name = paste0(path,"tree_", id , "_nchar", n_char, "_", case,  ".nex")
    ape::write.nexus(tree_fp_rho1,
                     file= tree_name)
    
    tree_w_rate_morph = tree_fp_rho1
    tree_w_rate_morph$edge.length = tree_w_rate_morph$edge.length * clock_rate_morph
    char_mat = sim_bin_char(tree = tree_w_rate_morph,
                            nchar = n_char)
    char_mat_name = paste0(path,"char_mat_", id , "_nchar", n_char, "_", case,  ".nex")
    ape::write.nexus.data(char_mat,
                          file = char_mat_name,
                          format = "standard")
    #plot(tree_fp_rho1)
    
    ## molecular alignment
    opts = paste0(mod, " -l", length_alignment)
    
    tree_w_rate_mol = tree_fp_rho1
    tree_w_rate_mol$edge.length = tree_w_rate_mol$edge.length * clock_rate_mol
    mol = sim_mol_char(tree = tree_w_rate_mol, extant_tips = extant_tips, opts = opts)
    ape::write.nexus.data(mol,
                          file = paste0(path, "mol_dat_", id, "_nchar", n_char, "_", case,  ".nex"))
  }
}


cat(paste0(n_succ, " successful trees\n"))
cat(paste0(n_fail, " failed trees\n"))

cat("Done \n")