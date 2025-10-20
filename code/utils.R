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