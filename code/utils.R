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


get_fossil_ages = function(tree, t_max){
  #' @title get ages of samples from sampled tree
  #' 
  #' @param tree sampled tree
  #' @param t_max length of simulation
  #' @param signif significant digits of ages
  #' 
  x = pmax(t_max -  ape::node.depth.edgelength(tree)[seq_len(length(tree$tip.label))], 0)
  df = data.frame(taxon = tree$tip.label,
                  min_age = x ,
                  max_age = x)
  return(df)
}

sim_bin_char = function(tree, rate, nchar){
  #' @title simulate traits along tree
  #' 
  #' @param tree sampled tree
  #' @param rate rate of evolution
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

subsample_fossils = function(fossils, n){
  #' randomly select a number of fossils
  ind = sample(seq_len(length(fossils$sp)), size = n)
  return(fossils[ind,])
}