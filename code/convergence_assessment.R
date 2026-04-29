ids = seq(1,40)
nchars = c(30, 100, 300, 1000)
analyses = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior")

check_file_existence = function(){

  cases = c("miller", "sinusoid")
  runs = c(1,2)
  
  names = c("id", "nchar", "analysis", "complete")
  df = data.frame(matrix(nrow = 0, ncol = length(names)))
  names(df) = names
  
  path = "data/fbd_strat"
  for (id in ids){
    for (nchar in nchars){
      for (case in cases){
        tree_sim = paste0(path, "/sim/tree_",id, "_nchar",nchar, "_", case, ".nex")
        tree_original = paste0(path, "/sim/original_tree_",id, "_nchar",nchar, "_", case, ".nex")
        original_fossils = paste0(path, "/sim/original_fossils_",id, "_nchar",nchar, "_", case, ".csv")
        mol_dat = paste0(path, "/sim/mol_dat_",id, "_nchar",nchar, "_", case, ".nex")
        fossils = paste0(path, "/sim/fossils_",id, "_nchar",nchar, "_", case, ".csv")
        char_mat = paste0(path, "/sim/char_mat_",id, "_nchar",nchar, "_", case, ".nex")
        log1 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_", case, "_run_", runs[1], ".log")
        tree1 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,  "_", case, "_run_", runs[2], ".log")
        log2 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_", case, "_run_", runs[1], ".log")
        tree2 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,  "_", case, "_run_", runs[2], ".log")
        
        df2 = data.frame(
          id = id,
          nchar = nchar,
          analysis = paste0("strat_",case),
          complete = all(file.exists(tree_sim, tree_original, original_fossils,
                                     mol_dat, fossils, char_mat,
                                     log1, log2,  tree1, tree2))
        )
        df = rbind(df, df2)
      }
    }
  }
  
  path = "data/fbd_base"
  for (id in ids){
    for (nchar in nchars){
      tree_sim = paste0(path, "/sim/tree_",id, "_nchar",nchar, ".nex")
      tree_original = paste0(path, "/sim/original_tree_",id, "_nchar",nchar, ".nex")
      original_fossils = paste0(path, "/sim/original_fossils_",id, "_nchar",nchar,  ".csv")
      mol_dat = paste0(path, "/sim/mol_dat_",id, "_nchar",nchar, ".nex")
      fossils = paste0(path, "/sim/fossils_",id, "_nchar",nchar,  ".csv")
      char_mat = paste0(path, "/sim/char_mat_",id, "_nchar",nchar,  ".nex")
      log1 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_run_", runs[1], ".log")
      tree1 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,   "_run_", runs[2], ".log")
      log2 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_run_", runs[1], ".log")
      tree2 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,   "_run_", runs[2], ".log")
      
      df2 = data.frame(
        id = id,
        nchar = nchar,
        analysis = paste0("base"),
        complete = all(file.exists(tree_sim, tree_original, original_fossils,
                                   mol_dat, fossils, char_mat,
                                   log1, log2,  tree1, tree2))
      )
      df = rbind(df, df2)
    }
  }
  
  
  path = "data/fbd_gap_est"
  for (id in ids){
    for (nchar in nchars){
      for (case in c("sinusoid")){
        tree_sim = paste0(path, "/sim/tree_",id, "_nchar",nchar, "_", case, ".nex")
        tree_original = paste0(path, "/sim/original_tree_",id, "_nchar",nchar, "_", case, ".nex")
        original_fossils = paste0(path, "/sim/original_fossils_",id, "_nchar",nchar, "_", case, ".csv")
        mol_dat = paste0(path, "/sim/mol_dat_",id, "_nchar",nchar, "_", case, ".nex")
        fossils = paste0(path, "/sim/fossils_",id, "_nchar",nchar, "_", case, ".csv")
        char_mat = paste0(path, "/sim/char_mat_",id, "_nchar",nchar, "_", case, ".nex")
        log1 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_", case, "_run_", runs[1], ".log")
        tree1 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,  "_", case, "_run_", runs[2], ".log")
        log2 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_", case, "_run_", runs[1], ".log")
        tree2 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,  "_", case, "_run_", runs[2], ".log")
        
        df2 = data.frame(
          id = id,
          nchar = nchar,
          analysis = paste0("gap_est"),
          complete = all(file.exists(tree_sim, tree_original, original_fossils,
                                     mol_dat, fossils, char_mat,
                                     log1, log2,  tree1, tree2))
        )
        df = rbind(df, df2)
      }
    }
  }
  
  
  path = "data/fbd_gap_prior"
  for (id in ids){
    for (nchar in nchars){
      for (case in c("sinusoid")){
        tree_sim = paste0(path, "/sim/tree_",id, "_nchar",nchar, "_", case, ".nex")
        tree_original = paste0(path, "/sim/original_tree_",id, "_nchar",nchar, "_", case, ".nex")
        original_fossils = paste0(path, "/sim/original_fossils_",id, "_nchar",nchar, "_", case, ".csv")
        mol_dat = paste0(path, "/sim/mol_dat_",id, "_nchar",nchar, "_", case, ".nex")
        fossils = paste0(path, "/sim/fossils_",id, "_nchar",nchar, "_", case, ".csv")
        char_mat = paste0(path, "/sim/char_mat_",id, "_nchar",nchar, "_", case, ".nex")
        log1 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_", case, "_run_", runs[1], ".log")
        tree1 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,  "_", case, "_run_", runs[2], ".log")
        log2 = paste0(path, "/rb_output/", "num_", id, "_nchar", nchar,  "_", case, "_run_", runs[1], ".log")
        tree2 = paste0(path, "/rb_output/", "tree_", id, "_nchar", nchar,  "_", case, "_run_", runs[2], ".log")
        
        df2 = data.frame(
          id = id,
          nchar = nchar,
          analysis = paste0("gap_prior"),
          complete = all(file.exists(tree_sim, tree_original, original_fossils,
                                     mol_dat, fossils, char_mat,
                                     log1, log2,  tree1, tree2))
        )
        df = rbind(df, df2)
      }
    }
  }
  return(df)  
}

df_existence = check_file_existence()

get_ess = function(id, nchar, analysis, burnin = 0){
  if (analysis == "base"){
    path = "data/fbd_base/"
    file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_run_1.log")
    file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_run_2.log")
  }
  if (analysis == "strat_miller"){
    path = "data/fbd_strat/"
    file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_miller_run_1.log")
    file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_miller_run_2.log")
  }
  if (analysis == "strat_sinusoid"){
    path = "data/fbd_strat/"
    file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_1.log")
    file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_2.log")
  }
  if (analysis == "gap_est"){
    path = "data/fbd_gap_est/"
    file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_1.log")
    file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_2.log")
  }
  if (analysis == "gap_prior"){
    path = "data/fbd_gap_prior/"
    file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_1.log")
    file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_2.log")
  }
  tr1 = RevGadgets::readTrace(path = file_run1, burnin = burnin)
  ess1 = as.mcmc(tr1[[1]]) |> effectiveSize() |> t() |> as.data.frame()
  tr2 = RevGadgets::readTrace(path = file_run2, burnin = burnin)
  ess2 = as.mcmc(tr2[[1]]) |> effectiveSize() |> t() |> as.data.frame()
  
  comb_trace = RevGadgets::readTrace(path = c(file_run1, file_run2), burnin = burnin) |> 
    RevGadgets::combineTraces()
  
  ess_full = as.mcmc(comb_trace$combined) |> effectiveSize() |> t() |> as.data.frame()
  ess_df = rbindlist(list(ess1, ess2, ess_full), fill = TRUE)
  ess_df$run = c("1", "2", "comb")
  ess_df$id = rep(id, 3)
  ess_df$nchar = rep(nchar, 3)
  ess_df$analysis = rep(analysis, 3)
  return(ess_df)
}

get_ess_of_existing_files = function(){
  df_ess = data.frame()
  for (id in ids){
    print(id)
    for (nchar in nchars){
      for (analysis in analyses){
        if (df_existence$complete[df_existence$id == id & df_existence$nchar == nchar & df_existence$analysis == analysis]){
          df_ess = rbindlist(list(df_ess, get_ess(id, nchar, analysis)), fill = TRUE)
        }
      }
    }
  }
  return(df_ess)
}

df_ess = get_ess_of_existing_files()

save(df_ess, file = "ess.RData")


ess_threshold = 200
df_converged = expand.grid(nchars = nchars, id = ids, analysis = analyses, run_successful = NA, converged = NA)
for (id1 in ids){
  for (nchar1 in nchars){
    for (analysis1 in analyses){
      run_success = df_existence$complete[df_existence$nchar == nchar1 & df_existence$id == id1 & df_existence$analysis == analysis1]
      df_converged$run_successful[df_converged$nchars == nchar1 & df_converged$id == id1 & df_converged$analysis == analysis1] =  run_success
      df_converged$converged[df_converged$nchars == nchar1 & df_converged$id == id1 & df_converged$analysis == analysis1] = run_success
      if (run_success){
        df = df_ess |>
          filter(analysis == analysis1 & nchar == nchar1 & id == id1) |> 
          select(where(~any(.x > 0.001 & !is.na(.x)))) |>
          select(-run, -id, -nchar, -analysis)
        converged = all(df> ess_threshold)
        df_converged$converged[df_converged$nchars == nchar1 & df_converged$id == id1 & df_converged$analysis == analysis1] = converged
      }
    }
  }
}

# proportion of converged
df_converged |>
  group_by(nchars, analysis) |>
  summarise(n_converged = mean(converged))

# 
# name = "data/fbd_base/rb_output/tree_29_nchar300_run_2.log"
# file.copy("data/fbd_base/rb_output/tree_29_nchar300_run_2.log", "data/fbd_base/rb_output/tree_29_nchar300_run_2.trees")
# a = rwty::load.trees("data/fbd_base/rb_output/tree_29_nchar300_run_2.trees", format = "revbayes")
# a
# 
# essSplitFreq(a$trees)
# plotEssSplits(a)
# 
# aa = checkConvergence(listname)
# 
# aa = checkConvergence(list_files = c("data/fbd_base/rb_output/tree_29_nchar300_run_2.trees"))
