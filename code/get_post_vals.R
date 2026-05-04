load("data/convergence_assessment.RData")
source("code/constants.R")

library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)
library(paleotree)
library(dispRity)
library(treeio)
library(dplyr)

ids = c(1:40)
nchars= c(30, 100, 300, 1000)
analyses = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior")
burnin = 0

get_divergence_times = function(map, ref){
  tol = 10^-5
  ages = suppressMessages(dateNodes(ref))
  extant_tips = ref$tip.label[ages < tol]
  df = data.frame(matrix(nrow = 0, ncol = 5))
  names(df) = c("tip1", "tip2", "true_age", "min_age", "max_age")
  if(length(extant_tips)< 2){
    warning("need at least 2 extant tips")
    return(df)
  }
  mrca_index = c()
  mrca_list = list()
  #print(extant_tips)
  k = 1
  for (i in 1:(length(extant_tips)-1)){
    for (j in (i+1):length(extant_tips)){
      tips = c(extant_tips[i], extant_tips[j])
      #print(tips)
      mrca = getMRCA(ref, tips)
      if (! mrca %in% mrca_index){
        mrca_index = c(mrca_index, mrca)
        mrca_list[[k]] = tips
        k = k+1
      }
    }
  }
  for (i in seq_along(mrca_index)){
    tips = mrca_list[[i]]
    true_age = ages[mrca_index[i]] |> unname()
    t_hdp = map@data$age_0.95_HPD[map@data$node == getMRCA(map@phylo, tips)][[1]]
    df2 = data.frame(tip1 = tips[1],
                     tip2 = tips[2],
                     true_age = true_age,
                     min_age = min(t_hdp),
                     max_age = max(t_hdp))
    df = rbind(df, df2)
  }
  df$covered = df$true_age >= df$min_age & df$true_age <= df$max_age
  return(df)
}


get_div_statistics = function(df){
  mean_prec = c((df$max_age - df$min_age)/df$true_age) |> mean()
  cov_freq = mean(df$covered)
  return(c("mean_prec" = mean_prec,
           "cov_freq" = cov_freq))
}

get_sa_prob = function(map, ref){
  tol = 10^-5
  extant_tips = ref$tip.label[suppressMessages(dateNodes(ref)) < tol]
  not_extant = ref$tip.label[! ref$tip.label %in% extant_tips ]
  not_extant_node = which(map@phylo$tip.label %in% not_extant)
  sas = map@data$sampled_ancestor[map@data$node %in% not_extant_node]
  sas[is.na(sas)] = 0
  names(sas) = map@phylo$tip.label[not_extant_node]
  
  sa_true_names = ref$tip.label[ref$edge[ref$edge.length == 0, 2]]
  sa_true = rep(0, length(not_extant))
  names(sa_true) = not_extant
  sa_true[sa_true_names] = 1
  
  return(list(sas, sa_true))
}

get_extant_tip_labels = function(tree, tol = 10^-5){
  return(tree$tip.label[dateNodes(tree) < tol])
}


mean_sa_cov_freq = function(l) {sapply(names(l[[1]]), function(x) abs(l[[2]][x] - l[[1]][x]) < 0.05) |> mean()}


tree_params = function(map, truth){
  df = get_divergence_times(map = map,
                            ref = truth)
  n = 3
  mean_counted = NA
  cov_freq_counted = NA
  if (length(df$true_age) >= 3){
    mean_counted = c((df$max_age - df$min_age)/df$true_age) |> mean()
    cov_freq_counted = mean(df$covered)
  }
  
  x = c("mean_prec_div_times" = c((df$max_age - df$min_age)/df$true_age) |> mean(),
        "cov_freq_div_times" = mean(df$covered),
        "mean_sa_cov_freq" = get_sa_prob(map = map,
                                         ref = truth) |>
          mean_sa_cov_freq(),
        "rf_dist" = phangorn::RF.dist(map@phylo, truth, rooted = TRUE, normalize = TRUE),
        "mean_counted" = mean_counted,
        "cov_freq_counted" = cov_freq_counted)
  return(x)
}



get_posterior_vals = function(){
  df_median = data.frame()
  
  for (id in ids){
    print(id)
    for (nchar in nchars){
      for (analysis in analyses){
        if (df_converged$converged[df_converged$nchars == nchar & df_converged$id == id & df_converged$analysis == analysis]){
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
          comb_trace = RevGadgets::readTrace(path = c(file_run1, file_run2), burnin = burnin) |> 
            RevGadgets::combineTraces()
          df_temp = as.data.frame(lapply(comb_trace$combined, function(x) quantile(x, p = c(0.05, 0.5, 0.95))))
          df_te = c(t(df_temp))
          names(df_te) = paste(rep(colnames(df_temp), times = nrow(df_temp)),rep(rownames(df_temp), each = ncol(df_temp)), sep = "_")
          df_te = df_te[order(names(df_te))]
          df_te["id"] = id
          df_te["nchar"] = nchar
          df_te["analysis"] = analysis
          df_median = data.table::rbindlist(list(df_median, as.data.frame(as.list(df_te))), fill = TRUE)
        }
      }
    }
  }
  return(df_median)
}

get_tree_statistics = function(){
  df_stat = data.frame()
  
  for (id in ids){
    print(id)
    for (nchar in nchars){
      for (analysis in analyses){
        if (df_converged$converged[df_converged$nchars == nchar & df_converged$id == id & df_converged$analysis == analysis]){
          if (analysis == "base"){
            path = "data/fbd_base/"
            map_path = paste0(path, "rb_output/",  "tree_", id, "_nchar", nchar, "_MAP.tre")
            truth_path = paste0(path, "sim/tree_", id, "_nchar", nchar , ".nex")
          }
          if (analysis == "strat_miller"){
            path = "data/fbd_strat/"
            map_path = paste0(path, "rb_output/",  "tree_", id, "_nchar", nchar, "_miller_MAP.tre")
            truth_path = paste0(path, "sim/tree_", id, "_nchar", nchar , "_miller.nex")
          }
          if (analysis == "strat_sinusoid"){
            path = "data/fbd_strat/"
            map_path = paste0(path, "rb_output/",  "tree_", id, "_nchar", nchar, "_sinusoid_MAP.tre")
            truth_path = paste0(path, "sim/tree_", id, "_nchar", nchar , "_sinusoid.nex")
          }
          if (analysis == "gap_est"){
            path = "data/fbd_gap_est/"
            map_path = paste0(path, "rb_output/",  "tree_", id, "_nchar", nchar, "_sinusoid_MAP.tre")
            truth_path = paste0(path, "sim/tree_", id, "_nchar", nchar , "_sinusoid.nex")
          }
          if (analysis == "gap_prior"){
            path = "data/fbd_gap_prior/"
            map_path = paste0(path, "rb_output/",  "tree_", id, "_nchar", nchar, "_sinusoid_MAP.tre")
            truth_path = paste0(path, "sim/tree_", id, "_nchar", nchar , "_sinusoid.nex")
          }
          if (!all(file.exists(map_path, truth_path))){stop()}
          
          truth = ape::read.nexus(truth_path)
          map = suppressWarnings(treeio::read.beast.newick(map_path)) 
          
          df_te = tree_params(map, truth) |> as.list() |> as.data.frame()
          df_te$id = id
          df_te$nchar = nchar
          df_te$analysis = analysis
          
          df_stat = data.table::rbindlist(list(df_stat, df_te))
        }
      }
    }
  }
  return(df_stat)
}

df_stat = get_tree_statistics()

df_median = get_posterior_vals()


df_median$ext_covered = as.numeric(df_median$extinction_rate_5.) < mu & as.numeric(df_median$extinction_rate_95.) > mu
df_median$spec_covered = as.numeric(df_median$speciation_rate_5.) < lambda & as.numeric(df_median$speciation_rate_95.) > lambda
df_median$orig_covered = as.numeric(df_median$origin_time_5.) < t_max & as.numeric(df_median$origin_time_95.) > t_max
df_median$branch_rates_mol_covered = as.numeric(df_median$branch_rates_mol_5.) < clock_rate_mol & as.numeric(df_median$branch_rates_mol_95.) > clock_rate_mol
df_median$branch_rates_morpho_covered = as.numeric(df_median$branch_rates_morpho_5.) < clock_rate_morph & as.numeric(df_median$branch_rates_morpho_95.) > clock_rate_morph

df_median$spec_rel_error = (as.numeric(df_median$speciation_rate_50.) - lambda)/lambda
df_median$ext_rel_error = (as.numeric(df_median$extinction_rate_50.) - mu)/ mu
df_median$origin_rel_error = (as.numeric(df_median$origin_time_50.) - t_max)/t_max
df_median$branch_rates_mol_rel_error = (as.numeric(df_median$branch_rates_mol_50.)- clock_rate_mol)/clock_rate_mol
df_median$branch_rates_morpho_rel_error = (as.numeric(df_median$branch_rates_morpho_50.)- clock_rate_morph)/clock_rate_morph

save(df_median, df_stat, file = "data/post_values.RData")

