library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)
library(paleotree)
library(dispRity)
library(treeio)
library(dplyr)

source("code/constants.R")


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
  
  x = c("mean_prec_div_times" = c((df$max_age - df$min_age)/df$true_age) |> mean(),
        "cov_freq_div_times" = mean(df$covered),
        "mean_sa_cov_freq" = get_sa_prob(map = map,
                                         ref = truth) |>
          mean_sa_cov_freq(),
        "rf_dist" = phangorn::RF.dist(map@phylo, truth, rooted = TRUE, normalize = TRUE))
  return(x)
}

id = 12
nchar = 1000
map_path = paste0("output/fbd_base/", "tree_", id, "_nchar", nchar, "_MAP.tre")
truth_path = paste0("data/sim/fbd_base/tree_", id, "_nchar", nchar , ".nex")

truth = ape::read.nexus(truth_path)
map = treeio::read.beast.newick(map_path) 

plot(map@phylo)
plot(truth)

tree_params(map, truth)

load("converged_runs.RData")


ids = c(1:40)
nchars= c(30, 100, 300, 1000)
analyses = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior")

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


scen_model_violation = c("base", "strat_miller", "strat_sinusoid")
scen_sampling_strategy = c("strat_sinusoid", "gap_est", "gap_prior")

plot_tree_stats_model_violation = function(){
  p1 = df_stat |>
    filter(!is.na(rf_dist)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = rf_dist, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    labs(title = "RF distance MAP - truth")
  
  p2 = df_stat |>
    filter(!is.na(mean_prec_div_times)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = mean_prec_div_times, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    ylim(0,2) +
    labs(title = "Mean precision divergence time")
  
  p3 = df_stat |>
    filter(!is.na(cov_freq_div_times)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = cov_freq_div_times, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    labs(title = "coverage frequency divergence times")
  
  p4 = df_stat |>
    filter(!is.na(mean_sa_cov_freq)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = mean_sa_cov_freq, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    ylim(0,1) + 
    labs(title = "mean SA coverage frequency")
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE)
}

p = plot_tree_stats_model_violation()
p
ggsave(filename = "figs/tree_statistics_model_violation.png",
       plot = plot_tree_stats_model_violation())

plot_tree_stats_sampling_strategy = function(){
  p1 = df_stat |>
    filter(!is.na(rf_dist)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = rf_dist, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    labs(title = "RF distance MAP - truth")
  
  p2 = df_stat |>
    filter(!is.na(mean_prec_div_times)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = mean_prec_div_times, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    ylim(0,2) +
    labs(title = "Mean precision divergence time")
  
  p3 = df_stat |>
    filter(!is.na(cov_freq_div_times)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = cov_freq_div_times, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    labs(title = "coverage frequency divergence times")
  
  p4 = df_stat |>
    filter(!is.na(mean_sa_cov_freq)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = mean_sa_cov_freq, color = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    ylim(0,1) + 
    labs(title = "mean SA coverage frequency")
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE)
}

p = plot_tree_stats_sampling_strategy()
p

ggsave(filename = "figs/tree_statistics_sampling_strategy.png",
       plot = plot_tree_stats_sampling_strategy())
