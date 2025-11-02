library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)

if (!dir.exists("figs/")) {
  dir.create("figs/", recursive = TRUE)
}

path_res = "output/fbd_base/"

source("code/constants.R")

# constants

n_rep = 20
ids = as.character(seq_len(n_rep))

alpha = 0.25
kappa = 5

case = "base"
# create empty dataframe for storage
colnames = c("id", "nchar","ext", "spec", "or_time", "psi", "br_rate_mol", "br_rate_morph", "kappa", "alpha", "case")
df = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df) = colnames

f = function(x){as.numeric(x[1] != x[2])}
covers_true_value = function(x, true_val){
  quantile(x - true_val, probs = c(0.05, 0.95)) |> unname() |> sign() |> f()
}

# read in data
for (id in ids){
  for (nchar in n_chars){
    tr1 = paste0(path_res,"num_", id, "_nchar", nchar,"_run_1.log")
    tr2 = paste0(path_res, "num_", id, "_nchar", nchar,"_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat(paste0("skipping file ", tr1, "\n"))
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ext_est = comb_trace[[1]]$extinction_rate |> covers_true_value(true_val = mu)
    spec_est = comb_trace[[1]]$speciation_rate |> covers_true_value(true_val = lambda)
    or_time_est = comb_trace[[1]]$origin_time |> covers_true_value(true_val = t_max)
    psi_est = comb_trace[[1]]$psi |> mean()
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> covers_true_value(true_val = clock_rate_morph)
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> covers_true_value(true_val = clock_rate_mol)
    kappa_est = comb_trace[[1]]$kappa |> covers_true_value(true_val = kappa)
    alpha_est = comb_trace[[1]]$alpha_mol |> covers_true_value(true_val = alpha)
    
    
    map = ape::read.tree(paste0(path_res, "tree_", id, "_nchar", nchar, "_MAP.tre"))
    truth = ape::read.nexus(paste0("data/sim/fbd_base/tree_", id, "_nchar", nchar , ".nex"))
    
    rf_dist = phangorn::RF.dist(map, truth, rooted = TRUE, normalize = TRUE)
    
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ext = ext_est,
                     spec = spec_est,
                     or_time = or_time_est,
                     psi = psi_est,
                     br_rate = br_rate_est,
                     rf_distance = rf_dist,
                     br_rate_mol = br_rate_mol_est,
                     alpha = alpha_est,
                     kappa = kappa_est,
                     case = "base")
    df = rbind(df, df2)
  }
}

case = "miller"
path_res = "output/fbd_strat/"

for (id in ids){
  for (nchar in n_chars){
    tr1 = paste0(path_res,"num_", id, "_nchar", nchar, "_", case, "_run_1.log")
    tr2 = paste0(path_res, "num_", id, "_nchar", nchar,"_", case, "_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat(paste0("skipping file ", tr1, "\n"))
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ext_est = comb_trace[[1]]$extinction_rate |> covers_true_value(true_val = mu)
    spec_est = comb_trace[[1]]$speciation_rate |> covers_true_value(true_val = lambda)
    or_time_est = comb_trace[[1]]$origin_time |> covers_true_value(true_val = t_max)
    psi_est = comb_trace[[1]]$psi |> mean()
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> covers_true_value(true_val = clock_rate_morph)
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> covers_true_value(true_val = clock_rate_mol)
    kappa_est = comb_trace[[1]]$kappa |> covers_true_value(true_val = kappa)
    alpha_est = comb_trace[[1]]$alpha_mol |> covers_true_value(true_val = alpha)
    
    
    map = ape::read.tree(paste0(path_res, "tree_", id, "_nchar", nchar, "_", case,  "_MAP.tre"))
    truth = ape::read.nexus(paste0("data/sim/fbd_strat/tree_", id, "_nchar", nchar ,"_", case,  ".nex"))
    
    rf_dist = phangorn::RF.dist(map, truth, rooted = TRUE, normalize = TRUE)
    
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ext = ext_est,
                     spec = spec_est,
                     or_time = or_time_est,
                     psi = psi_est,
                     br_rate = br_rate_est,
                     rf_distance = rf_dist,
                     br_rate_mol = br_rate_mol_est,
                     alpha = alpha_est,
                     kappa = kappa_est,
                     case = "miller")
    df = rbind(df, df2)
  }
}


case = "sinusoid"


for (id in ids){
  for (nchar in n_chars){
    tr1 = paste0(path_res,"num_", id, "_nchar", nchar, "_", case, "_run_1.log")
    tr2 = paste0(path_res, "num_", id, "_nchar", nchar,"_", case, "_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat(paste0("skipping file ", tr1, "\n"))
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ext_est = comb_trace[[1]]$extinction_rate |> covers_true_value(true_val = mu)
    spec_est = comb_trace[[1]]$speciation_rate |> covers_true_value(true_val = lambda)
    or_time_est = comb_trace[[1]]$origin_time |> covers_true_value(true_val = t_max)
    psi_est = comb_trace[[1]]$psi |> mean()
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> covers_true_value(true_val = clock_rate_morph)
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> covers_true_value(true_val = clock_rate_mol)
    kappa_est = comb_trace[[1]]$kappa |> covers_true_value(true_val = kappa)
    alpha_est = comb_trace[[1]]$alpha_mol |> covers_true_value(true_val = alpha)
    
    
    map = ape::read.tree(paste0(path_res, "tree_", id, "_nchar", nchar, "_", case,  "_MAP.tre"))
    truth = ape::read.nexus(paste0("data/sim/fbd_strat/tree_", id, "_nchar", nchar ,"_", case,  ".nex"))
    
    rf_dist = phangorn::RF.dist(map, truth, rooted = TRUE, normalize = TRUE)
    
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ext = ext_est,
                     spec = spec_est,
                     or_time = or_time_est,
                     psi = psi_est,
                     br_rate = br_rate_est,
                     rf_distance = rf_dist,
                     br_rate_mol = br_rate_mol_est,
                     alpha = alpha_est,
                     kappa = kappa_est,
                     case = "sinusoid")
    df = rbind(df, df2)
  }
}


## ext rate ##

ext_mat = matrix(data = NA, nrow = 5, ncol = 3, dimnames = list("nchar" = c(30, 100, 300, 1000, "all"), "case" = c("base", "miller", "sinusoid")))
for (case in c("base", "miller", "sinusoid")){
  df1 = df[df$case == case, ]
  ext_mat["all", case] = sum(df1[["ext"]])/ length(df1[["ext"]])
  for (nchar in c(30, 100, 300, 1000)){
    df2 = df1[df1$nchar == nchar,]
    ext_mat[as.character(nchar), case] = sum(df2[["ext"]])/ length(df2[["ext"]])
  }
}

get_coverage_ratio = function(name){
  mat = matrix(data = NA, nrow = 5, ncol = 3, dimnames = list("nchar" = c(30, 100, 300, 1000, "all"), "case" = c("base", "miller", "sinusoid")))
  for (case in c("base", "miller", "sinusoid")){
    df1 = df[df$case == case, ]
    mat["all", case] = sum(df1[[name]])/ length(df1[[name]])
    for (nchar in c(30, 100, 300, 1000)){
      df2 = df1[df1$nchar == nchar,]
      mat[as.character(nchar), case] = sum(df2[[name]])/ length(df2[[name]])
    }
  } 
  return(mat)
}

poss_names = names(df)[c(3,4,7,9,10, 11)]

ext_cov_ratio  = get_coverage_ratio(poss_names[1])
ext_cov_ratio

or_cov_ratio  = get_coverage_ratio(poss_names[2])
or_cov_ratio

morph_cov_ratio  = get_coverage_ratio(poss_names[3])
morph_cov_ratio

mol_cov_ratio  = get_coverage_ratio(poss_names[4])
mol_cov_ratio

alpha_cov_ratio  = get_coverage_ratio(poss_names[5])
alpha_cov_ratio

kappa_cov_ratio  = get_coverage_ratio(poss_names[6])
kappa_cov_ratio

# 
# # convert columns to appropriate types
# df$id = as.numeric(df$id)
# df$nchar = factor(df$nchar, levels = n_chars)
# df$ext_rel = (mu - df$ext)/mu
# df$spec_rel = (lambda - df$spec)/lambda
# df$or_time_rel = (t_max- df$or_time ) / t_max
# df$br_rate_rel = (clock_rate_morph - df$br_rate )/ clock_rate_morph
# df$br_rate_mol_rel = (clock_rate_mol - df$br_rate_mol)/clock_rate_mol
# df$alpha_rel = (alpha -df$alpha)/alpha
# df$kappa_rel = (kappa - df$kappa)/ kappa
# df$case = factor(df$case, levels = c("base", "miller", "sinusoid"))
# 
# p1 = df |> ggplot(aes(x = nchar, y = ext_rel, fill = case)) +
#   geom_boxplot() +
#   labs(title = "Extinction Rate", x = "# Characters", y = "Relative error")
# p1
# p2 = df |> ggplot(aes(x = nchar, y = spec_rel, fill = case)) +
#   geom_boxplot() +
#   labs(title = "Speciation Rate", x = "# Characters", y = "Relative error")
# p2
# p3 = df |> ggplot(aes(x = nchar, y = or_time_rel, fill = case)) +
#   geom_boxplot() +
#   labs(title = "Origin Time", x = "# Characters", y = "Relative error") +
#   ylim(c(0, -0.2))
# p3
# 
# p4 = df |> ggplot(aes(x = nchar, y = br_rate_rel, fill = case)) +
#   geom_boxplot() +
#   labs(title = "Branch Rate Morph", x = "# Characters", y = "Relative error") + 
#   ylim(c(-1, 1))
# p4
# 
# p5 = df |> ggplot(aes(x = nchar, y = rf_distance, fill = case)) +
#   geom_boxplot() + 
#   labs(title = "RF distance", x = "# Characters", y = "RF distance")
# p5
# 
# p6 = df |> ggplot(aes(x = nchar, y = br_rate_mol_rel, fill = case))+
#   geom_boxplot() + 
#   labs(title = "Mol clock rate", x = "# Characters", y = "rel error Mol clock rate")
# p6
# 
# p7 = df |> ggplot(aes(x = nchar, y = kappa_rel, fill = case))+
#   geom_boxplot() +
#   labs(title = "kappa", x = "# characters", y = "relative error")
# p7
# 
# p8 = df |> ggplot(aes(x = nchar, y = alpha_rel, fill = case))+
#   geom_boxplot() +
#   labs(title = "alpha", x = "# characters", y = "relative error")
# p8
# 
# fig = ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 3, nrow = 3, common.legend = TRUE)
# fig
# ggsave("figs/fbd_comp.png", plot = fig)
# 
# # df |> ggplot(aes(x = nchar, y = br_rate_mol_rel)) +
# #   geom_boxplot() +
# #   labs(title = "Branch Rate Mol", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = alpha_rel)) +
# #   geom_boxplot() +
# #   labs(title = "Alpha", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = kappa_rel)) +
# #   geom_boxplot() +
# #   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = psi1)) +
# #   geom_boxplot() +
# #   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = psi2)) +
# #   geom_boxplot() +
# #   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = psi3)) +
# #   geom_boxplot() +
# #   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = psi4)) +
# #   geom_boxplot() +
# #   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# # 
# # df |> ggplot(aes(x = nchar, y = psi5)) +
# #   geom_boxplot() +
# #   labs(title = "Kappa", x = "# Characters", y = "Relative error")