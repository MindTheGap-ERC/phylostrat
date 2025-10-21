library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)

if (!dir.exists("figs/fbd_base")) {
  dir.create("figs/fbd_base", recursive = TRUE)
}

path_res = "output/fbd_base/"

source("code/constants.R")

# constants

n_rep = 10
ids = as.character(seq_len(n_rep))

alpha = 0.25
kappa = 5


# create empty dataframe for storage
colnames = c("id", "nchar","ext", "spec", "or_time", "psi", "br_rate_mol", "br_rate_morph", "kappa", "alpha")
df = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df) = colnames

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
    
    ext_est = comb_trace[[1]]$extinction_rate |> mean() |> unname()
    spec_est = comb_trace[[1]]$speciation_rate |> mean()
    or_time_est = comb_trace[[1]]$origin_time |> mean()
    psi_est = comb_trace[[1]]$psi |> mean()
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> mean()
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> mean()
    kappa_est = comb_trace[[1]]$kappa |> mean()
    alpha_est = comb_trace[[1]]$alpha_mol |> mean()
    
    
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
                     kappa = kappa_est)
    df = rbind(df, df2)
  }
}

# convert columns to appropriate types
df$id = as.numeric(df$id)
df$nchar = factor(df$nchar, levels = n_chars)
df$ext_rel = (mu - df$ext)/mu
df$spec_rel = (lambda - df$spec)/lambda
df$or_time_rel = (t_max- df$or_time ) / t_max
df$br_rate_rel = (clock_rate_morph - df$br_rate )/ clock_rate_morph
df$br_rate_mol_rel = (clock_rate_mol - df$br_rate_mol)/clock_rate_mol
df$alpha_rel = (alpha -df$alpha)/alpha
df$kappa_rel = (kappa - df$kappa)/ kappa

p1 = df |> ggplot(aes(x = nchar, y = ext_rel)) +
  geom_boxplot() +
  labs(title = "Extinction Rate", x = "# Characters", y = "Relative error")
p1
p2 = df |> ggplot(aes(x = nchar, y = spec_rel)) +
  geom_boxplot() +
  labs(title = "Speciation Rate", x = "# Characters", y = "Relative error")
p2
p3 = df |> ggplot(aes(x = nchar, y = or_time_rel)) +
  geom_boxplot() +
  labs(title = "Origin Time", x = "# Characters", y = "Relative error")
p3

p5 = df |> ggplot(aes(x = nchar, y = br_rate_rel)) +
  geom_boxplot() +
  labs(title = "Branch Rate Morph", x = "# Characters", y = "Relative error")
p5
p6 = df |> ggplot(aes(x = nchar, y = rf_distance)) +
  geom_boxplot() + 
  labs(title = "RF distance", x = "# Characters", y = "RF distance")
p6

p7 = df |> ggplot(aes(x = nchar, y = br_rate_mol_rel))+
  geom_boxplot() + 
  labs(title = "Mol clock rate", x = "# Characters", y = "rel error Mol clock rate")
p7

p8 = df |> ggplot(aes(x = nchar, y = kappa_rel))+
  geom_boxplot() +
  labs(title = "kappa", x = "# characters", y = "relative error")
p8

p9 = df |> ggplot(aes(x = nchar, y = alpha_rel))+
  geom_boxplot() +
  labs(title = "alpha", x = "# characters", y = "relative error")
p9

fig = ggarrange(p1, p2, p3, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)
fig
ggsave("figs/fbd_base/fbd_base_summary.png", plot = fig)

# df |> ggplot(aes(x = nchar, y = br_rate_mol_rel)) +
#   geom_boxplot() +
#   labs(title = "Branch Rate Mol", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = alpha_rel)) +
#   geom_boxplot() +
#   labs(title = "Alpha", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = kappa_rel)) +
#   geom_boxplot() +
#   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = psi1)) +
#   geom_boxplot() +
#   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = psi2)) +
#   geom_boxplot() +
#   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = psi3)) +
#   geom_boxplot() +
#   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = psi4)) +
#   geom_boxplot() +
#   labs(title = "Kappa", x = "# Characters", y = "Relative error")
# 
# df |> ggplot(aes(x = nchar, y = psi5)) +
#   geom_boxplot() +
#   labs(title = "Kappa", x = "# Characters", y = "Relative error")