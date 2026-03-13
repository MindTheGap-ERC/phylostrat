library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)
library(ggplot2)

source("code/constants.R")
alpha = 0.25
kappa = 5


runs = 1:2 |> as.character()
reps = 1:4 |> as.character()
nchars = c(30, 100, 300, 1000) |> as.character()
scenarios = c("continuous", "big_gap", "short_gap")

missing_files = c()

colnames = c("id", "nchar","ext", "spec", "or_time", "psi", "br_rate_mol", "br_rate_morph", "kappa", "alpha", "case")
df = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df) = colnames

for (scenario in scenarios){
  for (nchar in nchars){
    for (rep in reps){
        if (scenario == "continuous"){
          path = "data/fbd_base/rb_output/"
          scen_add = ""
        }
        else{
          path = "data/fbd_strat/rb_output/"
          if (scenario == "big_gap"){
            scen_add = "sinusoid_"
          }
          else{
            scen_add = "miller_"
          }
        }
        file_name_1 = paste0(path, "num_", rep, "_nchar", nchar, "_",scen_add, "run_1", ".log")
        file_name_2 = paste0(path, "num_", rep, "_nchar", nchar, "_",scen_add, "run_2", ".log")
        if (!file.exists(file_name_1) | !file.exists(file_name_2)){
          missing_files = c(missing_files, file_name_1, file_name_2)
          break
        }
        comb_trace = RevGadgets::readTrace(path = c(file_name_1, file_name_2), burnin = 0) |> 
          RevGadgets::combineTraces()
        ext_est = comb_trace[[1]]$extinction_rate |> mean() |> unname()
        spec_est = comb_trace[[1]]$speciation_rate |> mean()
        or_time_est = comb_trace[[1]]$origin_time |> mean()
        br_rate_est = comb_trace[[1]]$branch_rates_morpho |> mean()
        br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> mean()
        kappa_est = comb_trace[[1]]$kappa |> mean()
        alpha_est = comb_trace[[1]]$alpha_mol |> mean()
        
        df2 = data.frame(id = rep,
                         nchar = nchar,
                         ext = ext_est,
                         spec = spec_est,
                         or_time = or_time_est,
                         br_rate = br_rate_est,
                         br_rate_mol = br_rate_mol_est,
                         alpha = alpha_est,
                         kappa = kappa_est,
                         case = scenario)
        df = rbind(df, df2)
    }
  }
}

df$id = as.numeric(df$id)
df$nchar = factor(df$nchar, levels = n_chars)
df$ext_rel = (mu - df$ext)/mu
df$spec_rel = (lambda - df$spec)/lambda
df$or_time_rel = (t_max- df$or_time ) / t_max
df$br_rate_rel = (clock_rate_morph - df$br_rate )/ clock_rate_morph
df$br_rate_mol_rel = (clock_rate_mol - df$br_rate_mol)/clock_rate_mol
df$alpha_rel = (alpha -df$alpha)/alpha
df$kappa_rel = (kappa - df$kappa)/ kappa
df$sc

p1 = df |> ggplot(aes(x = nchar, y = ext_rel, fill = case)) +
  geom_boxplot() +
  labs(title = "Extinction Rate", x = "# Characters", y = "Relative error")
p1
p2 = df |> ggplot(aes(x = nchar, y = spec_rel, fill = case)) +
  geom_boxplot() +
  labs(title = "Speciation Rate", x = "# Characters", y = "Relative error")
p2
p3 = df |> ggplot(aes(x = nchar, y = or_time_rel, fill = case)) +
  geom_boxplot() +
  labs(title = "Origin Time", x = "# Characters", y = "Relative error")
p3

p4 = df |> ggplot(aes(x = nchar, y = br_rate_rel, fill = case)) +
  geom_boxplot() +
  labs(title = "Branch Rate Morph", x = "# Characters", y = "Relative error")
p4
p5 = df |> ggplot(aes(x = nchar, y = rf_distance, fill = case)) +
  geom_boxplot() + 
  labs(title = "RF distance", x = "# Characters", y = "RF distance")
p5

p6 = df |> ggplot(aes(x = nchar, y = br_rate_mol_rel))+
  geom_boxplot() + 
  labs(title = "Mol clock rate", x = "# Characters", y = "rel error Mol clock rate")
p6

p7 = df |> ggplot(aes(x = nchar, y = kappa_rel))+
  geom_boxplot() +
  labs(title = "kappa", x = "# characters", y = "relative error")
p7

p8 = df |> ggplot(aes(x = nchar, y = alpha_rel))+
  geom_boxplot() +
  labs(title = "alpha", x = "# characters", y = "relative error")
p8

fig = ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 3, nrow = 3)
fig
ggsave(paste0("figs/fbd_strat/fbd_strat_summary_", case, ".png"), plot = fig)

