load("cont_params_values.RData")

library(dplyr)
library(ggplot2)

df = df_median 
df$analysis = factor(df$analysis)
df$nchar = factor(df$nchar, levels = n_chars, ordered = TRUE)



df |> ggplot(aes(x = nchar, y = spec_rel_error, fill = analysis)) +
  geom_boxplot() +
  geom_hline(yintercept = 0)

df |> ggplot(aes(x = nchar, y = ext_rel_error, fill = analysis)) +
  geom_boxplot() +
  geom_hline(yintercept = 0)

plot_base_vs_strat = function(){
  p1 = df |> 
    filter(analysis %in% c("base", "strat_miller", "strat_sinusoid")) |>
    ggplot(aes(x = nchar, y = ext_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error extinction")
  p2 = df |> 
    filter(analysis %in% c("base", "strat_miller", "strat_sinusoid")) |>
    ggplot(aes(x = nchar, y = spec_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error speciation")
  p3 = df |> 
    filter(analysis %in% c("base", "strat_miller", "strat_sinusoid")) |>
    ggplot(aes(x = nchar, y = origin_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error origin")
  
  p4 = df |> 
    filter(analysis %in% c("base", "strat_miller", "strat_sinusoid")) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error mol branch rates")
  
  p5 = df |> 
    filter(analysis %in% c("base", "strat_miller", "strat_sinusoid")) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error morpho branch rates")
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, common.legend = TRUE)
  return(p)
}

p = plot_base_vs_strat()
p
plot_sampling_strategies_comp = function(){
  selection = c("strat_miller", "gap_est", "gap_prior")
  p1 = df |> 
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = ext_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error extinction")
  p2 = df |> 
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = spec_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error speciation")
  p2
  p3 = df |> 
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = origin_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error origin")
  
  p4 = df |> 
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error mol branch rates")
  
  p5 = df |> 
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error morpho branch rates")
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, common.legend = TRUE)
}
p = plot_sampling_strategies_comp()
p

plot_coverage_freq_sampling = function(){
  selection = c("base", "strat_miller", "strat_sinusoid")
  p1 = df |> select(spec_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(spec_coverage_freq = mean(spec_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = spec_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p2 = df |> select(ext_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(ext_coverage_freq = mean(ext_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = ext_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p3 = df |> select(branch_rates_mol_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_mol_coverage_freq = mean(branch_rates_mol_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p4 = df |> select(branch_rates_morpho_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_morpho_coverage_freq = mean(branch_rates_morpho_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, common.legend = TRUE)
  return(p)
  
}

p = plot_coverage_freq_sampling()
p
ggsave(filename = "figs/coverage_freq_model_violation.png",
       plot = plot_coverage_freq_sampling())

plot_coverage_freq_sampling_strategies = function(){
  selection = c("strat_sinusoid", "gap_est", "gap_prior")
  p1 = df |> select(spec_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(spec_coverage_freq = mean(spec_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = spec_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p2 = df |> select(ext_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(ext_coverage_freq = mean(ext_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = ext_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p3 = df |> select(branch_rates_mol_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_mol_coverage_freq = mean(branch_rates_mol_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p4 = df |> select(branch_rates_morpho_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_morpho_coverage_freq = mean(branch_rates_morpho_covered), .groups = "drop") |>
    filter(analysis %in% selection) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, common.legend = TRUE)
  return(p)
  
}

p = plot_coverage_freq_sampling_strategies()
p
ggsave(filename = "figs/coverage_freq_sampling_info.png",
       plot = plot_coverage_freq_sampling_strategies())
