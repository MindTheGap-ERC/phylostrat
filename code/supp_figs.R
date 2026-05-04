source("code/constants.R")
load("data/post_values.RData")
load("data/convergence_assessment.RData")

library(dplyr)
library(scales)

scen_axis_labels = c("base" = "cFBD +\nCFS",
                     "strat_miller" = "cFBD +\nMSG",
                     "strat_sinusoid" = "cFBD +\nFLG",
                     "gap_est" = "sFBDw +\nFLG",
                     "gap_prior" = "sFBDs +\nFLG")
nchar_axis_label = "# Char"

plot_n_converged_runs = function(){
  p = df_converged |>
    group_by(analysis, nchars) |>
    summarise(n_converged = sum(converged), .groups = "drop") |>
    mutate(nchars = as.factor(nchars)) |>
    ggplot(aes(x = analysis, y = n_converged, color = nchars)) +
    geom_point(position = position_jitter(height = 0, width = 0.15)) +
    scale_y_continuous(
      breaks = breaks_width(1),
      labels = as.integer,
      limits = c(30, 40)
    ) +
    labs(title = "# Converged runs",
         x = "Analysis",
         y = "n",
         color = nchar_axis_label) +
    scale_x_discrete(labels = scen_axis_labels)
  return(p)
}
ggsave(filename = "figs/post_analysis/converged_runs.png",
       plot = plot_n_converged_runs())


#### Tree figures ####

scen_model_violation = c("base", "strat_miller", "strat_sinusoid")
scen_sampling_strategy = c("strat_sinusoid", "gap_est", "gap_prior")

plot_tree_stats_model_violation = function(){
  
  p1 = df_stat |>
    filter(!is.na(rf_dist)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = rf_dist, color = analysis)) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(color = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    labs(title = "RF distance MAP - truth",
         x = nchar_axis_label,
         y = "RF Distance")
    
  
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

ggsave(filename = "figs/post_analysis/tree_statistics_model_violation.png",
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

ggsave(filename = "figs/post_analysis/tree_statistics_sampling_strategy.png",
       plot = plot_tree_stats_sampling_strategy())

#### Posterior of continuous characters ####
df = df_median 
df$analysis = factor(df$analysis, levels = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior"), ordered = TRUE)
df$nchar = factor(df$nchar, levels = n_chars, ordered = TRUE)

plot_base_vs_strat = function(){
  p1 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = ext_rel_error, fill = analysis)) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = 0.9),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                                jitter.width = 0.2,
                                                jitter.height = 0.01),
                alpha = 0.6) +
    geom_hline(yintercept = 0) +
    labs(y = "Rel. error extinction",
         x = "# Char")
  p2 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = spec_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error speciation")
  p3 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = origin_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error origin")
  
  p4 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error mol branch rates")
  
  p5 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_rel_error, fill = analysis)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    labs(y = "relative error morpho branch rates")
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, common.legend = TRUE)
  return(p)
}

p = plot_base_vs_strat()
p
ggsave(filename = "figs/post_analysis/param_com_model_violations.png",
       plot = plot_base_vs_strat())


plot_sampling_strategies_comp = function(){
  selection = scen_sampling_strategy
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

ggsave(filename = "figs/post_analysis/param_comp_sampling_strategies.png",
       plot = plot_sampling_strategies_comp())

plot_coverage_freq = function(){
  ylim_min = 0.4
  selection = scen_model_violation
  p1 = df |> select(spec_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(spec_coverage_freq = mean(spec_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = spec_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9) +
    labs(x = nchar_axis_label,
         y = "Cov. freq.",
         title = "Speciation",
         color = "Analysis") +
    ylim(ylim_min, 1) +
    scale_color_discrete(labels = scen_axis_labels)
  
  p2 = df |> select(ext_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(ext_coverage_freq = mean(ext_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = ext_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9)+
    labs(x = nchar_axis_label,
         y = "Cov. freq.",
         title = "Extinction",
         color = "Analysis") +
    ylim(ylim_min, 1) +
    scale_color_discrete(labels = scen_axis_labels)
  
  p3 = df |> select(branch_rates_mol_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_mol_coverage_freq = mean(branch_rates_mol_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = branch_rates_mol_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9) +
    labs(x = nchar_axis_label,
         y = "Cov. freq.",
         title = "Mol. branch rate",
         color = "Analysis") +
    ylim(ylim_min, 1) +
    scale_color_discrete(labels = scen_axis_labels)
  
  p4 = df |> select(branch_rates_morpho_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_morpho_coverage_freq = mean(branch_rates_morpho_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_coverage_freq, color = analysis)) +
    geom_point() +
    geom_hline(yintercept = 0.9) +
    labs(x = nchar_axis_label,
         y = "Cov. freq.",
         title = "Morph. branch rate",
         color = "Analysis") +
    ylim(ylim_min, 1) +
    scale_color_discrete(labels = scen_axis_labels)
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "bottom")
  return(p)
  
}

ggsave(filename = "figs/post_analysis/coverage_freq.png",
       plot = plot_coverage_freq())


# doodles for comparing the ranges
d = df |>
  select(analysis, id, nchar, speciation_rate_5., speciation_rate_95., speciation_rate_50.) |>
  mutate(across(c(speciation_rate_5., speciation_rate_95., speciation_rate_50.), as.numeric)) |>
  mutate(span = speciation_rate_95. - speciation_rate_5.) |>
  group_by(analysis, nchar) |>
  mutate(order = row_number(span)) |>
  ungroup() |>
  group_by(analysis, nchar) |>
  arrange(order, .by_group = TRUE) |>
  ungroup() |>
  mutate(x = row_number())

d  |>
  ggplot(aes(x = x, y = speciation_rate_50. ,color = nchar, group = nchar)) +
  geom_pointrange(aes(
    ymin = speciation_rate_5.,
    ymax = speciation_rate_95.),
    position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = lambda)


dd = df |>
  select(analysis, id, nchar, speciation_rate_50., extinction_rate_50.) |>
  mutate(across(c(speciation_rate_50., extinction_rate_50.), as.numeric)) |>
  mutate(across(c(nchar, analysis), as.factor)) |>
  ggplot(aes(x = speciation_rate_50., y = extinction_rate_50., color = analysis, group = analysis)) +
  geom_point() +
  annotate("point", x = lambda, y = mu, size = 5, color = "black")

dd

