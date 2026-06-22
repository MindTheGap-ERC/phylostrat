source("code/constants.R")
load("data/post_values.RData")
load("data/convergence_assessment.RData")

library(ggplot2)
library(dplyr)
library(scales)
library(admtools)

fig_width_2col_cm = 18
fig_height_2_col_max_cm = 17


scen_axis_labels = c("base" = "cFBD +\nCFS",
                     "strat_miller" = "cFBD +\nMSG",
                     "strat_sinusoid" = "cFBD +\nFLG",
                     "gap_est" = "sFBDw +\nFLG",
                     "gap_prior" = "sFBDs +\nFLG")
nchar_axis_label = "# Morph. char"
scen_model_violation = c("base", "strat_miller", "strat_sinusoid")
scen_sampling_strategy = c("strat_sinusoid", "gap_est", "gap_prior")

plot_n_converged_runs = function(){
  p = df_converged |>
    group_by(analysis, nchars) |>
    summarise(n_converged = sum(converged), .groups = "drop") |>
    mutate(nchars = as.factor(nchars)) |>
    ggplot(aes(x = analysis, y = n_converged, color = nchars, shape = nchars)) +
    geom_point(position = position_jitter(height = 0, width = 0.1)) +
    scale_y_continuous(
      breaks = breaks_width(1),
      labels = as.integer,
      limits = c(30, 40)
    ) +
    labs(title = "Converged runs",
         x = "Analysis",
         y = "Converged (out of 40)",
         color = nchar_axis_label,
         shape = nchar_axis_label) +
    scale_x_discrete(labels = scen_axis_labels) +
    theme(legend.position = "bottom")
  return(p)
}
ggsave(filename = "figs/sm/converged_runs.png",
       plot = plot_n_converged_runs(),
       bg = "white",
       width = fig_width_2col_cm,
       height = 10,
       units = "cm")


#### Tree figures ####
df_stat$analysis = factor(df_stat$analysis,
                          levels = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior"),
                          ordered = TRUE)

plot_tree_stats_model_violation = function(){
  dodge = 0.6
  p1 = df_stat |>
    filter(!is.na(rf_dist)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = rf_dist, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(color = analysis, fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    labs(title = "RF distance MAP - truth",
         x = nchar_axis_label,
         y = "RF Distance [-]",
         fill = "Analysis") +
    guides(color = "none")  +
    ylim(0, 1)
    
  
  p2 = df_stat |>
    filter(!is.na(mean_prec_div_times)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = mean_prec_div_times, color = analysis, fill = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    labs(title = "Divergence times",
         x = nchar_axis_label,
         y = "Mean precision [-]",
         fill = "Analysis") +
    guides(color = "none") +
    ylim(0,3)
  
  p3 = df_stat |>
    filter(!is.na(cov_freq_counted) & !is.nan(cov_freq_counted)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = cov_freq_counted, color = analysis, fill = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    labs(title = "Divergence times",
         x = nchar_axis_label,
         y = "95% HDP coverage freq. [-]",
         fill = "Analysis") +
    guides(color = "none") +
    ylim(0,1) # +
    # geom_hline(yintercept = 0.95,
    #            linetype = "dashed")
  
  p4 = df_stat |>
    filter(!is.na(mean_sa_cov_freq)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = mean_sa_cov_freq, color = analysis, fill = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    ylim(0,1) + 
    labs(title = "SA identification",
         y = "95% id. freq. [-]",
         x = nchar_axis_label,
         fill = "Analysis") +
    guides(color = "none") #+
    # geom_hline(yintercept = 0.95,
    #            linetype = "dashed")
  
  p = ggpubr::ggarrange(p1, p2, p3, p4, 
                        ncol = 2,
                        nrow = 2,
                        common.legend = TRUE,
                        labels = LETTERS[1:4],
                        legend = "bottom")
}

ggsave(filename = "figs/ms/tree_statistics_model_violation.png",
       plot = plot_tree_stats_model_violation(),
       bg = "white",
       width = fig_width_2col_cm,
       height = 15,
       units = "cm")

plot_tree_stats_sampling_strategy = function(){
  dodge = 0.6
  p1 = df_stat |>
    filter(!is.na(rf_dist)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = rf_dist, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(color = analysis, fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    labs(title = "RF distance MAP - truth",
         x = nchar_axis_label,
         y = "RF Distance [-]",
         fill = "Analysis")+
    guides(color = "none")  +
    ylim(0, 1)
  
  
  p2 = df_stat |>
    filter(!is.na(mean_prec_div_times)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = mean_prec_div_times, color = analysis, fill = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    labs(title = "Divergence times",
         x = nchar_axis_label,
         y = "Mean precision [-]",
         fill = "Analysis") +
    guides(color = "none") +
    ylim(0,3)
  
  p3 = df_stat |>
    filter(!is.na(cov_freq_counted) & !is.nan(cov_freq_counted)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = cov_freq_counted, color = analysis, fill = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    labs(title = "Divergence times",
         x = nchar_axis_label,
         y = "95% HDP coverage freq. [-]",
         fill = "Analysis") +
    guides(color = "none") +
    ylim(0,1) #+
    # geom_hline(yintercept = 0.95,
    #            linetype = "dashed")
    # 
  p4 = df_stat |>
    filter(!is.na(mean_sa_cov_freq)) |>
    mutate(across(c(nchar, analysis), factor)) |>
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = mean_sa_cov_freq, color = analysis, fill = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    ylim(0,1) + 
    labs(title = "SA identification",
         y = "95% id. freq. [-]",
         x = nchar_axis_label,
         fill = "Analysis") +
    guides(color = "none")  #+
    # geom_hline(yintercept = 0.95,
    #       linetype = "dashed")
    # 
  p = ggpubr::ggarrange(p1, p2, p3, p4, 
                        ncol = 2,
                        nrow = 2,
                        common.legend = TRUE,
                        labels = LETTERS[1:4],
                        legend = "bottom")
}

ggsave(filename = "figs/ms/tree_statistics_sampling_strategy.png",
       plot = plot_tree_stats_sampling_strategy(),
       bg = "white",
       width = fig_width_2col_cm,
       height = 15,
       units = "cm")

#### Posterior of continuous characters ####
df = df_median 
df$analysis = factor(df$analysis,
                     levels = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior"),
                     ordered = TRUE)
df$nchar = factor(df$nchar,
                  levels = n_chars,
                  ordered = TRUE)

plot_base_vs_strat = function(){
  dodge = 0.6
  p1 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = ext_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Extinction rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p2 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = spec_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Speciation rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p3 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = origin_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Origin",
         fill = "Analysis") +
    guides(color = "none")
  
  p4 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Molecular clock rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p5 = df |> 
    filter(analysis %in% scen_model_violation) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_model_violation]) +
    scale_color_discrete(labels = scen_axis_labels[scen_model_violation]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Morphological clock rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p = ggpubr::ggarrange(p2, p1, p5, p4, p3,
                        ncol = 2,
                        nrow = 3,
                        common.legend = TRUE,
                        legend = "bottom",
                        labels = LETTERS[1:5])
  return(p)
}

p = plot_base_vs_strat()
p
ggsave(filename = "figs/ms/param_com_model_violations.png",
       plot = plot_base_vs_strat(),
       bg = "white",
       width = fig_width_2col_cm,
       height = fig_height_2_col_max_cm,
       units = "cm")


plot_sampling_strategies_comp = function(){
  dodge = 0.6
  p1 = df |> 
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = ext_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Extinction rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p2 = df |> 
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = spec_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Speciation rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p3 = df |> 
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = origin_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Origin",
         fill = "Analysis") +
    guides(color = "none") +
    ylim(0, 1)
  
  p4 = df |> 
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = branch_rates_mol_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Molecular clock rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p5 = df |> 
    filter(analysis %in% scen_sampling_strategy) |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_rel_error, fill = analysis, color = analysis)) +
    scale_fill_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    scale_color_discrete(labels = scen_axis_labels[scen_sampling_strategy]) +
    geom_violin(aes(fill = analysis, color = analysis),
                position = position_dodge(width = dodge),
                alpha = 0.4) +
    geom_jitter(position = position_jitterdodge(dodge.width = dodge,
                                                jitter.width = 0.2,
                                                jitter.height = 0.0),
                alpha = 0.6) +
    geom_hline(yintercept = 0,
               linetype = "dashed") +
    labs(y = "Rel. error [-]",
         x = nchar_axis_label,
         title = "Morphological clock rate",
         fill = "Analysis") +
    guides(color = "none")
  
  p = ggpubr::ggarrange(p2, p1, p5, p4, p3,
                        ncol = 2,
                        nrow = 3,
                        common.legend = TRUE,
                        legend = "bottom",
                        labels = LETTERS[1:5])
  return(p)
}
p = plot_sampling_strategies_comp()
p

ggsave(filename = "figs/ms/param_comp_sampling_strategies.png",
       plot = plot_sampling_strategies_comp(),
       bg = "white",
       width = fig_width_2col_cm,
       height = fig_height_2_col_max_cm,
       units = "cm")

plot_coverage = function(){
  ylim_min = 0.5
  p1 = df |> 
    select(spec_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(spec_coverage_freq = mean(spec_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = spec_coverage_freq, color = analysis, shape = analysis)) +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    labs(x = nchar_axis_label,
         y = "Coverage",
         title = "Speciation rate",
         color = "Analysis",
         shape = "Analysis") +
    ylim(ylim_min, 1) +
    scale_color_discrete(labels = scen_axis_labels) +
    scale_shape_discrete(labels = scen_axis_labels)
  
  p2 = df |> 
    select(ext_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(ext_coverage_freq = mean(ext_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = ext_coverage_freq, color = analysis, shape = analysis)) +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    geom_hline(yintercept = 0.9, linetype = "dashed")+
    labs(x = nchar_axis_label,
         y = "Coverage",
         title = "Extinction rate",
         color = "Analysis") +
    ylim(ylim_min, 1)  +
    scale_color_discrete(labels = scen_axis_labels) +
    scale_shape_discrete(labels = scen_axis_labels)
  
  p3 = df |> 
    select(branch_rates_mol_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_mol_coverage_freq = mean(branch_rates_mol_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = branch_rates_mol_coverage_freq, color = analysis, shape = analysis)) +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    labs(x = nchar_axis_label,
         y = "Coverage",
         title = "Mol. clock rate",
         color = "Analysis",
         shape = "Analysis") +
    ylim(ylim_min, 1) +
    scale_color_discrete(labels = scen_axis_labels) +
    scale_shape_discrete(labels = scen_axis_labels)
  
  p4 = df |> 
    select(branch_rates_morpho_covered, analysis, nchar) |>
    group_by(analysis, nchar) |>
    summarise(branch_rates_morpho_coverage_freq = mean(branch_rates_morpho_covered), .groups = "drop") |>
    ggplot(aes(x = nchar, y = branch_rates_morpho_coverage_freq, color = analysis, shape = analysis)) +
    geom_point(position = position_jitter(width = 0.05, height = 0)) +
    geom_hline(yintercept = 0.9, linetype = "dashed") +
    labs(x = nchar_axis_label,
         y = "Coverage",
         title = "Morph. clock rate",
         color = "Analysis",
         shape = "Analysis") +
    ylim(ylim_min, 1)  +
    scale_color_discrete(labels = scen_axis_labels) +
    scale_shape_discrete(labels = scen_axis_labels)
  
  p = ggpubr::ggarrange(p1, p2, p3, p4,
                        common.legend = TRUE,
                        legend = "bottom",
                        labels = LETTERS[1:4])
  return(p)
  
}

ggsave(filename = "figs/sm/coverage.png",
       plot = plot_coverage(),
       bg = "white",
       width = fig_height_2_col_max_cm,
       height = 14,
       units = "cm")

#### Strat Figure ####
df_strat = read.csv("data/strat/selected_adms.csv")

adm_miller = admtools::tp_to_adm(t = df_strat$t,
                                 h = df_strat$miller,
                                 T_unit = "Myr",
                                 L_unit = "m")
adm_sinusoid = admtools::tp_to_adm(t = df_strat$t,
                                   h = df_strat$sinusoid,
                                   T_unit = "Myr",
                                   L_unit = "m")

plot_strat_fig = function(){
  t_mod = df_strat$t[seq(1, length(df_strat$t), by = 1)]
  
  df2 = data.frame(t = rep(t_mod, 2),
                   h = c(time_to_strat(t_mod, adm_miller), time_to_strat(t_mod, adm_sinusoid)),
                   sl = c(rep("Miller et al.", length(t_mod)), rep("Sinusoid", length(t_mod))))
  p1 = df2 |> 
    ggplot(aes(x = t, y = h, color = sl)) +
    geom_line(linewidth = 1.5) +
    ggtitle("Age-depth models") +
    labs(x = "Elapsed simulation time [Myr]",
         y = "Stratigraphic height [m]",
         color = "Sea level curve") +
    theme(legend.position = c(0.2 ,0.8))
  
  df3 = data.frame(dur = c(get_hiat_duration(adm_miller),
                           get_hiat_duration(adm_sinusoid)),
                   sl = factor(c( rep("Miller et al.", get_hiat_no(adm_miller)),
                                  rep("Sinusoid", get_hiat_no(adm_sinusoid)))))
  p2 = df3 |> ggplot(aes(x = dur, fill = sl)) +
    geom_histogram(alpha=0.6, position = "identity", color = "#e9ecef", bins = 20) +
    scale_x_log10() +
    labs(x = "Gap duration [Myr]",
         y = "No. of gaps",
         fill = "Sea level curve") +
    theme(legend.position = c(0.8, 0.8)) +
    ggtitle("Gap duration")
  
  p = ggpubr::ggarrange(p1, p2, 
                        ncol = 2, 
                        nrow = 1, 
                        labels = LETTERS[1:2])
  
  return(p)
}

p = plot_strat_fig()
ggsave(filename = "figs/ms/strat_info.png",
       plot = plot_strat_fig(),
       bg = "white",
       width = fig_width_2col_cm,
       height = 7,
       units = "cm")

#### Precision vs. accuracy ####
df_median$ext_hdi_rel = (as.numeric(df_median$extinction_rate_95.) - as.numeric(df_median$extinction_rate_5.))/mu
df_median$spec_hdi_rel = (as.numeric(df_median$speciation_rate_95.) - as.numeric(df_median$speciation_rate_5.))/lambda

df_median$analysis = factor(df_median$analysis)
df_median$nchar = factor(df_median$nchar)


plot_acc_vs_prec_ext = function(){
  plot_comp = function(nchars){
    p = df_median |> 
      filter(nchar == nchars) |>
      ggplot(aes(x = ext_rel_error,
                 y = ext_hdi_rel,
                 shape = ext_covered, 
                 color = analysis)) +
      geom_point() +
      scale_color_discrete(labels = scen_axis_labels,
                           name = "Case")  +
      labs(title = paste(nchars , "Morph. Char."),
           x = "Rel. error ext. rate [-]",
           y = "Rel. HDI width ext. rate [-]") +
      guides(shape = "none") +
      ylim(range(df_median$ext_hdi_rel)) +
      xlim(range(df_median$ext_rel_error))
    return(p)
  }

  p = ggpubr::ggarrange(plot_comp("30"),
                        plot_comp("100"),
                        plot_comp("300"),
                        plot_comp("1000"),
                        common.legend = TRUE,
                        labels = LETTERS[1:4],
                        legend = "bottom")
  return(p)
}

plot_acc_vs_prec_spec = function(){
  plot_comp = function(nchars){
    p = df_median |> 
      filter(nchar == nchars) |>
      ggplot(aes(x = spec_rel_error,
                 y = spec_hdi_rel,
                 shape = spec_covered, 
                 color = analysis)) +
      geom_point() +
      scale_color_discrete(labels = scen_axis_labels,
                           name = "Case")  +
      labs(title = paste(nchars , "Morph. Char."),
           x = "Rel. error spec. rate [-]",
           y = "Rel. HDI width spec. rate [-]") +
      guides(shape = "none") +
      ylim(range(df_median$spec_hdi_rel)) +
      xlim(range(df_median$spec_rel_error))
    return(p)
  }
  
  p = ggpubr::ggarrange(plot_comp("30"),
                        plot_comp("100"),
                        plot_comp("300"),
                        plot_comp("1000"),
                        common.legend = TRUE,
                        labels = LETTERS[1:4],
                        legend = "bottom")
  return(p)
}


ggsave(filename = "figs/sm/ext_acc_vs_prec.png",
       plot = plot_acc_vs_prec_ext(),
       bg = "white",
       width = fig_width_2col_cm,
       height = fig_width_2col_cm,
       units = "cm")

ggsave(filename = "figs/sm/spec_acc_vs_prec.png",
       plot = plot_acc_vs_prec_spec(),
       bg = "white",
       width = fig_width_2col_cm,
       height = fig_width_2col_cm,
       units = "cm")
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

