library(admtools)
library(ggplot2)
df = read.csv("data/strat/selected_adms.csv")

adm_miller = admtools::tp_to_adm(t = df$t, h = df$miller, T_unit = "Myr", L_unit = "m")
adm_sinusoid = admtools::tp_to_adm(t = df$t, h = df$sinusoid, T_unit = "Myr", L_unit = "m")


plot_adms = function(){
  df2 = data.frame(t = c(df$t, df$t),
                   h = c(time_to_strat(df$t, adm_miller), time_to_strat(df$t, adm_sinusoid)),
                   sl = c(rep("Miller et al.", length(df$t)), rep("Sinusoid", length(df$t))))
  p = df2 |> 
    ggplot(aes(x = t, y = h, color = sl)) +
    geom_line(linewidth = 1.5) +
    ggtitle("Age-depth models") +
    labs(x = "Elapsed simulation time [Myr]",
         y = "Stratigraphic height [m]",
         color = "Sea level curve") +
    theme(legend.position = c(0.1 ,0.9))
  ggsave("figs/adms.png", plot = p)
  return(p)
}

p1 = plot_adms()

plot_gap_durations = function(){
  df2 = data.frame(dur = c(get_hiat_duration(adm_miller), get_hiat_duration(adm_sinusoid)),
                   sl = factor(c( rep("Miller et al.", get_hiat_no(adm_miller)), rep("Sinusoid", get_hiat_no(adm_sinusoid)))))
  p = df2 |> ggplot(aes(x = dur, fill = sl)) +
    geom_histogram(alpha=0.6, position = "identity", color = "#e9ecef", bins = 20) +
    scale_x_log10() +
    labs(x = "Gap duration [Myr]",
         y = "No. of gaps",
         fill = "Sea level curve") +
    theme(legend.position = c(0.9, 0.9)) +
    ggtitle("Gap duration")
  ggsave("figs/gap_duration.png", plot = p)
  return(p)
}
p2 = plot_gap_durations()
          
p3 = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, labels = LETTERS[1:2])
ggsave("figs/strat_info.png",
       plot = p3)
