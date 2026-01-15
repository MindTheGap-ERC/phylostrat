#### Stratigraphic context ####
source("code/constants.R")
miller_data = read.csv("data/strat/miller_2020_adm.csv")
sinusoid_data = read.csv("data/strat/sinusoid_adm.csv")
adm_miller = tp_to_adm(t = miller_data$time..Myr.- min(miller_data$time..Myr.), h = miller_data$adm_1..m.)
adm_sinusoid = tp_to_adm(t = sinusoid_data$time..Myr., h = sinusoid_data$adm_1..m.)

df = data.frame(t = sinusoid_data$time..Myr., sinusoid = sinusoid_data$adm_1..m., miller = miller_data$adm_1..m.)

plot(df$t, df$sinusoid)
plot(df$t, df$miller)

write.csv(df, "data/strat/selected_adms.csv")

plot(adm_miller, lty_destr = 0)
plot(adm_sinusoid, lty_destr = 0)

get_hiat_list(adm_sinusoid)
get_hiat_duration(adm_sinusoid)

h_list = get_hiat_list(adm_sinusoid)
t = c()
for (i in seq_along(h_list) ){
  if (h_list[[i]]["end"] - h_list[[i]]["start"]> 0.5){
    t = c(t, h_list[[i]]["start"], h_list[[i]]["end"])
  }
}

timeline = t_max - rev(t) |> unname()
#timeline

#cat("timeline is", timeline)

wd_data = read.csv("data/strat/sinusoid_wd.csv")
plot(x = wd_data$time..Myr., 
     y = pmax(wd_data$wd_6..m., 0), 
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Water depth [m]",
     main = "Water depth",
     lwd = 4)

library(StratPal)
f = snd_niche(40, 10)
g = approxfun(wd_data$time..Myr., pmax(wd_data$wd_6..m., 0))
t = wd_data$time..Myr.
plot(t, f(g(t)),
     type = "l",
     xlab = "Time [Myr]",
     ylab = "Collection probability",
     main = "Collection probability",
     lwd = 4)

plot(adm_sinusoid, lty_destr = 0, lwd_acc = 4)
mtext("Time [Myr]", side = 1, line = 2.7)
mtext("Stratigraphic Height [m]", side = 2, line = 2.7)
title(main = "Age-depth model 2 km from shore")

avg_sed = max_height(adm_sinusoid)/ max_time(adm_sinusoid)
adm_const = tp_to_adm(t = c(0, max_time(adm_sinusoid)), h = c(0, max_height(adm_sinusoid)))
heights = seq(0, max_height(adm_sinusoid))
plot(heights, strat_to_time(heights, adm_sinusoid) - strat_to_time(heights, adm_const),
     xlab = "Stratigraphic height [m]",
     ylab = "Age error [Myr]",
     type = "l",
     main = "Age Error",
     lwd = 4)
