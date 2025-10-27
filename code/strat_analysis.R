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




