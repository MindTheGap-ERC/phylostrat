#### Stratigraphic context ####
miller_data = read.csv("data/strat/miller_2020_adm.csv")
sinusoid_data = read.csv("data/strat/sinusoid_adm.csv")
adm_miller = tp_to_adm(t = miller_data$time..Myr.- min(miller_data$time..Myr.), h = miller_data$adm_3..m.)
adm_sinusoid = tp_to_adm(t = sinusoid_data$time..Myr., h = sinusoid_data$adm_1..m.)