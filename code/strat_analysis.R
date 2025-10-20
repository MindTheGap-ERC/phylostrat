miller = read.csv("data/strat/miller_2020_adm.csv")
sinusoid = read.csv("data/strat/sinusoid_adm.csv")
plot(miller$time..Myr., miller$adm_2..m.)
library(admtools)
adm = tp_to_adm(t = miller$time..Myr., h = miller$adm_3..m.)

plot(adm, col_destr = NULL )

adm_2 = tp_to_adm(t = sinusoid$time..Myr., h = sinusoid$adm_1..m.)

plot(adm_2, col_destr = NULL)

get_incompleteness(adm)
get_incompleteness(adm_2)


sin_sac = read.csv("data/strat/sinusoid_sc.csv")

plot(sin_sac$time..Myr., sin_sac$sc_1_f1..m.)
