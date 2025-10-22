df = read.csv("data/strat/selected_adms.csv")

adm_miller = admtools::tp_to_adm(t = df$t, h = df$miller, T_unit = "Myr", L_unit = "m")
adm_sinusoid = admtools::tp_to_adm(t = df$t, h = df$sinusoid, T_unit = "Myr", L_unit = "m")
