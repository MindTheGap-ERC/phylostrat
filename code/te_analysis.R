library(ggplot2)
library(RevGadgets)

if (!dir.exists("figs/te")) {
  dir.create("figs/te", recursive = TRUE)
}

## constants
nchars = c(100, 300, 1000)
cases = c("cont", "inc_A", "inc_B")
ids = as.character(1:10)
path = "test/te_phylo/"
app1 = "skyline_A"
app2 = "skyline_A_full"
true_ex = 0.3
true_or = 1

## standard scenarios

colnames = c("id", "nchar","case", "ex_mean", "or_mean")
df_cont = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df_cont) = colnames
for (id in ids){
  for (nchar in nchars){
    for (case in cases){
      cat(paste0(id, nchar, case))
      tr1 = paste0(path, "res_const_te_", case, "_" ,id, "_", nchar, "_run_1.log")
      tr2 = paste0(path, "res_const_te_", case, "_" ,id, "_", nchar, "_run_2.log")
      
      if(!file.exists(tr1) | !file.exists(tr2)){
        cat("skipping file")
        next
      }
      
      comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
        RevGadgets::combineTraces()
      
      ex_mean = comb_trace[[1]]$extinction_rate |> mean()
      or_mean = comb_trace[[1]]$speciation_rate |> mean()
      rm(comb_trace)
      df2 = data.frame(id = id,
                       nchar = nchar,
                       ex_mean = ex_mean,
                       or_mean = or_mean,
                       case = case)
      df_cont = rbind(df_cont, df2)
      print(df_cont)
      
    }
  }
}
df_cont$ex_re = (df_cont$ex_mean - true_ex)/true_ex
df_cont$or_re = (df_cont$or_mean - true_or)/true_or
df_cont$nchar = factor(df_cont$nchar, levels = nchars)
df_cont$case = factor(df_cont$case, levels = cases)

save(df_cont, file = "data/df_cont.RData")

load("data/df_cont.RData")

# plot relative error by nchar and case
p_ex = df_cont |>
  ggplot(aes(x = nchar, y = ex_re, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.5))
p_ex
p_or = df_cont |>
  ggplot(aes(x = nchar, y = or_re, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.5))
p_or
# save plot
ggsave("figs/te/or_cont.png", plot = p_or)
ggsave("figs/te/ex_cont.png", plot = p_ex)

#### skyline ####
colnames = c("id", "nchar","ex_mean", "or_mean")
df_skyline = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df_skyline) = colnames
for (id in ids){
  for (nchar in nchars){
    cat(paste0(id, nchar))
    tr1 = paste0(path, "res_", app1, "_te_" ,id, "_", nchar, "_run_1.log")
    tr2 = paste0(path, "res_", app1, "_te_" ,id, "_", nchar, "_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat("skipping file\n")
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ex_mean = comb_trace[[1]]$extinction_rate |> mean()
    or_mean = comb_trace[[1]]$speciation_rate |> mean()
    rm(comb_trace)
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ex_mean = ex_mean,
                     or_mean = or_mean)
    df_skyline = rbind(df_skyline, df2)
    print(df_skyline)
  }
}

df_skyline$ex_re = (df_skyline$ex_mean - true_ex)/true_ex
df_skyline$or_re = (df_skyline$or_mean - true_or)/true_or
df_skyline$nchar = factor(df_skyline$nchar, levels = nchars)

p_skyline_ex = df_skyline |>
  ggplot(aes(x = nchar, y = ex_re)) +
  geom_boxplot(position = position_dodge(width = 0.5))
p_skyline_ex
p_skyline_or = df_skyline |>
  ggplot(aes(x = nchar, y = or_re)) +
  geom_boxplot(position = position_dodge(width = 0.5))
p_skyline_or


ggsave("figs/te/or_skyline.png", plot = p_skyline_or)
ggsave("figs/te/ex_skyline.png", plot = p_skyline_ex)

save(df_skyline, file = "data/df_skyline.RData")
load("data/df_skyline.RData")


## skyline full
colnames = c("id", "nchar","ex_mean", "or_mean", "psi1", "psi2", "psi3", "psi4", "psi5")
df_skyline_full = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df_skyline_full) = colnames
for (id in ids){
  for (nchar in nchars){
    cat(paste0(id, nchar))
    tr1 = paste0(path, "res_", app2, "_te_" ,id, "_", nchar, "_run_1.log")
    tr2 = paste0(path, "res_", app2, "_te_" ,id, "_", nchar, "_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat("skipping file\n")
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ex_mean = comb_trace[[1]]$extinction_rate |> mean()
    or_mean = comb_trace[[1]]$speciation_rate |> mean()
    psi1 = comb_trace[[1]]$`psi[1]` |> mean()
    psi2 = comb_trace[[1]]$`psi[2]`  |> mean()
    psi3 = comb_trace[[1]]$`psi[3]`  |> mean()
    psi4 = comb_trace[[1]]$`psi[4]`  |> mean()
    psi5 = comb_trace[[1]]$`psi[5]`  |> mean()
    rm(comb_trace)
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ex_mean = ex_mean,
                     or_mean = or_mean,
                     psi1 = psi1,
                     psi2 = psi2,
                     psi3 = psi3,
                     psi4 = psi4,
                     psi5 = psi5)
    df_skyline_full = rbind(df_skyline_full, df2)
    print(df_skyline_full)
  }
}

df_skyline_full$ex_re = (df_skyline_full$ex_mean - true_ex)/true_ex
df_skyline_full$or_re = (df_skyline_full$or_mean - true_or)/true_or
df_skyline_full$nchar = factor(df_skyline_full$nchar, levels = nchars)

p_skyline_full_ex = df_skyline_full |>
  ggplot(aes(x = nchar, y = ex_re)) +
  geom_boxplot(position = position_dodge(width = 0.5))
p_skyline_full_ex
p_skyline_full_or = df_skyline_full |>
  ggplot(aes(x = nchar, y = or_re)) +
  geom_boxplot(position = position_dodge(width = 0.5))
p_skyline_full_or

ggsave("figs/te/or_skyline_full.png", plot = p_skyline_full_or)
ggsave("figs/te/ex_skyline_full.png", plot = p_skyline_full_ex)

save(df_skyline_full, file = "data/df_skyline_full.RData")

## combine data sets
load("data/df_skyline_full.RData")
load("data/df_skyline.RData")
load("data/df_cont.RData")
library(ggplot2)
df = df_cont

df2 = df_skyline
df2$case = factor(rep("skyline_A", nrow(df2)))
df = rbind(df, df2)

df3 = df_skyline_full
df3 = df3[ , c("id", "nchar", "ex_mean", "or_mean", "ex_re", "or_re")]
df3$case = factor(rep("skyline_A_full", nrow(df3)))

df = rbind(df, df3)
p1 = df |>
  ggplot(aes(x = nchar, y = ex_re, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  ylab("Relative error extinction rate")
p1
p2 = df |>
  ggplot(aes(x = nchar, y = or_re, fill = case)) +
  geom_boxplot(position = position_dodge(width = 0.5)) +
  ylab("Relative error speciation rate")
p2

ggsave("figs/te/or_combined.png", plot = p2)
ggsave("figs/te/ex_combined.png", plot = p1)


df4 = data.frame(psi = c(df_skyline_full$psi1, df_skyline_full$psi2, df_skyline_full$psi3, df_skyline_full$psi4, df_skyline_full$psi5),
                 int =  c(rep("1", length(df_skyline_full$psi1)),
                          rep("2", length(df_skyline_full$psi1)),
                          rep("3", length(df_skyline_full$psi1)),
                          rep("4", length(df_skyline_full$psi1)),
                          rep("5", length(df_skyline_full$psi1))))

p3 = df4 |> ggplot(aes(x = int, y = psi)) +
  geom_boxplot() +
  xlab("Interval")

ggsave("figs/te/psi_intervals.png", plot = p3)
