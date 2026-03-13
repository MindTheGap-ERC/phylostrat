library(RevGadgets)
library(coda)

nchar = "30"
run= "1"
rep = "1"
folder = "fbd_base"
case = ""

path = paste0("data/", folder, "/rb_output/num_",rep,"_nchar", nchar, case, "_run_", run, ".log")

tr = RevGadgets::readTrace(path = c(path)) |>
  as.mcmc()
aa = effectiveSize(tr[[1]]) |> names()

nchars = c("30", "100", "300", "1000")
reps = 1:40 |> as.character()
runs = 1:2 |> as.character()

names = c(aa[-1], "converged", "run", "rep", "nchar", "analysis", "sampling")
df = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df) = names
missing_files = c()

for (nchar in nchars){
  for (rep in reps){
    for (run in runs){
      for (sampling in c("constant", "long_gap", "short_gap")){
        if (sampling == "constant"){
          folder = "fbd_base"
          case = ""
        }
        if (sampling == "long_gap"){
          folder = "fbd_strat"
          case = "_sinusoid"
        }
        if (sampling == "short_gap"){
          folder = "fbd_strat"
          case = "_miller"
        }
        path = paste0("data/", folder, "/rb_output/num_",rep,"_nchar", nchar, case, "_run_", run, ".log")
        if (file.exists(path)){
          tr = RevGadgets::readTrace(path = c(path))
          tr_quant = as.mcmc(tr)
          aa = effectiveSize(tr_quant[[1]])
          aa = aa[-1]
          aa = c(aa, "converged" = all(aa>200),
                 "run" = run,
                 "rep" = rep, 
                 "nchar" = nchar, 
                 "analysis" = "base",
                 "path" = path,
                 "sampling" = sampling,
                 "priors" = "none")
          df2 = as.data.frame(as.list(aa))
          df = rbind(df, df2)
        }
        if (! file.exists(path)){
          missing_files = c(missing_files, path)
        }
      }
    }
  }
}


# gappy inference
folder = "fbd_gap_est"
case = "_sinusoid"

path = paste0("data/", folder, "/rb_output/num_",rep,"_nchar", nchar, case, "_run_", run, ".log")

tr = RevGadgets::readTrace(path = c(path)) |>
  as.mcmc()
aa = effectiveSize(tr[[1]]) |> names()
names = c(aa[-1], "converged", "run", "rep", "nchar", "analysis", "priors")
df3 = data.frame(matrix(nrow = 0, ncol = length(names)))
names(df3) = names

for (nchar in nchars){
  for (rep in reps){
    for (run in runs){
      for (priors in c("strong", "weak")){
        if (priors  == "strong"){
          folder = "fbd_gap_prior"
          case = "_sinusoid"
        }
        if (priors == "weak"){
          folder = "fbd_gap_est"
          case = "_sinusoid"
        }
        path = paste0("data/", folder, "/rb_output/num_",rep,"_nchar", nchar, case, "_run_", run, ".log")
        if (file.exists(path)){
          tr = RevGadgets::readTrace(path = c(path))
          tr_quant = as.mcmc(tr)
          aa = effectiveSize(tr_quant[[1]])
          aa = aa[-1]
          if (priors == "strong"){
            aa_sel = aa[! names(aa) %in% c("psi[2]", "psi[4]")]
          }
          if (priors == "weak"){
            aa_sel = aa
          }
          aa = c(aa, "converged" = all(aa_sel>200),
                 "run" = run,
                 "rep" = rep, 
                 "nchar" = nchar, 
                 "analysis" = "base",
                 "path" = path,
                 "priors" = priors,
                 "sampling" = "long_gap")
          df2 = as.data.frame(as.list(aa))
          df3 = rbind(df3, df2)
        }
        if (! file.exists(path)){
          missing_files = c(missing_files, path)
        }
      }
    }
  }
}

df$ID = seq_along(df$Posterior)
df$ID
df3$ID = max(df$ID) + seq_along(df3$Posterior)
df_total = merge(df, df3, by = 'ID', all = TRUE)
