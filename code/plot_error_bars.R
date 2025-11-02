library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)

if (!dir.exists("figs/")) {
  dir.create("figs/", recursive = TRUE)
}

path_res = "output/fbd_base/"

source("code/constants.R")

# constants

n_rep = 20
ids = as.character(seq_len(n_rep))

alpha = 0.25
kappa = 5

case = "base"
# create empty dataframe for storage
colnames = c("id"  ,            "nchar"  ,         "ext_min"     ,    "ext_max"    ,     "spec_min" ,      
             "spec_max"    ,    "br_rate_min"  ,   "br_rate_max"   ,  "br_rate_mol_min", "br_rate_mol_max",
              "alpha_min"  ,     "alpha_max"    ,   "kappa_min"   ,    "kappa_max"  ,     "case"  )
df = data.frame(matrix(nrow = 0, ncol = length(colnames)))
colnames(df) = colnames



norm_data = function(x, true_val){
  quantile(x - true_val, probs = c(0.05, 0.95)) |> unname() 
}
list_base = list()
# read in data
for (id in ids){
  for (nchar in n_chars){
    tr1 = paste0(path_res,"num_", id, "_nchar", nchar,"_run_1.log")
    tr2 = paste0(path_res, "num_", id, "_nchar", nchar,"_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat(paste0("skipping file ", tr1, "\n"))
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ext_est = comb_trace[[1]]$extinction_rate |> norm_data(true_val = mu)
    spec_est = comb_trace[[1]]$speciation_rate |> norm_data(true_val = lambda)
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> norm_data(true_val = clock_rate_morph)
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> norm_data(true_val = clock_rate_mol)
    kappa_est = comb_trace[[1]]$kappa |> norm_data(true_val = kappa)
    alpha_est = comb_trace[[1]]$alpha_mol |> norm_data(true_val = alpha)

    
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ext_min = ext_est |> min(),
                     ext_max = ext_est |> max(),
                     spec_min = spec_est |> min(),
                     spec_max = spec_est |> max(),
                     br_rate_min = br_rate_est |> min(),
                     br_rate_max = br_rate_est |> max(),
                     br_rate_mol_min = br_rate_mol_est |> min(),
                     br_rate_mol_max = br_rate_mol_est |> max(),
                     alpha_min = alpha_est|> min(),
                     alpha_max = alpha_est|> max(),
                     kappa_min = kappa_est|> min(),
                     kappa_max = kappa_est|> max(),
                     case = "base")
    df = rbind(df, df2)
  }
}



case = "miller"
path_res = "output/fbd_strat/"

for (id in ids){
  for (nchar in n_chars){
    tr1 = paste0(path_res,"num_", id, "_nchar", nchar, "_", case, "_run_1.log")
    tr2 = paste0(path_res, "num_", id, "_nchar", nchar,"_", case, "_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat(paste0("skipping file ", tr1, "\n"))
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ext_est = comb_trace[[1]]$extinction_rate |> norm_data(true_val = mu)
    spec_est = comb_trace[[1]]$speciation_rate |> norm_data(true_val = lambda)
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> norm_data(true_val = clock_rate_morph)
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> norm_data(true_val = clock_rate_mol)
    kappa_est = comb_trace[[1]]$kappa |> norm_data(true_val = kappa)
    alpha_est = comb_trace[[1]]$alpha_mol |> norm_data(true_val = alpha)
    
    
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ext_min = ext_est |> min(),
                     ext_max = ext_est |> max(),
                     spec_min = spec_est |> min(),
                     spec_max = spec_est |> max(),
                     br_rate_min = br_rate_est |> min(),
                     br_rate_max = br_rate_est |> max(),
                     br_rate_mol_min = br_rate_mol_est |> min(),
                     br_rate_mol_max = br_rate_mol_est |> max(),
                     alpha_min = alpha_est|> min(),
                     alpha_max = alpha_est|> max(),
                     kappa_min = kappa_est|> min(),
                     kappa_max = kappa_est|> max(),
                     case = "miller")
    df = rbind(df, df2)
  }
}


case = "sinusoid"


for (id in ids){
  for (nchar in n_chars){
    tr1 = paste0(path_res,"num_", id, "_nchar", nchar, "_", case, "_run_1.log")
    tr2 = paste0(path_res, "num_", id, "_nchar", nchar,"_", case, "_run_2.log")
    
    if(!file.exists(tr1) | !file.exists(tr2)){
      cat(paste0("skipping file ", tr1, "\n"))
      next
    }
    
    comb_trace = RevGadgets::readTrace(path = c(tr1, tr2), burnin = 0) |> 
      RevGadgets::combineTraces()
    
    ext_est = comb_trace[[1]]$extinction_rate |> norm_data(true_val = mu)
    spec_est = comb_trace[[1]]$speciation_rate |> norm_data(true_val = lambda)
    br_rate_est = comb_trace[[1]]$branch_rates_morpho |> norm_data(true_val = clock_rate_morph)
    br_rate_mol_est = comb_trace[[1]]$branch_rates_mol |> norm_data(true_val = clock_rate_mol)
    kappa_est = comb_trace[[1]]$kappa |> norm_data(true_val = kappa)
    alpha_est = comb_trace[[1]]$alpha_mol |> norm_data(true_val = alpha)
    
    
    df2 = data.frame(id = id,
                     nchar = nchar,
                     ext_min = ext_est |> min(),
                     ext_max = ext_est |> max(),
                     spec_min = spec_est |> min(),
                     spec_max = spec_est |> max(),
                     br_rate_min = br_rate_est |> min(),
                     br_rate_max = br_rate_est |> max(),
                     br_rate_mol_min = br_rate_mol_est|> min(),
                     br_rate_mol_max = br_rate_mol_est |> max(),
                     alpha_min = alpha_est|> min(),
                     alpha_max = alpha_est|> max(),
                     kappa_min = kappa_est|> min(),
                     kappa_max = kappa_est|> max(),
                     case = "sinusoid")
    
    df = rbind(df, df2)
  }
}


cols = c("base" = "black", "miller" = "blue", "sinusoid" = "red")
lwds = c("30" = 1, "100" = 2, "300" = 3, "1000" = 4)


plot_error_bars = function(param){
  y_min = min(df[,paste0(param, "_min")])
  y_max = max(df[,paste0(param, "_max")])
  plot(NULL, ylim = c(y_min, y_max), xlim = c(0, length(df$id)),
       ylab = paste0(param, "normalized to mean"), xlab = NULL,
       main = param)
  i = 0
  for (case in c("base", "miller", "sinusoid")){
    for (nchar in c(30, 100, 300, 1000)){
      mins = df[df$case == case & df$nchar == nchar,paste0(param, "_min")]
      maxs = df[df$case == case & df$nchar == nchar,paste0(param, "_max")]
      for (j in seq_along(mins)){
        lines(x = c(i, i), y = c(mins[j], maxs[j]), col = cols[case], lwd = lwds[as.character(nchar)])
        i = i + 1
      }
    }
  }
  legend("topleft", col = cols, legend = names(cols), lwd = 1)
  legend("topright", lwd = lwds, legend = names(lwds), col = "black")
  
}
poss_names = c("ext", "spec", "br_rate", "br_rate_mol", "alpha", "kappa")
for (name in poss_names){
  png(paste0("figs/error_bars", name, ".png"), width = 1000, height = 1000)
  plot_error_bars(name)
  dev.off()
}

