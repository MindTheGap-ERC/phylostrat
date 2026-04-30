load("converged_runs.RData")

ids = c(1:40)
nchars= c(30, 100, 300, 1000)
analyses = c("base", "strat_miller", "strat_sinusoid", "gap_est", "gap_prior")
burnin = 0

get_median_posterior_vals = function(){
  df_median = data.frame()
  
  for (id in ids){
    for (nchar in nchars){
      for (analysis in analyses){
        if (df_converged$converged[df_converged$nchars == nchar & df_converged$id == id & df_converged$analysis == analysis]){
          if (analysis == "base"){
            path = "data/fbd_base/"
            file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_run_1.log")
            file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_run_2.log")
          }
          if (analysis == "strat_miller"){
            path = "data/fbd_strat/"
            file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_miller_run_1.log")
            file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_miller_run_2.log")
          }
          if (analysis == "strat_sinusoid"){
            path = "data/fbd_strat/"
            file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_1.log")
            file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_2.log")
          }
          if (analysis == "gap_est"){
            path = "data/fbd_gap_est/"
            file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_1.log")
            file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_2.log")
          }
          if (analysis == "gap_prior"){
            path = "data/fbd_gap_prior/"
            file_run1 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_1.log")
            file_run2 = paste0(path,"rb_output/num_", id, "_nchar", nchar,"_sinusoid_run_2.log")
          }
          
          comb_trace = RevGadgets::readTrace(path = c(file_run1, file_run2), burnin = burnin) |> 
            RevGadgets::combineTraces()
          df_temp = as.data.frame(lapply(comb_trace$combined, function(x) quantile(x, p = c(0.05, 0.5, 0.95))))
          df_te = c(t(df_temp))
          names(df_te) = paste(rep(colnames(df_temp), times = nrow(df_temp)),rep(rownames(df_temp), each = ncol(df_temp)), sep = "_")
          df_te = df_te[order(names(df_te))]
          df_te["id"] = id
          df_te["nchar"] = nchar
          df_te["analysis"] = analysis
          df_median = data.table::rbindlist(list(df_median, as.data.frame(as.list(df_te))), fill = TRUE)
        }
      }
    }
  }
  return(df_median)
}

df_median = get_median_posterior_vals()

save(df_median, file = "median_post_values.RData")
