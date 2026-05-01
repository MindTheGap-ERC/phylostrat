library(ggplot2)
library(RevGadgets)
library(phangorn)
library(ggpubr)
library(paleotree)
library(dispRity)
library(treeio)

source("code/constants.R")
id = 12
nchar = 1000
map_path = paste0("output/fbd_base/", "tree_", id, "_nchar", nchar, "_MAP.tre")
truth_path = paste0("data/sim/fbd_base/tree_", id, "_nchar", nchar , ".nex")

truth = ape::read.nexus(truth_path)
map = treeio::read.beast.newick(map_path) 

plot(map@phylo)
plot(truth)

get_divergence_times = function(map, ref){
  tol = 10^-5
  ages = suppressMessages(dateNodes(ref))
  extant_tips = ref$tip.label[ages < tol]
  df = data.frame(matrix(nrow = 0, ncol = 5))
  names(df) = c("tip1", "tip2", "true_age", "min_age", "max_age")
  if(length(extant_tips)< 2){
    warning("need at least 2 extant tips")
    return(df)
  }
  mrca_index = c()
  mrca_list = list()
  #print(extant_tips)
  k = 1
  for (i in 1:(length(extant_tips)-1)){
    for (j in (i+1):length(extant_tips)){
      tips = c(extant_tips[i], extant_tips[j])
      #print(tips)
      mrca = getMRCA(ref, tips)
      if (! mrca %in% mrca_index){
        mrca_index = c(mrca_index, mrca)
        mrca_list[[k]] = tips
        k = k+1
      }
    }
  }
  for (i in seq_along(mrca_index)){
    tips = mrca_list[[i]]
    true_age = ages[mrca_index[i]] |> unname()
    t_hdp = map@data$age_0.95_HPD[map@data$node == getMRCA(map@phylo, tips)][[1]]
    df2 = data.frame(tip1 = tips[1],
                     tip2 = tips[2],
                     true_age = true_age,
                     min_age = min(t_hdp),
                     max_age = max(t_hdp))
    df = rbind(df, df2)
  }
  df$covered = df$true_age >= df$min_age & df$true_age <= df$max_age
  return(df)
}
df = get_divergence_times(map, truth)

get_div_statistics = function(df){
  mean_prec = c((df$max_age - df$min_age)/df$true_age) |> mean()
  cov_freq = mean(df$covered)
  return(c("mean_prec" = mean_prec,
           "cov_freq" = cov_freq))
}

get_sa_prob = function(map, ref){
  tol = 10^-5
  extant_tips = ref$tip.label[suppressMessages(dateNodes(ref)) < tol]
  not_extant = ref$tip.label[! ref$tip.label %in% extant_tips ]
  not_extant_node = which(map@phylo$tip.label %in% not_extant)
  sas = map@data$sampled_ancestor[map@data$node %in% not_extant_node]
  sas[is.na(sas)] = 0
  names(sas) = map@phylo$tip.label[not_extant_node]
  
  sa_true_names = ref$tip.label[ref$edge[ref$edge.length == 0, 2]]
  sa_true = rep(0, length(not_extant))
  names(sa_true) = not_extant
  sa_true[sa_true_names] = 1
  
  return(list(sas, sa_true))
}

get_extant_tip_labels = function(tree, tol = 10^-5){
  return(tree$tip.label[dateNodes(tree) < tol])
}

l = get_sa_prob(map, truth)

all(names(l[[1]]) %in% names(l[[2]]))

mean_sa_cov_freq = function(l) {sapply(names(l[[1]]), function(x) l[[2]][x] - l[[1]][x] < 0.05) |> mean()}
mean_sa_cov_freq(l)

rf_dist = phangorn::RF.dist(map@phylo, truth, rooted = TRUE, normalize = TRUE)

tree_params = function(map, truth){
  df = get_divergence_times(map = map,
                            ref = truth)
  
  x = c("mean_prec_div_times" = c((df$max_age - df$min_age)/df$true_age) |> mean(),
        "cov_freq_div_times" = mean(df$covered),
        "mean_sa_cov_freq" = get_sa_prob(map = map,
                                         ref = truth) |>
          mean_sa_cov_freq(),
        "rf_dist" = phangorn::RF.dist(map@phylo, truth, rooted = TRUE, normalize = TRUE))
  return(x)
}

tree_params(map, truth)
map@data$age_0.95_HPD[[3]]
map@info

map@data$index

tnd = dateNodes(map@phylo)[3]
map@data$age_

plot(truth)
plot(map@phylo)

dateNodes(map@phylo)

a = tree.age(map@phylo)
sort(a$ages)

plot(truth, show.node.label = TRUE) 

plot(drop.tip(truth, "t39_1"))
?drop.tip
truth$edge.length


tree = truth@phylo

ages = dateNodes(tree)




#### Determine correct identification of SAs



extant_tips = get_extant_tip_labels(ref)
not_extant = ref$tip.label[! ref$tip.label %in% extant_tips ]

not_extant_node = which(map@phylo$tip.label %in% not_extant)
sas = map@data$sampled_ancestor[map@data$node %in% not_extant_node]

true_foss = c("t19_1", "t25_1", "t32_3")
map@data$sampled_ancestor[map@data$node %in% which(map@phylo$tip.label %in% true_foss)]








tr = l[[2]][order(names(l[[2]]))]
obs = l[[1]][order(names(l[[1]]))]

sum(abs(tr - obs) <= 0.05) / length(tr)

df$min_adjusted = df$min_age - df$true_age
df$max_adjusted = df$max_age - df$true_age

plot(NULL,
     xlim = c(0, 12),
     ylim = c(-1, 1))
for (i in seq_along(df$tip1)){
  lines(x = c(i,i), y = c(df$min_adjusted[i], df$max_adjusted[i]))
}
lines(c(0, 12), c(0,0))

remove_sa = function(x){
  is_sa = x$edge.length == 0 #& truth$edge[,2] > Ntip(truth)
  sa_nodes = x$edge[is_sa, 2]
  phy_no_sa =  drop.tip(x, sa_nodes)
}

plot(map@phylo |> remove_sa())
plot(truth |> remove_sa())

map@phylo |> remove_sa() |> tree.age()
truth |> remove_sa() |> tree.age()

map@phylo |> remove_sa() |> dateNodes()

df
mean(df$max_age - df$min_age)

rf_dist = phangorn::RF.dist(map@phylo, truth, rooted = TRUE, normalize = TRUE)
