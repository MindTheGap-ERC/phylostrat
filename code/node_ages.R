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

# read in data
id = 11
nchar = 100
library(paleotree)
library(dispRity)
    
    map = ape::read.tree(paste0(path_res, "tree_", id, "_nchar", nchar, "_MAP.tre"))
    truth = ape::read.nexus(paste0("data/sim/fbd_base/tree_", id, "_nchar", nchar , ".nex"))
    
    dateNodes(map)
    
a = tree.age(map)
sort(a$ages)

plot(truth, show.node.label = TRUE) 

plot(drop.tip(truth, "t39_1"))
?drop.tip
truth$edge.length



remove_sa = function(x){
  is_sa = x$edge.length == 0 #& truth$edge[,2] > Ntip(truth)
  sa_nodes = x$edge[is_sa, 2]
  phy_no_sa =  drop.tip(x, sa_nodes)
}

plot(map |> remove_sa())
plot(truth |> remove_sa())

map |> remove_sa() |> tree.age()
truth |> remove_sa() |> tree.age()

map |> remove_sa() |> dateNodes()


aa = readTrees("output/fbd_base/tree_10_nchar30_run_1.log", tree_name = "fbd_tree", n_cores = 4)
