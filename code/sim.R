install.packages("TreeSim")
install.packages("FossilSim")
install.packages("ape")
install.packages("geiger")
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)
browseVignettes("FossilSim")

set.seed(1)
t = TreeSim::sim.bd.age(age = 2, numbsim = 1, lambda = 2, mu = 01, frac = 1)[[1]]
plot.phylo(t)
f = sim.fossils.poisson(rate =  2, tree = t)
plot(f, tree = t)
# create sampled ancestor tree
SAt = SAtree.from.fossils(tree = t, fossils = f)

# remove unsampled lineages from the tree
SA = sampled.tree.from.combined(SAt$tree)


plot(SA)
par = matrix(c(-1, 1, 1, -1), nrow = 2, ncol = 2)
char_mat = geiger::sim.char(SA, par = par, model = "discrete", nsim = 10)

char_mat = apply(char_mat, 1, function(x) x -1)
ape::write.nexus.data(char_mat, file = "char_mat.nex", format = "standard")
