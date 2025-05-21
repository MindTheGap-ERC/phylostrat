# install.packages("TreeSim")
# install.packages("FossilSim")
# install.packages("ape")
# install.packages("geiger")
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)


set.seed(1)
SA = FossilSim::sim.fbd.age(age = 2, numbsim = 1, lambda = 1, mu = 0, psi = 2, frac = 1, complete = TRUE)[[1]]

ape::write.nexus(SA, file = "data/test_tree.nex")
plot(SA)
par = matrix(c(-1, 1, 1, -1), nrow = 2, ncol = 2)
char_mat = geiger::sim.char(SA, par = par, model = "discrete", nsim = 1000)

char_mat = apply(char_mat, 1, function(x) x -1)
ape::write.nexus.data(t(char_mat), file = "data/char_mat.nex", format = "standard")
