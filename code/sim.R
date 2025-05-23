# install.packages("TreeSim")
# install.packages("FossilSim")
# install.packages("ape")
# install.packages("geiger")
library(FossilSim)
library(ape)
library(geiger)
library(TreeSim)


set.seed(1)
SA = FossilSim::sim.fbd.age(age = 2, # 2 myr record
                            numbsim = 1, # one tree simulated
                            lambda = 1, # origination rate
                            mu = 0.5, # extinction rate
                            psi = 2, # sampling rate
                            frac = 1, # sampling fraction after `age` years
                            complete = FALSE, # return complete tree? if true, returns unsampled lineages 
                            mrca = TRUE)[[1]] # start from one linage of the mrcs

plot(SA)
axis(1)
# write example tree
ape::write.nexus(SA, file = "data/test_tree.nex")

tree = TreeSim::sim.bd.age(age = 2,
                           numbsim = 1,
                           lambda = 1,
                           mu = 0.5,
                           mrca = TRUE,
                           complete = FALSE)[[1]]
f = sim.fossils.poisson(rate = 2, tree = tree)

plot(f, tree = tree)

sa = SAtree.from.fossils(tree, f, rho = 1)

sampled.tree.from.combined(tree = sa$tree) |> plot()
# simulate character matrix
par = matrix(c(-1, 1, 1, -1), nrow = 2, ncol = 2)
char_mat = geiger::sim.char(SA, par = par, model = "discrete", nsim = 1000)
char_mat = apply(char_mat, 1, function(x) x -1)
ape::write.nexus.data(t(char_mat), file = "data/char_mat.nex", format = "standard")


# export ages  
x = pmax(2 -  ape::node.depth.edgelength(SA)[seq_len(length(SA$tip.label))], 0)
df = data.frame(taxon = SA$tip.label, min_age = x |> signif(digits = 3), max_age = x |> signif(digits = 3))
write.table(df, file = "data/taxa.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


