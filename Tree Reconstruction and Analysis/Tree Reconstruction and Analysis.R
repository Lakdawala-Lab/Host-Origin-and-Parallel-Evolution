# Title: Tree Reconstruction and Analysis
# AUTHORS

# Jennifer E. Jones <jenny@jennyjonesphd.com>
# Seema S. Lakdawala <lakdawala@pitt.edu>

# Load libraries.

library(ape)
library(phangorn)
library(TreeDist)

sessionInfo() # Need DECIPHER >= v2.20.0

# Set working directory.

setwd("<<PATH TO Parallel-Evolution>>/") # Needs trailing slash.
# Example "~/Desktop/Host-Origin-and-Parallel-Evolution-main/"

# RECONSTRUCTING PHYLOGENETIC TREES

# Analysis of avian H9 virus gene sequences is demonstrated here.

# First perform model testing on alignments for maximum-likelihood method. Alignments can be loaded from
# https://github.com/Lakdawala-Lab/Host-Origin-and-Parallel-Evolution/Data/Post-processing Alignments

files <- c("./H9 PB2 avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 PB1 avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 PA avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 HA avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 NP avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 NA avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 M avian vRNA OTU 95 MSA.fasta.gz",
           "./H9 NS avian vRNA OTU 95 MSA.fasta.gz")

# Pre-populate empty lists for use later in the analysis.

length <- length(files)
fit <- vector(mode = "list", length = length) # A list of the fitted trees.
treeBS <- vector(mode = "list", length = length) # A list of the bootstrapped trees.
mtsort <- data.frame("LogLik" = vector(mode = "character", length = length), # Sorting the results from model testing.
                     "AIC" = vector(mode = "character", length = length),
                     "BIC" = vector(mode = "character", length = length))
mt2<-vector(mode = "list", length = length) # Sorting the results from model testing.
mt2<-lapply(mt2, function(x) x = vector(mode = "list", length = 3))

# Files are first read into ape and then passed to phangorn as a phyDat object.

vRNA <- lapply(files, read.dna, format = "fasta")
vRNApd <- lapply(vRNA, phyDat, type = "DNA", levels = NULL)

# Calculate the distance matrix from each alignment and build an initial neighbor-joining tree.

d <- lapply(vRNApd, dist.ml, model = "JC69")
vRNA_NJ <- lapply(d, NJ)

# Compute the likelihood of the starting tree using the pml function.

for (i in seq_along(fit)) {
  fit[[i]] <- pml(vRNA_NJ[[i]], vRNApd[[i]])
  print(fit[[i]])
}

# Perform model testing on each multiple sequence alignment.

mt <- lapply(vRNApd, modelTest)

# Choose the most appropriate model for each tree based on the lowest AIC.

for (i in seq_along(mt2)) {
  mt2[[i]][[1]] <- mt[[i]][order(mt[[i]]$logLik, decreasing = T),]
  mtsort$LogLik[i] <- (mt2[[i]][[1]]$Model[1])
  
  mt2[[i]][[2]] <- mt[[i]][order(mt[[i]]$AIC, decreasing = F),]
  mtsort$AIC[i] <- (mt2[[i]][[2]]$Model[1])
  
  mt2[[i]][[3]] <- mt[[i]][order(mt[[i]]$BIC, decreasing = F),]
  mtsort$BIC[i] <- (mt2[[i]][[3]]$Model[1])
}

# Optimize each tree for the best model. In this case, GTR was used for all trees.

fit2 <- lapply(fit, optim.pml, model = "GTR", optNni=TRUE, rearrangement = "stochastic")

# Where additional parameters are recommended during model testing (e.g., '+G', '+I', or '+G+I'), 
# extract the specified parameters and update the tree. Individual parameters
# will vary for each tree and should be updated according to the model testing results.

w <- which(mtsort$AIC == "GTR+G")
fit2G <- fit2[w]
mtsortG <- mtsort[w,]
mtG <- mt[w]

w <- which(mtsort$AIC == "GTR+I")
fit2I <- fit2[w]
mtsortI <- mtsort[w,]
mtI <- mt[w]

w <- which(mtsort$AIC == "GTR+G+I")
fit2GI <- fit2[w]
mtsortGI <- mtsort[w,]
mtGI <- mt[w]

for (i in seq_along(fit2G)) {
  env <- attr(mtG[[i]], "env")
  fit2G[[i]] <- update(fit2G[[i]],
                       k = (eval(get(mtsortG[i,][[2]], env), env))$k, 
                       shape = (eval(get(mtsortG[i,][[2]], env), env))$shape)
}

for (i in seq_along(fit2I)) {
  env <- attr(mtI[[i]], "env")
  fit2I[[i]] <- update(fit2I[[i]],
                       inv = (eval(get(mtsortI[i,][[2]], env), env))$inv)
}

for (i in seq_along(fit2GI)) {
  env <- attr(mtGI[[i]], "env")
  fit2GI[[i]] <- update(fit2GI[[i]],
                        k = (eval(get(mtsortGI[i,][[2]], env), env))$k, 
                        shape = (eval(get(mtsortGI[i,][[2]], env), env))$shape,
                        inv = (eval(get(mtsortGI[i,][[2]], env), env))$inv)
}

# Optimize the trees with the specified parameters.

fit3G <- lapply(fit2G, optim.pml, model = "GTR", optGamma=TRUE, optNni=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
fit3I <- lapply(fit2I, optim.pml, model = "GTR", optNni=TRUE, optInv=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
fit3GI <- lapply(fit2GI, optim.pml, model = "GTR", optGamma=TRUE, optNni=TRUE, optInv=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

# Apply bootstraps (1,000 were used here) and test for support.

bsG <- lapply(fit3G, bs=1000, optNni=TRUE, optGamma=TRUE, multicore=TRUE, control = pml.control(trace=0))
bsI <- lapply(fit3I, bootstrap.pml, bs=1000, optNni=TRUE, optInv=TRUE, multicore=TRUE, control = pml.control(trace=0))
bsGI <- lapply(fit3GI, bootstrap.pml, bs=1000, optNni=TRUE, optGamma=TRUE, optInv=TRUE, multicore=TRUE, control = pml.control(trace=0))

# Merge the fitted and bootstrapped tree lists for plotting and export.

treeBSg <- vector(mode = "list", length = length(fit3G)) # A list of the bootstrapped trees.
for (i in seq_along(treeBSg)) {
  treeBSg[[i]] <- plotBS(midpoint(fit3G[[i]]$tree), bsG[[i]], p = 50, type = "p")
  write.tree(treeBSg[[i]], file = paste(substring(files[i], 1,  nchar(files[i]) - nchar(".fasta.gz")), " tree", sep=""))
}

treeBSi <- vector(mode = "list", length = length(fit2I)) # A list of the bootstrapped trees.
for (i in seq_along(treeBSi)) {
  treeBSi[[i]] <- plotBS(midpoint(fit3I[[i]]$tree), bsI[[i]], p = 50, type = "p")
  write.treei(treeBSgi[[i]], file = paste(substring(files[i], 1,  nchar(files[i]) - nchar(".fasta.gz")), " tree", sep=""))
}

treeBSgi <- vector(mode = "list", length = length(fit2GI)) # A list of the bootstrapped trees.
for (i in seq_along(treeBSgi)) {
  treeBSgi[[i]] <- plotBS(midpoint(fit3GI[[i]]$tree), bsGI[[i]], p = 50, type = "p")
  write.tree(treeBSgi[[i]], file = paste(substring(files[i], 1,  nchar(files[i]) - nchar(".fasta.gz")), " tree", sep=""))
}

# QUANTIFICATION OF TREE DISTANCES

files <- c("./H9 PB2 avian vRNA OTU 95 MSA tree", 
           "./H9 PB1 avian vRNA OTU 95 MSA tree",
           "./H9 PA avian vRNA OTU 95 MSA tree",
           "./H9 HA avian vRNA OTU 95 MSA tree",
           "./H9 NP avian vRNA OTU 95 MSA tree",
           "./H9 NA avian vRNA OTU 95 MSA tree",
           "./H9 M avian vRNA OTU 95 MSA tree",
           "./H9 NS avian vRNA OTU 95 MSA tree")

length <- length(files)
trees <- vector(mode = "list", length = length)

trees <- lapply(files, read.tree)
trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
cid <- TreeDistance(trees)
rownames(cid) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(cid) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
write.table(cid, file = "./Avian H9 OTU 95 CID.csv", sep=",", quote = FALSE, row.names = T)
