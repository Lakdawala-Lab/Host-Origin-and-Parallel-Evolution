# Title: Processing and selecting sequences
# AUTHORS

# Jennifer E. Jones <jenny@jennyjonesphd.com>
# Seema S. Lakdawala <lakdawala@pitt.edu>
  
# This analysis requires the DECIPHER package. This package and its dependencies
# can be downloaded from Bioconductor at http://bioconductor.org/packages/release/bioc/html/DECIPHER.html.

# Load libraries.

library(DECIPHER) # Used for subsetting, alignment, and clustering of sequences.

sessionInfo() # Need DECIPHER >= v2.20.0

# Set working directory.

setwd("<<path to FASTA file directory>>/") # Needs trailing slash.
# Example "~/Desktop/Host-Origin-and-Parallel-Evolution-main/"

# FASTA FILE IMPORT AND CLEANUP

# Load FASTA sequences for analysis. Datasets analyzed in this manuscript can be found at 
# https://github.com/Lakdawala-Lab/Host-Origin-and-Parallel-Evolution/Data/Pre-processed FASTA Files.
# Avian H9 viruses are used as an example here.

subtype <- "H9"
segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
nSegments <- length(segments)

files <- paste(paste0("./", subtype), segments, "segment.fasta.gz")

# FASTA file QC.

length <- length(files)
mRNA <- vector(mode = "list", length = length)
vRNA <- vector(mode = "list", length = length)
strains <- vector(mode = "list", length = length)

for (i in seq_along(vRNA)) {
  mRNA[[i]] <- readDNAStringSet(files[i]) # IAV fasta files are downloaded as positive-sense (mRNA) sequences. 
  vRNA[[i]] <- reverseComplement(mRNA[[i]]) # Convert sequence files to the genomic (vRNA) sequences.

  # Extract the strain names from FASTA files. This step is required to ensure accurate concatenation of full-length genomes.
  strains[[i]] <- gsub(".+\\|Strain Name:(.+?)\\|Segment.+", "\\1", names(vRNA[[i]])) 
  
  # Remove duplicate entries from FASTA files.
  d <- which(!duplicated(strains[[i]]))
  vRNA[[i]] <- vRNA[[i]][d]
}

# Selecting strains from specific hosts.

hosts <- gsub(".+\\|Host:", "\\1", names(vRNA[[1]]))
unique(hosts) # Human, Environment, Equine, Ferret, mink, Mink, Weasel, Unknown will all be removed.

# Avian hosts can be subset further into aquatic bird hosts or landfowl hosts as described in the manuscript.

hosts <- vector(mode = "list", length = length)
for (i in seq_along(vRNA)) {
  hosts[[i]] <- gsub(".+\\|Host:", "\\1", names(vRNA[[i]]))
  w <- c(which(hosts[[i]] == "Human"), which(hosts[[i]] == "Environment"),
         which(hosts[[i]] == "Equine"), which(hosts[[i]] == "Ferret"), which(hosts[[i]] == "mink"),
         which(hosts[[i]] == "Mink"), which(hosts[[i]] == "Weasel"), which(hosts[[i]] == "Unknown"))
  vRNA[[i]] <- vRNA[[i]][-w]
  names(vRNA[[i]]) <- gsub(".+\\|Strain Name:(.+?)\\|Segment.+", "\\1", names(vRNA[[i]])) 
}

# Remove FASTA files with incomplete genomic sequences.

t <- table(unlist(sapply(vRNA, function(x) unique(names(x)))))
t <- names(t)[t == nSegments]
vRNA <- lapply(vRNA, `[`, t)

# Assemble full-length genomes by concatenation.

vRNA[[9]] <- do.call(xscat, vRNA)
names(vRNA[[9]]) <- names(vRNA[[1]])

# Align sequences.

vRNA <- lapply(vRNA, AlignSeqs, processors = NULL) 

# Drop sequences with ambiguities (N, R, Y, W, S, M, or K).

ambig <- lapply(as.character(vRNA[[9]]), strsplit, "N")
ambigN <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigN[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigN <- which(!is.na(ambigN))

ambig <- lapply(as.character(vRNA[[9]]), strsplit, "R") 
ambigR <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigR[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigR <- which(!is.na(ambigR))

ambig <- lapply(as.character(vRNA[[1]]), strsplit, "Y") 
ambigY <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigY[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigY <- which(!is.na(ambigY))

ambig <- lapply(as.character(vRNA[[9]]), strsplit, "W") 
ambigW <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigW[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigW <- which(!is.na(ambigW))

ambig <- lapply(as.character(vRNA[[9]]), strsplit, "S") 
ambigS <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigS[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigS <- which(!is.na(ambigS))

ambig <- lapply(as.character(vRNA[[9]]), strsplit, "M") 
ambigM <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigM[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigM <- which(!is.na(ambigM))

ambig <- lapply(as.character(vRNA[[9]]), strsplit, "K") 
ambigK <- numeric()
for (i in c(1 : length(ambig))) {
  # runs uptil the length of inner lists at ith indices
  for (j in c(1: length(ambig[[i]])))
  { 
    for (k in c(1: length(ambig[[i]][[j]])))
      if (k > 1) {
        ambigK[i] <- i
        #  cat("Strain", i, "Ambiguities", k, "; ")
      }
  }
}
ambigK <- which(!is.na(ambigK))
w <- unique(c(ambigN, ambigY, ambigR, ambigW, ambigS, ambigM, ambigK))
vRNA <- lapply(vRNA, `[`, -w)

# CLUSTERING INTO OPERATIONAL TAXONOMIC UNITS (OTUs) AND SEQUENCE SELECTION

# Calculate distances between full-length concatenated genomic sequences 
# and compare clusters with different cutoffs ranging from 86-99% sequence identity.

d <- DistanceMatrix(vRNA[[9]], type="dist", correction="JC", processors = NULL)
otu <- IdClusters(d, method = "NJ", cutoff = c(0.01, 0.02, 0.03, 0.04, 0.05, 
                                               0.06, 0.07, 0.08, 0.09, 0.1, 
                                               0.11, 0.12, 0.13, 0.14), 
                  type = "clusters", myXStringSet = vRNA[[9]], processors = NULL)

# View the number of clusters in each species tree with each cutoff.

sapply(otu, max)

# Choose the desired cutoff (in this case, 95% sequence identity was selected).

otu <- otu[, 5, drop=FALSE]

# Write clustering data to working directory.

write.table(otu1, file = "Avian H9 Clusters.csv", sep=",")

# Manually choose representative sequences from each cluster. Then, subset the alignments 
# of the sequences selected. The sequences analyzed here are indicated in the 
# 'Avian H9 Strains Analyzed.csv' file, which can be found in the following sub-folder:
# https://github.com/Lakdawala-Lab/Host-Origin-and-Parallel-Evolution/Data/Pre-processed FASTA Files.

id <- read.csv("Avian H9 Strains Analyzed.csv")
id <- id[2:length(id)]

for (i in seq_along(files)) {
  for (j in seq_along(id)) {
    m <- match(names(vRNA[[i]]), id[, j])
    w <- which(!is.na(m))
    writeXStringSet(vRNA[[i]][w], file= paste0(substring(files[i], 1, nchar(files[i]) - nchar("sequences.fasta.gz")), "avian vRNA MSA OTU 95 ", j, ".fasta.gz", sep=""))
  }
}

# Fully processed FASTA files are now ready for tree reconstruction and analysis. Source code to complete this analysis
# is provided in the 'Tree Reconstruction and Analysis' sub-folder.





