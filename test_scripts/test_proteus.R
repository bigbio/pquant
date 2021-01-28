#! /usr/bin/R

library(proteus)
setwd(getwd())  # The path where you store this R script
evi <- read.csv("./evidence.csv")
meta <- read.csv("./metadata.csv")
pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)
prodat.med <- normalizeData(prodat)

# PX from the meta's condition, and you can only compare two conditions once.
res <- limmaDE(prodat.med, conditions=c("P1","P2"))

# draw volcano plot
plotVolcano(res)