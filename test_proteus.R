library(proteus)
evi <- read.csv("D:/dataset/R downstream analysis/proteus/ourdata/evidence.csv")
meta <- read.csv("D:/dataset/R downstream analysis/proteus/ourdata/metadata.csv")
setwd("D:/dataset/R downstream analysis/proteus/ourdata/")
pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)
prodat.med <- normalizeData(prodat)

# PX from the meta's condition, and you can only compare two conditions once.
res <- limmaDE(prodat.med, conditions=c("P1","P2"))

# draw volcano plot
plotVolcano(res)