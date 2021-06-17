library('MSstats', warn.conflicts = F, quietly = T, verbose = F)

setwd('D:/dataset/R downstream analysis/pquant/data')
load("preShiny.RData")

setwd('D:/dataset/R downstream analysis/pquant/')
source("shiny-app/app.R")

pquant_shiny(DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons)