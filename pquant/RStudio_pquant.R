### get a "/pquant/data/" path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./data/")
fileData <- read.csv('out_msstats.csv')

source("../RStudio_prePquant.R")
prePquant(fileData)

# *** if you already have a "preShiny.RData", execute codes from here ***
load("preShiny.RData")

source("../shiny-app/app.R")
pquant_shiny(DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons)
