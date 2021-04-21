# testing consensus output between MSstats and Proteus

library(VennDiagram)

source("R/functions.R")

project_id <- 'UPS1'
output_dir <- "./output"

#### 1. Importing and pre-processing quant data
path_quant_data <- "datasets/out_msstats.csv"

# import quant data (optput from quantification pipeline)
data <- import_quant_data(path_quant_data)

# format condition column (keep original)
data <- format_condition(data, condition_field = "Condition")

# format long protein groups names (keep original)
data <- format_protein_ids(data, protein_field = "ProteinName")


#### 2. Processing quant data with Proteus
# format quant data to Proteus
path_proteus <- quantdata2proteus(data)

# generate Proteus evidence file
data.proteus <- readEvidenceFile(path_proteus,
                                 data.cols = updateEvidenceColumns(),
                                 measure.cols = updateMeasureColumns())

# generate metadata file from evidence file
meta <- getMetaData(data.proteus)

# generate peptide table
pepdat <- makePeptideTable(data.proteus, meta)

# generate protein table
prodat <- makeProteinTable(pepdat)

# normalize intensities
prodat.med <- normalizeData(prodat)

# run differencial expression analysis (only two conditions)
res.proteus <- limmaDE(prodat.med, conditions = c('C1', 'C2'))


#### 3. Processing quant data with MSstats

# process quantitification raw data (normalization and imputation)
data.msstats <- dataProcess(raw = data,
                            normalization = 'equalizeMedians',
                            summaryMethod = 'TMP',
                            censoredInt = "NA",
                            cutoffCensored = "minFeature",
                            MBimpute = TRUE,
                            maxQuantileforCensored=0.999)

# set condition combination
# TODO: function to generate contrast matrix for comparitions
comparison <- matrix(c(-1,1,0,0,0,0,0,0,0),nrow=1)
row.names(comparison) <- "C1-C2"

# run differencial expression analysis (only two conditions)
res.msstats <- groupComparison(contrast.matrix = comparison, data = data.msstats)


##### 4. Consensus analysis

de.proteus <- getDifferentialExpressionResults(object = res.proteus, method = "proteus")
de.msstats <- getDifferentialExpressionResults(object = res.msstats, method = "msstats")


# Venn up-regulated (consensus)
proteus.up.set <- subset(de.proteus, logFC > 0 & p_value_adj < 0.05)$protein
msstats.up.set <- subset(de.msstats, logFC > 0 & p_value_adj < 0.05)$protein

venn.diagram(
  x = list(as.character(proteus.up.set), as.character(msstats.up.set)),
  main = 'Up-regulated (logFC > 0, Adjusted P < 0.05)',
  category.names = c("Proteus" , "MSstats"),
  filename = 'venn_diagramm_up.png',
  fill=c('blue', 'grey'),
  output=TRUE
)

# Venn Up-regulated (consensus)
proteus.up.set <- subset(de.proteus, logFC > 0 & p_value_adj < 0.05)$protein
msstats.up.set <- subset(de.msstats, logFC > 0 & p_value_adj < 0.05)$protein

venn.diagram(
  x = list(as.character(proteus.up.set), as.character(msstats.up.set)),
  main = 'Up-regulated (logFC > 0, Adjusted P < 0.05)',
  category.names = c("Proteus" , "MSstats"),
  height = 800,
  width = 800,
  resolution = 150,
  filename = 'venn_diagramm_up.png',
  fill=c('blue', 'grey'),
  output=TRUE
)

# Venn Down-regulated (consensus)
proteus.down.set <- subset(de.proteus, logFC < 0 & p_value_adj < 0.05)$protein
msstats.down.set <- subset(de.msstats, logFC < 0 & p_value_adj < 0.05)$protein

venn.diagram(
  x = list(as.character(proteus.down.set), as.character(msstats.down.set)),
  main = 'Down-regulated (logFC < 0, Adjusted P < 0.05)',
  category.names = c("Proteus" , "MSstats"),
  height = 800,
  width = 800,
  resolution = 150,
  filename = 'venn_diagramm_down.png',
  fill=c('blue', 'grey'),
  output=TRUE
)



