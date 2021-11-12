## test Proteus

source("R/functions.R")

project_id <- 'UPS1'
output_dir <- paste0("output/", "proteus-", project_id)

path_quant_data <- "datasets/out_msstats.csv"

# format quant data to Proteus
path_proteus <- quantdata2proteus(path_quant_data)

# generate Proteus evidence file
evi <- readEvidenceFile(path_proteus,
                        data.cols = updateEvidenceColumns(),
                        measure.cols = updateMeasureColumns())

# format condition column (keep original)
evi <- format_condition(evi, condition_field = "condition")

# format long protein groups names (keep original)
evi <- format_protein_ids(evi, protein_field = "protein")

# generate metadata file from evidence file
meta <- getMetaData(evi)

# generate peptide table
pepdat <- makePeptideTable(evi, meta)

# generate protein table
prodat <- makeProteinTable(pepdat)

# normalize intensities
prodat.med <- normalizeData(prodat)

# run differencial expression analysis (only two conditions)
res <- limmaDE(prodat.med, conditions = c('C1', 'C2'))

# visualization

# shini volcano plot
# plotVolcano_live(prodat.med, res)

# vocano plot
pdf(paste0(output_dir, "_volcano.pdf"))
proteus::plotVolcano(res)
dev.off()

# FID plot
pdf(paste0(output_dir, "_fid.pdf"))
proteus::plotFID(prodat.med, pair = c('C1', 'C2'))
dev.off()

# P-value distribution plot
pdf(paste0(output_dir, "_pval_distribution.pdf"))
proteus::plotPdist(res)
dev.off()


