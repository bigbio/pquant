## test MSstats

source("R/functions.R")

project_id <- 'UPS1'
output_dir <- paste0("output/", project_id)

path_quant_data <- "datasets/out_msstats.csv"

# import quant data (optput from quantification pipeline)
data <- import_quant_data(path_quant_data)

# format condition column (keep original)
data <- format_condition(data, condition_field = "Condition")

# format long protein groups names (keep original)
data <- format_protein_ids(data, protein_field = "ProteinName")

# process quantitification raw data (normalization and imputation)
data.processed <- dataProcess(raw = data,
                              normalization = 'equalizeMedians',
                              summaryMethod = 'TMP',
                              censoredInt = "NA",
                              cutoffCensored = "minFeature",
                              MBimpute = TRUE,
                              maxQuantileforCensored=0.999)

# use type="QCplot" with all proteins
# change the upper limit of y-axis=35
# set up the size of pdf
dataProcessPlots(data=data.processed,
                 type="QCplot",
                 ylimUp=35,
                 width=8,
                 height=6,
                 text.size=4,
                 address=output_dir)


# use type="Profileplot"
dataProcessPlots(data = data.processed,
                 type="Profileplot",
                 ylimUp=35,
                 featureName="NA",
                 width=8,
                 height=6,
                 text.size=4,
                 address=output_dir)

# use type="Conditionplot"
dataProcessPlots(data = data.processed,
                 type="Conditionplot",
                 width=8,
                 height=6,
                 text.size=4,
                 address=output_dir)

# testing two possible condition combination
# TODO: function to generate contrast matrix for comparitions
comparison1 <- matrix(c(-1,1,0,0,0,0,0,0,0),nrow=1)
comparison2 <- matrix(c(0,-1,1,0,0,0,0,0,0),nrow=1)
comparison <- rbind(comparison1, comparison2)

row.names(comparison) <- c("C1-C2", "C2-C3")
data.comparisons <- groupComparison(contrast.matrix = comparison, data = data.processed)


# normal quantile-quantile plots
modelBasedQCPlots(data=data.comparisons, type="QQPlots",
                  text.size=4, address=output_dir)

# residual plots
modelBasedQCPlots(data=data.comparisons, type="ResidualPlots",
                 text.size=4, address=output_dir)

# volcanoPlots
groupComparisonPlots(data = data.comparisons$ComparisonResult, type = 'VolcanoPlot',
                     ProteinName=F, text.size=4, address=output_dir)

# Heatmaps
groupComparisonPlots(data = data.comparisons$ComparisonResult, type = 'Heatmap',
                     text.size = 4, address = output_dir)

# ComparisonPlot
groupComparisonPlots(data=data.comparisons$ComparisonResult, type="ComparisonPlot",
                     width=8, height=6, address=output_dir)
