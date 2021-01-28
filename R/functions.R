# Useful function for parsing quantification data

library(dplyr)
library(proteus)
library(MSstats)

# import CSV dataset
import_quant_data <- function(path) {
  data <- read.csv(path, stringsAsFactors = F)
  col_names <- names(data)
  expected_cols <-
    c(
      "ProteinName",
      "PeptideSequence",
      "PrecursorCharge",
      "FragmentIon",
      "ProductCharge",
      "IsotopeLabelType",
      "Condition",
      "BioReplicate",
      "Run",
      "Intensity",
      "Reference"
    )
  if (!all(col_names %in% expected_cols)) {
    stop(paste0(
      "Misformated input quantification data. Expected columns: ",
      paste(expected_cols, collapse = ",")
    ))
  }
  return(data)
}

# rename conditions name
format_condition <- function(df, condition_field = 'Condition') {
  if(!condition_field %in% names(df)) stop("Condition column not found...")
  conds <- as.character(df[, condition_field])
  unique_conds <- unique(conds)
  conds_formated <- factor(conds,
                           levels = unique_conds,
                           labels = paste0("C", 1:length(unique_conds)))
  df[, condition_field] <- as.character(conds_formated)
  df[, paste0(condition_field, '_original')] <-
    conds # keep original name conditions
  return(df)
}

# count char occurences
count_char <- function(string, char = ";") {
  char_vector <- unlist(strsplit(string, "", fixed = T))
  return(length(char_vector[char_vector == char]))
}

# format protein groups IDs
format_protein_ids <- function(df, protein_field = 'ProteinName', separator = ';') {
  if(!protein_field %in% names(df)) stop("ProteinName column not found...")
    proteins <- as.character(df[, protein_field])
    unique_proteins <- unique(proteins)

    # format (only) protein groups
    protein_labels <- vector()
    for(index in 1:length(unique_proteins)){
        p <- unique_proteins[index]
        if (count_char(p, char = separator) >= 1) {
           p <- unlist(strsplit(p, separator, fixed = T))[1]
           p <- paste0(p, "-PG:", index)
           protein_labels[[index]] <- p
        } else {
           protein_labels[[index]] <- p
        }
     }

    proteins_formated <- factor(proteins, levels = unique_proteins, labels = protein_labels)
    df[, protein_field] <- as.character(proteins_formated)
    df[, paste0(protein_field, '_original')] <- proteins
  return(df)
}

# rename evidence columns
updateEvidenceColumns <- function(){

  evidence_cols <- proteus::evidenceColumns

  # defaults
  evidence_cols$sequence <- "PeptideSequence"
  evidence_cols$modified_sequence <- "ModifiedSequence"
  evidence_cols$modifications <- "Modifications"
  evidence_cols$protein_group <- "ProteinGroup"
  evidence_cols$protein <- "ProteinName"
  evidence_cols$experiment <- "Run"
  evidence_cols$condition <- "Condition"
  evidence_cols$charge <- "PrecursorCharge"
  evidence_cols$reverse <- "Reverse"
  evidence_cols$contaminant <- "Contaminant"

  # other information
  evidence_cols$fragment_ion <- "FragmentIon"
  evidence_cols$isotope_label_type <- "IsotopeLabelType"
  evidence_cols$replicate <- "BioReplicate"
  evidence_cols$raw_mz <- "Reference"

  return(evidence_cols)
}

# rename measure column
updateMeasureColumns <- function(){
  measure_cols <- proteus::measureColumns
  measure_cols$intensity <- "Intensity"
  return(measure_cols)
}

# format quant data to compatibility with Proteus
quantdata2proteus <- function(path){
  df <- read.csv(path, stringsAsFactors = F)
  # TODO:
  # Set to missing the information not exported by the current pipeline.
  # This is a temporary solution, better to export this information and
  # make the output from the current pipeline compatible with Proteus
  df['ModifiedSequence'] <- NA
  df['Modifications'] <- NA
  df['ProteinGroup'] <- NA
  df['Reverse'] <- NA
  df['Contaminant'] <- NA
  output_path <- paste0(path, 'proteus_formatted.tsv.gz')
  write.table(df, file = output_path, sep = '\t', row.names = F, quote = F)
  message(paste0("Exported Proteus-formatted evidence file to: ", output_path))
  return(output_path)
}

# generate metadata from evidence file (Proteus)
getMetaData <- function(df) {
  key_cols <- c('experiment', 'condition', 'replicate')
  if(!all(key_cols %in% names(df))) stop(paste0("Column not found. Expected: ", paste0(key_cols, collapse = ",")))
  df <- subset.data.frame(df, select = c('experiment', 'condition', 'replicate'))
  meta <- dplyr::group_by(df, experiment, condition, replicate) %>% summarise()
  meta <- dplyr::mutate(meta, sample=paste0(condition,'-',experiment),
                              measure='Intensity')
  return(meta)
}






