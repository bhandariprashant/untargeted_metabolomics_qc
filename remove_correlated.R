#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(collapse)
library(arrow)
library(data.table)
library(optparse)

# Set up command line options
option_list <- list(
  make_option(c("-i", "--input"), 
              type="character", 
              help="Input RDS file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Please provide an input file path using -i or --input")
}

# Main workflow with basic error handling
tryCatch({
  # Read input data
  dat1 <- readRDS(opt$input)
  
  # Extract metabolite columns
  metabolite_cols <- grep("^X", names(dat1), value = TRUE)
  metabolite_data <- dat1[, metabolite_cols]
  
  # Calculate correlation matrix
  cor_matrix <- cor(metabolite_data)
  
  # Convert to pairs format
  cor_pairs <- which(upper.tri(cor_matrix), arr.ind = TRUE)
  cor_results <- data.frame(
    metabolite1 = rownames(cor_matrix)[cor_pairs[,1]],
    metabolite2 = colnames(cor_matrix)[cor_pairs[,2]],
    correlation = cor_matrix[cor_pairs]
  )
  
  high_cor <- cor_results[cor_results$correlation > 0.9, ]
  
  # Extract masses and identify metabolites to remove
  high_cor$mass1 <- as.numeric(str_extract(high_cor$metabolite1, "(?<=X)\\d+\\.\\d+"))
  high_cor$mass2 <- as.numeric(str_extract(high_cor$metabolite2, "(?<=X)\\d+\\.\\d+"))
  
  saveRDS(high_cor,
          paste0(basename(opt$input),"_high_cor.rds"))
  
  metabolites_to_remove <- unique(
    ifelse(high_cor$mass1 >= high_cor$mass2,
           as.character(high_cor$metabolite1),
           as.character(high_cor$metabolite2))
  )
  
  # Filter dataset and save results
  dat1_filtered <- dat1[, !names(dat1) %in% metabolites_to_remove]
  
  saveRDS(metabolites_to_remove,
          paste0(basename(opt$input),"_metabolites_removed.rds"))
  
  write_parquet(dat1_filtered, 
                paste0(basename(opt$input),"_metabolites_filtered.parquet"))
  
}, error = function(e) {
  stop(paste("Error in processing:", e$message))
})