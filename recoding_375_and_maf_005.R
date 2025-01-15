#!/usr/bin/env Rscript
# Load required packages
library(tidyverse)
library(collapse)
library(arrow)
library(data.table)
library(optparse)

# Set up command line options - we only need the input parameter
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

# Read the input file instead of directly assigning the path
dat1 <- read_parquet(opt$input)

# Rest of your processing code remains exactly the same
# First, let's identify all the metabolite columns (those starting with "X")
metabolite_columns <- names(dat1)[startsWith(names(dat1), "X")]

# Step 1: Calculate the presence percentage for each identifier and each metabolite
presence_summary <- dat1 %>%
  group_by(identifier) %>%
  summarise(
    across(
      all_of(metabolite_columns),
      list(
        total = ~n(),
        non_zero = ~sum(. != 0),
        percentage = ~(sum(. != 0) / n()) * 100
      )
    )
  )

valid_identifiers_matrix <- sapply(metabolite_columns, function(col) {
  percentage_col <- paste0(col, "_percentage")
  presence_summary[[percentage_col]] >= 37.5
})

identifier_matches <- match(dat1$identifier, presence_summary$identifier)
keep_values <- valid_identifiers_matrix[identifier_matches, ]

dat1_processed <- dat1
dat1_processed[metabolite_columns] <- Map(
  function(col, keep) ifelse(keep, dat1[[col]], 0),
  metabolite_columns,
  data.frame(keep_values)
)

filter_metabolites <- function(data, threshold = 0.05) {
  metabolite_cols <- grep("^X", colnames(data), value = TRUE)
  total_unique_ids <- length(unique(data$identifier))
  
  calculate_nonzero_proportion <- function(col_name) {
    nonzero_identifiers <- data[data[[col_name]] > 0, ]$identifier |> unique() |> length()
    prop_nonzero <- nonzero_identifiers / total_unique_ids
    cat(sprintf("Column %s: %d/%d (%.2f%%) identifiers have non-zero values\n", 
                col_name, nonzero_identifiers, total_unique_ids, prop_nonzero * 100))
    return(prop_nonzero > threshold)
  }
  
  cols_to_keep <- vapply(metabolite_cols, calculate_nonzero_proportion, logical(1))
  metabolites_to_keep <- metabolite_cols[cols_to_keep]
  
  all_cols_to_keep <- union(
    setdiff(colnames(data), metabolite_cols),
    metabolites_to_keep
  )
  
  cat(sprintf("\nKept %d out of %d metabolite columns (those with >%g%% non-zero values)\n", 
              sum(cols_to_keep), length(metabolite_cols), threshold * 100))
  
  data[, all_cols_to_keep, drop = FALSE]
}

filtered_df <- filter_metabolites(dat1_processed)

# Write back to the same input file, replacing it
write_parquet(filtered_df, opt$input)