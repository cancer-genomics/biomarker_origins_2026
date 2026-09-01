#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript split_cached_tf_matrices.R index access_dir out_parent_dir")
}

index        <- as.integer(args[1])
access_dir   <- args[2]
out_parent   <- args[3]

cat("=== TF Cache Splitter ===\n")
cat("Index: ", index, "\n")
cat("Access Dir: ", access_dir, "\n")
cat("Output Parent Dir: ", out_parent, "\n\n")

# -------------------------------------------------------------------------
# 1. Find all cached matrix files
# -------------------------------------------------------------------------
cache_files <- sort(list.files(
  access_dir,
  pattern = "gc_corrected_tf_footprint_matrix_cache\\.rds$",
  full.names = TRUE
))

if (length(cache_files) == 0) {
  stop("No cached matrix files found in access_dir: ", access_dir)
}

if (index < 1 || index > length(cache_files)) {
  stop("Index out of range: ", index, " (n=", length(cache_files), ")")
}

cache_file <- cache_files[index]
cat("Processing file: ", basename(cache_file), "\n")

# -------------------------------------------------------------------------
# 2. Extract sample ID
# -------------------------------------------------------------------------
sample_id <- sub("\\..*$", "", basename(cache_file))
sample_dir <- file.path(out_parent, sample_id)

if (!dir.exists(sample_dir)) {
  dir.create(sample_dir, recursive = TRUE)
}

# If directory already contains metadata, assume completed
metadata_path <- file.path(sample_dir, "metadata.rds")
if (file.exists(metadata_path)) {
  cat("✓ Sample already processed, skipping:", sample_id, "\n")
  quit(save="no")
}

# -------------------------------------------------------------------------
# 3. Load cached matrix file
# -------------------------------------------------------------------------
cat("Loading full cached matrix file... (this may take time)\n")
cache <- readRDS(cache_file)

# -------------------------------------------------------------------------
# 4. Save metadata
# -------------------------------------------------------------------------
cat("Saving metadata...\n")
saveRDS(cache$metadata, metadata_path)

# -------------------------------------------------------------------------
# 5. Save each TF matrix independently
# -------------------------------------------------------------------------
tf_list  <- cache$matrices
tf_names <- names(tf_list)

cat("Saving", length(tf_names), "TF matrices...\n")

for (tf in tf_names) {
  outfile <- file.path(sample_dir, paste0("tf_", tf, ".rds"))
  saveRDS(tf_list[[tf]], outfile)
}

cat("✓ Completed sample:", sample_id, "\n")