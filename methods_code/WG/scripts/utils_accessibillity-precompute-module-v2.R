#!/usr/bin/env Rscript
# ===================================================================
# TF Footprint Pre-computation Pipeline (v2 – optimized + matrix output)
#
# Key features:
#   ✓ Supports .starch (via BEDOPS helper)
#   ✓ Preprocessed fragments REQUIRED (must contain weight column)
#   ✓ Saves EITHER:
#        (1) large combined footprint matrix (old mode), OR
#        (2) per-TF matrices in a per-sample folder (new mode)
#
# Usage:
#   Rscript tf_footprint_precompute_v2.R index fragDir tfbsDir outDir \
#           window bin_size frag_min frag_max platform save_mode
#
# save_mode ∈ {"combined", "per_tf"}
# ===================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(rtracklayer)
  library(data.table)
})

HG19_GENOME <- BSgenome.Hsapiens.UCSC.hg19
VALID_CHROMOSOMES <- paste0("chr", c(1:22, "X", "Y"))

# ======================== BEDOPS Helper ========================
BEDOPS_HELPER <- "/dcs04/scharpf/data/nvulpesc/tools/bedops_R.R"
if (!file.exists(BEDOPS_HELPER)) {
  stop("BEDOPS helper script not found: ", BEDOPS_HELPER)
}
source(BEDOPS_HELPER)
if (!exists("unstarch_to_granges")) {
  stop("Function unstarch_to_granges() missing after sourcing BEDOPS helper.")
}

# ======================== ARGUMENTS ============================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: Rscript tf_footprint_precompute_v2.R index fragDir tfbsDir outDir window bin_size frag_min frag_max platform save_mode")
}

index     <- as.numeric(args[1])
fragDir   <- args[2]
tfbsDir   <- args[3]
outDir    <- args[4]
window    <- as.numeric(args[5])
bin_size  <- as.numeric(args[6])
frag_min  <- as.numeric(args[7])
frag_max  <- as.numeric(args[8])
platform  <- tolower(args[9])
save_mode <- tolower(args[10])   # NEW ARGUMENT

if (!platform %in% c("hiseq","novaseq"))
  stop("Platform must be 'hiseq' or 'novaseq'")

if (!save_mode %in% c("combined", "per_tf"))
  stop("save_mode must be 'combined' or 'per_tf'")

# ======================== FRAGMENT LOADER ======================
load_fragments_preprocessed <- function(frag_path, frag_min, frag_max) {

  cat("Loading fragments:", basename(frag_path), "\n")

  frag <- if (grepl("\\.starch$", frag_path, ignore.case = TRUE)) {
    unstarch_to_granges(frag_path, tmpdir = tempdir())
  } else {
    readRDS(frag_path)
  }

  if (!length(frag)) stop("Fragment file empty: ", frag_path)

  # Normalize sequence names
  suppressWarnings({
    seqlevelsStyle(frag) <- "UCSC"
    common <- intersect(seqlevels(frag), VALID_CHROMOSOMES)
    if (length(common) > 0) {
      frag <- keepSeqlevels(frag, common, pruning.mode = "coarse")
      seqinfo(frag) <- keepSeqlevels(seqinfo(HG19_GENOME),
                                     seqlevels(frag),
                                     pruning.mode = "coarse")
    }
  })

  frag <- frag[seqnames(frag) %in% VALID_CHROMOSOMES]
  frag <- frag[width(frag) >= frag_min & width(frag) <= frag_max]

  if (!length(frag))
    stop("No fragments remain after filtering.")

  if (!("weight" %in% names(mcols(frag)))) {
    stop("Preprocessed fragments must contain a 'weight' column.")
  }

  return(frag)
}

# ======================== TF LOADER ============================
load_tf <- function(tf_name, tfbs_dir) {
  tf_clean <- toupper(tf_name)
  fp <- file.path(tfbsDir, paste0("tf_", tf_clean, "_5k.bed"))
  if (!file.exists(fp)) stop("TF file missing: ", fp)
  gr <- import(fp)
  seqlevelsStyle(gr) <- "UCSC"
  gr
}

# ======================== WINDOWS ==============================
compute_tf_windows <- function(tf_sites, window) {
  centers <- round(start(tf_sites) + (width(tf_sites) - 1) / 2)
  chr <- as.character(seqnames(tf_sites))
  chr_len <- seqlengths(HG19_GENOME)[chr]
  chr_len[is.na(chr_len)] <- max(seqlengths(HG19_GENOME), na.rm = TRUE)

  start_pos <- pmax(1L, centers - window)
  end_pos   <- pmin(centers + window, chr_len)

  data.table(chr = chr, start=start_pos, end=end_pos)
}

# ======================== COVERAGE BINS ========================
matrix_coverage <- function(cov, tf_windows, bin_size, positions) {

  n_sites <- nrow(tf_windows)
  n_bins  <- length(positions)

  M <- matrix(0, nrow=n_sites, ncol=n_bins)

  for (i in seq_len(n_sites)) {
    chr <- tf_windows$chr[i]
    chr_cov <- cov[[chr]]
    if (is.null(chr_cov)) next

    s <- tf_windows$start[i]
    e <- tf_windows$end[i]

    if (s > length(chr_cov)) next
    e <- min(e, length(chr_cov))
    if (s >= e) next

    full_width <- e - s + 1L
    full_bins  <- full_width %/% bin_size
    if (full_bins == 0) next

    bin_starts <- s + (0:(full_bins - 1)) * bin_size
    bin_ends   <- bin_starts + bin_size - 1

    v <- Views(chr_cov, IRanges(bin_starts, bin_ends))
    bin_sums <- as.numeric(viewSums(v))

    M[i, 1:full_bins] <- bin_sums
  }

  return(M)
}

# ======================== MAIN ================================

frag_paths <- sort(list.files(
  fragDir,
  pattern = "\\.(rds|starch)$",
  full.names = TRUE,
  ignore.case = TRUE
))

if (index < 1 || index > length(frag_paths))
  stop("Index out of range. Found: ", length(frag_paths))

curr_path <- frag_paths[index]
sample_id <- gsub("\\.(rds|starch)$", "", basename(curr_path))

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Output file (for combined mode)
outfile_combined <- file.path(outDir, paste0(sample_id, "_tf_footprint_matrix_cache.rds"))

# Output folder (for per-TF mode)
sample_out_dir <- file.path(outDir, sample_id)

cat("\n============================================================\n")
cat(" TF Footprint Matrix Precomputation (v2)\n")
cat(" Sample:   ", sample_id, "\n")
cat(" Mode:     ", save_mode, "\n")
cat(" Platform: ", platform, "\n")
cat("============================================================\n\n")

# -------- Skip processed sample --------------------------------
if (save_mode == "combined" && file.exists(outfile_combined)) {
  cat("✓ Combined output exists, skipping.\n")
  quit(save = "no")
}

if (save_mode == "per_tf" && dir.exists(sample_out_dir)) {
  cat("✓ Per-TF folder already exists, skipping.\n")
  quit(save = "no")
}

# -------- Load fragments --------------------------------------
cat("Step 1: Load fragments\n")
frag <- load_fragments_preprocessed(curr_path, frag_min, frag_max)

# -------- Coverage --------------------------------------------
cat("Step 2: Weighted coverage\n")
cov <- coverage(frag, weight = mcols(frag)$weight)

# -------- Load TFs --------------------------------------------
cat("Step 3: Discover TF files\n")
tf_files <- list.files(tfbsDir, pattern="^tf_.*_5k\\.bed$")
tf_names <- tolower(gsub("^tf_(.*)_5k\\.bed$", "\\1", tf_files))

# -------- Precompute positions --------------------------------
positions <- seq(-window, window, by=bin_size)

# -------- Process TFs -----------------------------------------
cat("Step 4: Processing TFs\n")

TF_matrices <- list()

for (i in seq_along(tf_names)) {

  tf <- tf_names[i]
  if (i %% 50 == 0) cat("  [", i, "/", length(tf_names), "] ", tf, "\n", sep="")

  tf_sites <- tryCatch(load_tf(tf, tfbsDir), error=function(e) NULL)
  if (is.null(tf_sites) || !length(tf_sites)) next

  tf_windows <- compute_tf_windows(tf_sites, window)
  M <- matrix_coverage(cov, tf_windows, bin_size, positions)

  TF_matrices[[tf]] <- M
}

# -------- Metadata --------------------------------------------
metadata <- list(
  sample_id    = sample_id,
  platform     = platform,
  window       = window,
  bin_size     = bin_size,
  frag_range   = c(frag_min, frag_max),
  tf_count     = length(TF_matrices),
  positions    = positions,
  creation_time = Sys.time()
)

# ===================================================================
# SAVE MODE OPTION
# ===================================================================

if (save_mode == "combined") {

  # ---------------- SAVE ONE LARGE FILE -------------------------
  output <- list(
    matrices = TF_matrices,
    metadata = metadata
  )

  saveRDS(output, outfile_combined)
  cat("✓ Saved combined cache: ", outfile_combined, "\n")

} else if (save_mode == "per_tf") {

  # ---------------- SAVE PER-TF FILES ---------------------------
  dir.create(sample_out_dir, recursive=TRUE, showWarnings=FALSE)

  # Save metadata
  saveRDS(metadata, file.path(sample_out_dir, "metadata.rds"))

  # Save each TF matrix separately
  for (tf in names(TF_matrices)) {
    outfile_tf <- file.path(sample_out_dir, paste0("tf_", tf, ".rds"))
    saveRDS(TF_matrices[[tf]], outfile_tf)
  }

  cat("✓ Saved ", length(TF_matrices), " TF matrices to ", sample_out_dir, "\n", sep="")
}

cat("============================================================\n\n")
