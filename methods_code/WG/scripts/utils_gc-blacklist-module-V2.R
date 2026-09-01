#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(data.table)
})

options(jupyter.rich_display = FALSE)

# ================================================================
# CONSTANTS / PATHS
# ================================================================
HG19_GENOME <- BSgenome.Hsapiens.UCSC.hg19

BLACKLIST_FILE    <- file.path("PlasmaTools", "data", "filters.hg19.rda")
GC_TARGET_HISEQ   <- file.path("reference_files", "target-hiseq.tsv")
GC_TARGET_NOVASEQ <- file.path("reference_files", "target-novaseq.tsv")

if (!file.exists(BLACKLIST_FILE)) {
  stop("Blacklist file not found: ", BLACKLIST_FILE)
}
load(BLACKLIST_FILE)  # loads filters.hg19
if (!exists("filters.hg19")) {
  stop("filters.hg19 object not found after loading blacklist")
}

BEDOPS_HELPER <- "/dcs04/scharpf/data/nvulpesc/tools/bedops_R.R"
if (!file.exists(BEDOPS_HELPER)) {
  stop("BEDOPS helper missing: ", BEDOPS_HELPER)
}
source(BEDOPS_HELPER)
if (!exists("unstarch_to_granges"))
  stop("Function unstarch_to_granges() missing after sourcing BEDOPS helper")

# ================================================================
# FUNCTIONS
# ================================================================

select_gc_target <- function(platform) {
  platform <- tolower(platform)
  switch(
    platform,
    "hiseq"   = GC_TARGET_HISEQ,
    "novaseq" = GC_TARGET_NOVASEQ,
    stop("platform must be 'hiseq' or 'novaseq'")
  )
}

# ----------------------------------------------------------------
# Load fragments (NO fragment-length filter)
# ----------------------------------------------------------------
load_fragments <- function(fp) {
  cat("Loading:", basename(fp), "\n"); flush.console()

  frag <- if (grepl("\\.starch$", fp, ignore.case = TRUE)) {
    cat("  Detected .starch format\n")
    unstarch_to_granges(fp, tmpdir = tempdir())
  } else if (grepl("\\.bed$", fp, ignore.case = TRUE)) {
    import(fp, format = "BED")
  } else {
    readRDS(fp)
  }

  initial_n <- length(frag)
  cat("  Initial fragments:", initial_n, "\n"); flush.console()

  seqlevelsStyle(frag) <- "UCSC"

  keep_chrs <- paste0("chr", c(1:22, "X", "Y"))
  frag <- frag[seqnames(frag) %in% keep_chrs]
  post_chr <- length(frag)
  cat("  After chr filter:", post_chr, "\n"); flush.console()

  if (!length(frag)) stop("No fragments remain after chr filtering")

  list(
    frag = frag,
    qc   = list(initial = initial_n, post_chrom = post_chr)
  )
}

# ----------------------------------------------------------------
# Blacklist filter
# ----------------------------------------------------------------
apply_blacklist_filter <- function(frag, blacklist) {
  cat("  Applying blacklist filter...\n"); flush.console()
  hits <- findOverlaps(frag, blacklist, select = "first")
  keep <- is.na(hits)
  frag <- frag[keep]

  cat("  After blacklist:", length(frag), "fragments\n"); flush.console()

  list(frag = frag, n = length(frag))
}

# ----------------------------------------------------------------
# GC CORRECTION — memory-optimized implementation
# ----------------------------------------------------------------
apply_gc_correction <- function(frag, GC_TARGET) {

  cat("  Computing GC content (memory-optimized)...\n"); flush.console()

  # Preallocate GC fraction vector
  gc_frac <- numeric(length(frag))

  # Process chromosome-by-chromosome to avoid a huge DNAStringSet
  chrs <- sort(unique(as.character(seqnames(frag))))

  for (chr in chrs) {
    idx_chr <- which(as.character(seqnames(frag)) == chr)
    f_chr   <- frag[idx_chr]

    # Extract chromosome sequence once (no copying)
    chr_seq <- HG19_GENOME[[chr]]

    # Views are zero-copy windows on the chromosome sequence
    v <- Views(chr_seq, start = start(f_chr), end = end(f_chr))

    # GC count per fragment: G + C
    gc_counts_chr <- rowSums(letterFrequency(v, c("G", "C")))
    gc_frac[idx_chr] <- gc_counts_chr / width(f_chr)
  }

  # GC fraction range filter
  keep_gc <- gc_frac >= 0.20 & gc_frac <= 0.80
  if (!any(keep_gc)) stop("All fragments removed by GC-range filtering")

  frag    <- frag[keep_gc]
  gc_frac <- gc_frac[keep_gc]

  gc_bin <- round(gc_frac, 2)

  dt <- data.table(
    chr = as.character(seqnames(frag)),
    gc  = gc_bin,
    idx = seq_along(frag)
  )

  # Observed GC distribution
  obs <- dt[, .N, by = .(chr, gc)]
  setnames(obs, "N", "n")

  # Merge expected vs observed
  setkey(obs, chr, gc)
  setkey(GC_TARGET, chr, gc)

  merged <- GC_TARGET[obs, nomatch = NA]
  merged <- merged[!is.na(target)]
  merged[, weight := target / n]

  # Join weights to fragments
  setkey(dt, chr, gc)
  setkey(merged, chr, gc)

  weighted_dt <- merged[dt]
  keep_final  <- !is.na(weighted_dt$weight)

  if (!any(keep_final))
    stop("All fragments lost after GC weighting — check GC target file")

  keep_idx  <- weighted_dt$idx[keep_final]
  frag      <- frag[keep_idx]
  w         <- weighted_dt$weight[keep_final]
  gc_final  <- weighted_dt$gc[keep_final]

  mcols(frag)$gc_bin  <- gc_final
  mcols(frag)$weight  <- w

  cat("  After GC correction:", length(frag), "fragments retained\n"); flush.console()

  list(
    frag            = frag,
    n_after_gc      = sum(keep_gc),
    n_after_weight  = length(frag),
    gc_bins         = length(unique(gc_final)),
    weight_summary  = summary(w)
  )
}

# ================================================================
# MAIN — index-aware
# ================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript utils_gc-blacklist-module-V2.R <index> <fragDir> <outDir> <platform> <frag_min> <frag_max>")
}

index     <- as.numeric(args[1])
fragDir   <- args[2]
outDir    <- args[3]
platform  <- args[4]
frag_min  <- as.numeric(args[5])   # unused but kept for interface consistency
frag_max  <- as.numeric(args[6])   # unused but kept for interface consistency

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Load blacklist
blacklist <- reduce(filters.hg19)

# Load GC target
gc_target_file <- select_gc_target(platform)
GC_TARGET <- fread(gc_target_file)

# Normalize GC_TARGET (matches your screenshot: seqnames / gc / gcmed)
tgt_cols <- names(GC_TARGET)
if ("seqnames" %in% tgt_cols) setnames(GC_TARGET, "seqnames", "chr")
if ("chrom"    %in% tgt_cols) setnames(GC_TARGET, "chrom",    "chr")
if ("gcmed"    %in% tgt_cols) setnames(GC_TARGET, "gcmed",    "target")

GC_TARGET[, chr    := as.character(chr)]
GC_TARGET[, gc     := round(as.numeric(gc), 2)]
GC_TARGET[, target := as.numeric(target)]
GC_TARGET <- GC_TARGET[gc >= 0.20 & gc <= 0.80]

# Get all fragment files
frag_files <- list.files(
  fragDir,
  pattern     = "\\.(rds|starch|bed)$",
  full.names  = TRUE,
  ignore.case = TRUE
)

frag_files <- sort(frag_files)

if (length(frag_files) == 0)
  stop("No fragment files found in: ", fragDir)

if (index < 1 || index > length(frag_files))
  stop("Index ", index, " out of range (1:", length(frag_files), ")")

infile      <- frag_files[index]
sample_name <- gsub("\\.(rds|starch|bed)$", "", basename(infile), ignore.case = TRUE)

outfile <- file.path(outDir, paste0(sample_name, ".gc_corrected.rds"))
if (file.exists(outfile)) {
  cat("✓ Output exists, skipping sample:", sample_name, "\n")
  quit(status = 0, save="no")
}

cat("\n==============================\n")
cat("QC + GC Correction for sample index:", index, "\n")
cat("Sample:", sample_name, "\n")
cat("File:", infile, "\n")
cat("==============================\n\n")
flush.console()

# ================================================================
# PROCESSING STEPS
# ================================================================
x    <- load_fragments(infile)
frag <- x$frag
qc   <- x$qc
gc()

bl   <- apply_blacklist_filter(frag, blacklist)
frag <- bl$frag
qc$after_blacklist <- bl$n
gc()

gc_out <- apply_gc_correction(frag, GC_TARGET)
frag   <- gc_out$frag

qc$after_gc_filter  <- gc_out$n_after_gc
qc$after_gc_weight  <- gc_out$n_after_weight
qc$gc_bins          <- gc_out$gc_bins
qc$weight_summary   <- capture.output(gc_out$weight_summary)
gc()

# ================================================================
# SAVE CLEANED FRAGMENTS
# ================================================================
outfile <- file.path(outDir, paste0(sample_name, ".gc_corrected.rds"))
saveRDS(frag, outfile)

# ================================================================
# SAVE QC TEXT REPORT
# ================================================================
qcfile <- file.path(outDir, paste0(sample_name, ".qc.txt"))
conn <- file(qcfile, "wt")

writeLines(c(
  "==================== SAMPLE QC REPORT ====================",
  paste("Sample:", sample_name),
  paste("Input file:", infile),
  paste("Platform:", platform),
  "",
  "------ Fragment Counts ------",
  paste("Initial fragments:",        qc$initial),
  paste("After chromosome filter:", qc$post_chrom),
  paste("After blacklist:",         qc$after_blacklist),
  paste("After GC range filter:",   qc$after_gc_filter),
  paste("After GC weight assignment:", qc$after_gc_weight),
  "",
  "------ GC Information ------",
  paste("GC bins represented:", qc$gc_bins),
  "",
  "------ Weight Summary ------",
  qc$weight_summary,
  "",
  "------ Width Summary ------",
  paste(capture.output(summary(width(frag))), collapse = "\n"),
  "",
  "==================== END QC REPORT ======================="
), conn)

close(conn)

cat("\n✓ Saved GC-corrected file:", outfile)
cat("\n✓ Saved QC report:", qcfile)
cat("\n✓ Final fragment count:", length(frag), "\n\n")
flush.console()
