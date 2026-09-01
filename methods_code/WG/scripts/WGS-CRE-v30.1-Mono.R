#!/usr/bin/env Rscript
# ===================================================================
# cfDNA CRE Fragmentomics Pipeline (v30.1 – GC-aware, mode-dependent)
#  - Modes:
#       "full"         : blacklist + fragment-level GC correction
#       "nogc":         : blacklist only, unit weights
#       "preprocessed" : skip blacklist & GC; use precomputed weights
#  - Weighted coverage & weighted hit counts for all features
#  - Summit vs Range–aware (per-file mode)
#  - Multi-directory CRE discovery (.rds > .bed dedup)
#  - Deterministic sample order; directory-aware CRE names
#  - Endmotif helper; flank-collapse warning
#  - Saves nested results + optional flat table
# ===================================================================

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(rtracklayer)
  library(GenomicRanges)
  library(Biostrings)
  library(transport)
  library(IRanges)
  library(yaml)
  library(zoo)
  library(data.table)
})

# Genome reference for hg19
HG19_GENOME <- BSgenome.Hsapiens.UCSC.hg19

# ============================ REFERENCES ============================
# Assume this script is run from scripts/ directory
BLACKLIST_FILE    = file.path("PlasmaTools", "data", "filters.hg19.rda")
GC_TARGET_HISEQ   = file.path("reference_files", "target-hiseq.tsv")
GC_TARGET_NOVASEQ = file.path("reference_files", "target-novaseq.tsv")

if (!file.exists(BLACKLIST_FILE)) {
  stop("Blacklist file not found: ", BLACKLIST_FILE)
}
load(BLACKLIST_FILE)  # loads `filters.hg19` (GRanges)

if (!exists("filters.hg19")) {
  stop("filters.hg19 object not found after loading ", BLACKLIST_FILE)
}

# Absolute path to bedops helper (for .starch support)
BEDOPS_HELPER <- "/dcs04/scharpf/data/nvulpesc/tools/bedops_R.R"
if (!file.exists(BEDOPS_HELPER)) {
  stop("BEDOPS helper script not found at: ", BEDOPS_HELPER)
}
source(BEDOPS_HELPER)
if (!exists("unstarch_to_granges")) {
  stop("Function `unstarch_to_granges()` not found after sourcing BEDOPS helper: ", BEDOPS_HELPER)
}

select_gc_reference <- function(platform) {
  platform <- tolower(platform)
  if (platform == "hiseq") {
    if (!file.exists(GC_TARGET_HISEQ))
      stop("HiSeq GC target file not found: ", GC_TARGET_HISEQ)
    GC_TARGET_HISEQ
  } else if (platform == "novaseq") {
    if (!file.exists(GC_TARGET_NOVASEQ))
      stop("NovaSeq GC target file not found: ", GC_TARGET_NOVASEQ)
    GC_TARGET_NOVASEQ
  } else {
    stop("Invalid platform: ", platform, " (must be 'hiseq' or 'novaseq')")
  }
}

# ============================ ARGS =============================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11)
  stop("Usage: Rscript WGS-CRE-v30.1.R index fragDir creDirs outDir center_bp flank_left flank_right frag_min frag_max mode platform [enabled_modules]
  • mode ∈ {full, nogc, preprocessed}")

index      <- as.numeric(args[1])
fragDir    <- args[2]
creDirs    <- args[3]
outDir     <- args[4]
center_bp  <- as.numeric(args[5])
flank_left <- args[6]   # e.g., "-3000:-2750"
flank_right<- args[7]   # e.g., "2750:3000"
frag_min   <- as.numeric(args[8])
frag_max   <- as.numeric(args[9])
mode       <- tolower(args[10])  # "full", "nogc", or "preprocessed"
platform   <- tolower(args[11])
enabled_arg<- ifelse(length(args) >= 12, args[12], "all")

mode <- match.arg(mode, c("full","nogc","preprocessed"))

gc_correction   <- (mode == "full")
use_preprocessed<- (mode == "preprocessed")

if (!platform %in% c("hiseq","novaseq")) {
  stop("Invalid platform: ", platform, " (must be 'hiseq' or 'novaseq')")
}

default_modules <- c("relcov","swfse","fld","qchasm","endmotif","qchasm_outward","halo")
enabled_modules <- if (tolower(enabled_arg) == "all") {
  default_modules
} else {
  intersect(strsplit(enabled_arg, ",")[[1]], default_modules)
}
if (length(enabled_modules) == 0)
  stop("No valid modules enabled. Choose from: ", paste(default_modules, collapse = ", "))

cat("\nEnabled modules:", paste(enabled_modules, collapse = ", "), "\n")
cat("Mode:", mode, "\n")
cat("Fragment-level GC correction:", gc_correction, " (platform:", platform, ")\n")
cat("Use preprocessed fragments:", use_preprocessed, "\n")

# ============================ CONSTANTS ========================
EPSILON <- 1e-8
MIN_FRAGMENTS_FOR_STATS <- 10L
VALID_CHROMOSOMES <- paste0("chr", c(1:22,"X","Y","M"))

# ============================ HELPERS ==========================
parse_window <- function(win) {
  parts <- as.numeric(strsplit(gsub("\\s","",win),":")[[1]])
  if (length(parts) != 2 || any(!is.finite(parts)))
    stop("Invalid window: ", win, " (expected like -300:-200)")
  if (parts[1] > parts[2]) parts <- rev(parts)
  parts
}

safe_import_cre <- function(fp) {
  if (!file.exists(fp)) {
    cat("Warning: CRE file missing:", fp, "\n")
    return(GRanges())
  }
  gr <- tryCatch({
    if (grepl("\\.rds$", fp, ignore.case = TRUE)) {
      readRDS(fp)
    } else {
      import(fp, format = "BED")
    }
  }, error = function(e) {
    cat("Warning: import failed:", fp, ":", e$message, "\n")
    return(GRanges())
  })

  if (length(gr) > 0) {
    suppressWarnings({
      seqlevelsStyle(gr) <- "UCSC"
      common_seqs <- intersect(seqlevels(gr), seqlevels(HG19_GENOME))
      if (length(common_seqs) > 0) {
        gr <- keepSeqlevels(gr, common_seqs, pruning.mode = "coarse")
        seqinfo(gr) <- keepSeqlevels(seqinfo(HG19_GENOME), seqlevels(gr), pruning.mode = "coarse")
        gr <- trim(gr)
      }
    })
  }
  gr
}

harmonize_to_frag <- function(cre.gr, frag.gr) {
  if (!length(cre.gr)) return(cre.gr)
  suppressWarnings(seqlevelsStyle(cre.gr) <- seqlevelsStyle(frag.gr))
  keep <- intersect(seqlevels(cre.gr), seqlevels(frag.gr))
  if (!length(keep)) return(GRanges())
  cre.gr <- keepSeqlevels(cre.gr, keep, pruning.mode = "coarse")
  suppressWarnings(seqinfo(cre.gr) <- seqinfo(frag.gr)[seqlevels(cre.gr)])
  trim(cre.gr[width(cre.gr) > 0])
}

safe_center_window <- function(gr, width_bp) {
  if (!length(gr)) return(GRanges())
  if (width_bp %% 2 == 0) width_bp <- width_bp + 1L
  trim(resize(gr, width = width_bp, fix = "center"))
}

safe_flank_window <- function(gr, fr) {
  if (!length(gr)) return(GRanges())
  w   <- abs(diff(fr)) + 1L
  off <- round(mean(fr))
  trim(GenomicRanges::shift(resize(gr, width = w, fix = "center"), off))
}

total_flank_width <- function(w) {
  flankL <- if (length(w$flankL)) w$flankL else GRanges()
  flankR <- if (length(w$flankR)) w$flankR else GRanges()
  all_flanks <- c(flankL, flankR)
  if (!length(all_flanks)) return(0L)
  flank_union <- reduce(all_flanks)
  sum(width(flank_union))
}

validate_windows <- function(w, p) {
  issues <- character(0)
  if (length(w$center) && length(w$flankL) && any(overlapsAny(w$center, w$flankL)))
    issues <- c(issues, "Center overlaps left flank")
  if (length(w$center) && length(w$flankR) && any(overlapsAny(w$center, w$flankR)))
    issues <- c(issues, "Center overlaps right flank")
  if (length(issues)) {
    cat("⚠️  Window validation warnings:\n")
    for (i in issues) cat("   -", i, "\n")
  }
}

cache_overlaps <- function(frag.gr, w) {
  list(
    center = if (length(w$center)) overlapsAny(frag.gr, w$center) else logical(length(frag.gr)),
    flankL = if (length(w$flankL)) overlapsAny(frag.gr, w$flankL) else logical(length(frag.gr)),
    flankR = if (length(w$flankR)) overlapsAny(frag.gr, w$flankR) else logical(length(frag.gr))
  )
}

mean_profile_across_windows <- function(cov, gr, width_bp, chunk = 5000L) {
  if (!length(gr)) return(rep(NA_real_, width_bp))
  if (width_bp %% 2 == 0) width_bp <- width_bp + 1L
  sum_profile <- numeric(width_bp)
  n <- 0L

  for (chr in intersect(names(cov), as.character(seqlevels(gr)))) {
    g <- gr[seqnames(gr) == chr]
    if (!length(g)) next
    v   <- cov[[chr]]
    rng <- ranges(g)

    idx <- split(seq_along(g), ceiling(seq_along(g) / chunk))
    for (ids in idx) {
      starts <- start(rng)[ids]
      ends   <- end(rng)[ids]
      starts <- pmax(1L, starts)
      ends   <- pmin(length(v), ends)
      vw     <- Views(v, start = starts, end = ends)
      m      <- as.matrix(vw)
      if (!is.matrix(m) || ncol(m) != width_bp) next
      sum_profile <- sum_profile + colSums(m, na.rm = TRUE)
      n <- n + nrow(m)
    }
  }

  if (!n) return(rep(NA_real_, width_bp))
  sum_profile / n
}

cov_profile_viewmeans <- function(cov, gr, width_bp) {
  if (!length(gr)) return(NA_real_)
  if (width_bp %% 2 == 0) width_bp <- width_bp + 1L
  mu <- mean_profile_across_windows(cov, gr, width_bp)
  if (length(mu) != width_bp || all(!is.finite(mu))) return(NA_real_)
  half <- width_bp %/% 2
  names(mu) <- as.character(seq.int(-half, half))
  mu
}

is_summit_file <- function(gr) {
  length(gr) > 0 && all(width(gr) == 1L)
}

# -------- Blacklist-only helper (for mode = "nogc") --------
apply_blacklist_only <- function(frag.gr, blacklist) {
  cat("  - Removing fragments overlapping blacklist regions (mode='nogc')...\n")
  keep <- !overlapsAny(frag.gr, blacklist)
  frag.gr <- frag.gr[keep]
  cat("    Remaining fragments after blacklist filter:", length(frag.gr), "\n")
  if (!length(frag.gr)) stop("No fragments remain after blacklist filtering.")
  frag.gr
}

# -------- Endmotif helper --------
get_5prime_kmer_intervals <- function(gr, k = 4) {
  is_minus <- as.logical(strand(gr) == "-")
  s <- start(gr)
  e <- end(gr)
  s5 <- ifelse(is_minus, pmax(s, e - k + 1L), s)
  e5 <- ifelse(is_minus, e, pmin(e, s + k - 1L))
  GRanges(seqnames = seqnames(gr),
          ranges   = IRanges(start = s5, end = e5),
          strand   = strand(gr))
}

# -------- Fragment-level GC correction (Cristiano-style) --------
apply_fragment_gc_correction <- function(
  frag.gr,
  genome,
  blacklist,
  target_file,
  by_chrom = TRUE
) {

  cat("  - Removing fragments overlapping blacklist regions...\n")
  keep <- !overlapsAny(frag.gr, blacklist)
  frag.gr <- frag.gr[keep]
  cat("    Remaining fragments after blacklist filter:", length(frag.gr), "\n")
  if (!length(frag.gr)) stop("No fragments remain after blacklist filtering.")

  # --------------------------------------------------------
  # Compute GC fraction per fragment
  # --------------------------------------------------------
  cat("  - Computing per-fragment GC fraction...\n")
  seqs <- getSeq(genome, frag.gr)
  gc_mat <- letterFrequency(seqs, letters = c("G","C"), as.prob = FALSE)
  gc_counts <- rowSums(gc_mat)
  gc_frac <- gc_counts / width(frag.gr)

  # Cristiano range filter (0.20–0.80)
  keep_gc <- !is.na(gc_frac) & gc_frac >= 0.20 & gc_frac <= 0.80
  frag.gr <- frag.gr[keep_gc]
  gc_frac <- gc_frac[keep_gc]

  if (!length(frag.gr)) stop("No fragments remain after GC range filtering.")
    
  mcols(frag.gr)$gc_frac <- gc_frac
  gc_bin <- round(gc_frac, 2)
  mcols(frag.gr)$gc_bin <- gc_bin

  # --------------------------------------------------------
  # Load and normalize GC target table
  # --------------------------------------------------------
  cat("  - Loading GC target table:", target_file, "\n")
  target <- fread(target_file)

  # Normalize column names from Cristiano reference file
  if ("seqnames" %in% names(target)) setnames(target, "seqnames", "chr")
  if ("chrom"    %in% names(target)) setnames(target, "chrom",    "chr")
  if ("gcmed"    %in% names(target)) setnames(target, "gcmed",    "target")

  if (!all(c("chr","gc","target") %in% names(target))) {
    stop("GC target file must contain columns 'chr', 'gc', 'target'.")
  }

  # Numeric normalization
  target[, chr    := as.character(chr)]
  target[, gc     := round(as.numeric(gc), 2)]
  target[, target := as.numeric(target)]

  # Keep GC bins in Cristiano range
  target <- target[gc >= 0.20 & gc <= 0.80]

  # Handle duplicate GC bins after rounding
  dups <- target[, .N, by = .(chr, gc)][N > 1]
  if (nrow(dups) > 0) {
    cat("⚠️  Duplicate GC bins detected — merging\n")
    target <- target[, .(target = sum(target, na.rm = TRUE)), by = .(chr, gc)]
  }

  # --------------------------------------------------------
  # Build fragment GC table
  # --------------------------------------------------------
  dt <- data.table(
    chr = as.character(seqnames(frag.gr)),
    gc  = gc_bin
  )

  # Normalize chr prefix if necessary
  if (any(grepl("^chr", dt$chr)) && !any(grepl("^chr", target$chr))) {
    target[, chr := paste0("chr", chr)]
  }

  # --------------------------------------------------------
  # Observed GC-bin frequencies
  # --------------------------------------------------------
  obs <- dt[, .N, by = .(chr, gc)]
  setnames(obs, "N", "n")

  # Merge observed vs expected
  merged <- merge(obs, target, by = c("chr","gc"), all.x = TRUE)
  merged <- merged[!is.na(target)]

  if (!nrow(merged))
    stop("No overlapping GC bins between fragment data and GC target.")

  # Weight = expected / observed
  merged[, weight := target / n]

  # --------------------------------------------------------
  # Assign per-fragment weights
  # --------------------------------------------------------
  key_obs <- paste(dt$chr, dt$gc)
  key_ref <- paste(merged$chr, merged$gc)
  idx_map <- match(key_obs, key_ref)
  weight_vec <- merged$weight[idx_map]

  # Drop fragments with NA weights
  keep_w <- !is.na(weight_vec)
  frag.gr <- frag.gr[keep_w]
  weight_vec <- weight_vec[keep_w]
  gc_bin <- gc_bin[keep_w]

  if (!length(frag.gr))
    stop("All fragments lost after weight assignment.")

  mcols(frag.gr)$gc_bin <- gc_bin
  mcols(frag.gr)$weight <- weight_vec

  cat("    Fragments retained after GC weighting:", length(frag.gr), "\n")
  return(frag.gr)
}

# ----------------------------------------------------------------------
# Load fragments (raw or preprocessed)
# ----------------------------------------------------------------------
load_fragments <- function(fp, frag_min, frag_max, use_preprocessed = FALSE) {
  cat("Loading fragments:", basename(fp), "\n")

  frag.gr <- if (grepl("\\.starch$", fp, ignore.case = TRUE)) {
    cat("  Detected .starch input\n")
    unstarch_to_granges(fp, tmpdir = tempdir())
  } else {
    readRDS(fp)
  }

  cat("  Initial fragment count:", length(frag.gr), "\n")

  # Seqlevel normalization
  suppressWarnings({
    seqlevelsStyle(frag.gr) <- "UCSC"
    common <- intersect(seqlevels(frag.gr), VALID_CHROMOSOMES)
    if (length(common) > 0) {
      frag.gr <- keepSeqlevels(frag.gr, common, pruning.mode = "coarse")
      seqinfo(frag.gr) <- keepSeqlevels(seqinfo(HG19_GENOME), seqlevels(frag.gr), pruning.mode = "coarse")
    }
  })

  # Chromosome filter
  frag.gr <- frag.gr[seqnames(frag.gr) %in% VALID_CHROMOSOMES]
  cat("  After chr filter:", length(frag.gr), "\n")

  # Fragment length filtering
  frag.gr <- frag.gr[width(frag.gr) >= frag_min & width(frag.gr) <= frag_max]
  cat("  After size filter:", length(frag.gr), "\n")

  if (!length(frag.gr)) stop("No fragments remain after filtering.")

  # If preprocessed mode → require weight column
  if (use_preprocessed && !("weight" %in% names(mcols(frag.gr)))) {
    stop("use_preprocessed=TRUE, but this fragment file lacks a 'weight' column.")
  }

  return(frag.gr)
}

# ======================= FEATURE FUNCTIONS (from v20 → v30.1) =======================

compute_rel_cov <- function(frag.gr, wins, cov, params, overlaps_cache = NULL) {
  W <- params$center
  width_bp <- 2 * params$Wmax + 1L
  mu <- cov_profile_viewmeans(cov, wins$profile, width_bp)

  rel_cov_profile <- NA_real_
  if (is.numeric(mu) && any(is.finite(mu))) {
    x_axis <- as.integer(names(mu))
    idx_center <- which(abs(x_axis) <= W)
    idx_flank  <- which((x_axis >= params$flank_left[1] & x_axis <= params$flank_left[2]) |
                        (x_axis >= params$flank_right[1] & x_axis <= params$flank_right[2]))
    if (length(idx_center) > 0 && length(idx_flank) > 0) {
      rel_cov_profile <- mean(mu[idx_center], na.rm = TRUE) /
        (mean(mu[idx_flank], na.rm = TRUE) + EPSILON)
    }
  }

  weights <- if ("weight" %in% names(mcols(frag.gr))) mcols(frag.gr)$weight else rep(1, length(frag.gr))

  if (!is.null(overlaps_cache)) {
    c_hits <- sum(weights[overlaps_cache$center])
    flank_idx <- overlaps_cache$flankL | overlaps_cache$flankR
    f_hits <- sum(weights[flank_idx])
  } else {
    c_log <- overlapsAny(frag.gr, wins$center)
    flank_log <- overlapsAny(frag.gr, wins$flankL) | overlapsAny(frag.gr, wins$flankR)
    c_hits <- sum(weights[c_log])
    f_hits <- sum(weights[flank_log])
  }

  # Summit — center is ±bp around summits
  if (!length(wins$center)) {
    center_bp_total <- 0L
  } else if (all(width(wins$center) == 1L)) {
    center_bp_total <- length(wins$center) * (2 * params$center + 1L)
  } else {
    # Range-mode CRE widths
    center_bp_total <- sum(width(wins$center))
    if (center_bp_total <= 0) {
      center_bp_total <- length(wins$center) * mean(width(wins$center))
    }
  }

  flank_bp_total <- total_flank_width(wins)

  rel_cov_density <- if (center_bp_total > 0 && flank_bp_total > 0) {
    (c_hits / center_bp_total) / ((f_hits / flank_bp_total) + EPSILON)
  } else NA_real_

  setNames(c(rel_cov_profile, rel_cov_density),
           c("rel_cov_profile", "rel_cov_density"))
}


compute_swfse <- function(frag.gr, wins, params, overlaps_cache = NULL) {
  win  <- params$frag_window
  step <- params$frag_step
  bins <- seq(params$frag_min, params$frag_max - win, by = step)
  res  <- numeric(length(bins))
  names(res) <- paste0("SWFSE_w", bins, "-", bins + win, "bp")

  flank_width_total <- total_flank_width(wins)

  if (!length(wins$center)) {
    center_bp_total <- 0L
  } else if (all(width(wins$center) == 1L)) {
    center_bp_total <- length(wins$center) * (2 * params$center + 1L)
  } else {
    center_bp_total <- sum(width(wins$center))
    if (center_bp_total <= 0) {
      center_bp_total <- length(wins$center) * mean(width(wins$center))
    }
  }

  if (flank_width_total <= 0 || !length(wins$center) || center_bp_total <= 0) {
    res[] <- NA_real_
  } else {
    weights_all <- if ("weight" %in% names(mcols(frag.gr))) mcols(frag.gr)$weight else rep(1, length(frag.gr))

    for (i in seq_along(bins)) {
      r1 <- bins[i]
      r2 <- bins[i] + win
      sel <- width(frag.gr) >= r1 & width(frag.gr) < r2
      if (!any(sel)) {
        res[i] <- NA_real_
        next
      }
      frag_subset <- frag.gr[sel]
      w_sub <- weights_all[sel]

      c_hits <- sum(w_sub[overlapsAny(frag_subset, wins$center)])
      f_hits <- sum(w_sub[overlapsAny(frag_subset, wins$flankL)]) +
                sum(w_sub[overlapsAny(frag_subset, wins$flankR)])

      res[i] <- log2((c_hits / (center_bp_total) + EPSILON) /
                     (f_hits / (flank_width_total) + EPSILON))
    }
  }

  mids <- bins + win / 2
  if (all(is.na(res))) {
    res_sum <- c(NA_real_, NA_real_, NA_real_)
  } else {
    rng <- max(res, na.rm = TRUE) - min(res, na.rm = TRUE)
    res_norm <- if (rng > 0) (res - min(res, na.rm = TRUE)) / rng else rep(0, length(res))
    slope <- suppressWarnings(tryCatch(coef(lm(res ~ mids))[2], error = function(e) NA_real_))
    ent   <- -sum(res_norm * log2(res_norm + EPSILON), na.rm = TRUE)
    res_sum <- c(mean(res, na.rm = TRUE), slope, ent)
  }

  names(res_sum) <- c(
    paste0("SWFSE_mean_", params$frag_min, "-", params$frag_max),
    paste0("SWFSE_slope_", params$frag_min, "-", params$frag_max),
    paste0("SWFSE_entropy_", params$frag_min, "-", params$frag_max)
  )

  c(res, res_sum)
}


compute_fld <- function(frag.gr, wins, params, overlaps_cache = NULL) {

  weighted_wasserstein <- function(x, wx, y, wy, p = 1) {
    wx <- wx / sum(wx)
    wy <- wy / sum(wy)
    ox <- order(x); x <- x[ox]; wx <- wx[ox]
    oy <- order(y); y <- y[oy]; wy <- wy[oy]
    Fx <- cumsum(wx)
    Fy <- cumsum(wy)
    qs <- sort(unique(c(Fx, Fy)))
    xq <- approx(Fx, x, qs, rule = 2)$y
    yq <- approx(Fy, y, qs, rule = 2)$y
    mean(abs(xq - yq)^p)^(1/p)
  }

  if (!is.null(overlaps_cache)) {
    c_sel <- frag.gr[overlaps_cache$center]
    f_sel <- frag.gr[overlaps_cache$flankL | overlaps_cache$flankR]
  } else {
    c_sel <- frag.gr[overlapsAny(frag.gr, wins$center)]
    f_sel <- c(frag.gr[overlapsAny(frag.gr, wins$flankL)],
               frag.gr[overlapsAny(frag.gr, wins$flankR)])
  }

  if (length(c_sel) < MIN_FRAGMENTS_FOR_STATS ||
      length(f_sel) < MIN_FRAGMENTS_FOR_STATS) {
    return(setNames(rep(NA_real_, 5L),
                    c("wasserstein1d", "wasserstein2d",
                      "FLD_skewness", "FLD_kurtosis", "FLD_entropy")))
  }

  wc <- mcols(c_sel)$weight
  wf <- mcols(f_sel)$weight

  c_frags <- width(c_sel)
  f_frags <- width(f_sel)

  w1 <- weighted_wasserstein(c_frags, wc, f_frags, wf, p = 1)
  w2 <- weighted_wasserstein(c_frags, wc, f_frags, wf, p = 2)

  pt <- wc / sum(wc)
  ent <- -sum(pt * log2(pt + EPSILON))
  m  <- weighted.mean(c_frags, wc)
  sdv <- sqrt(weighted.mean((c_frags - m)^2, wc))
  sdv <- ifelse(sdv == 0, EPSILON, sdv)

  sk <- weighted.mean((c_frags - m)^3, wc) / (sdv^3 + EPSILON)
  ku <- weighted.mean((c_frags - m)^4, wc) / (sdv^4 + EPSILON)

  setNames(c(w1, w2, sk, ku, ent),
           c("wasserstein1d", "wasserstein2d",
             "FLD_skewness", "FLD_kurtosis", "FLD_entropy"))
}


compute_qchasm <- function(frag.gr, wins, cov, params, bin = "global") {
  W <- params$center
  width_bp <- 2 * params$Wmax + 1L
  mu <- cov_profile_viewmeans(cov, wins$profile, width_bp)
  nm <- sprintf("qChasm_q25_±%dbp_%s", W, bin)

  if (!is.numeric(mu) || all(!is.finite(mu))) return(setNames(NA_real_, nm))

  x_axis <- as.integer(names(mu))
  idx_center <- which(abs(x_axis) <= W)
  idx_flank  <- which((x_axis >= params$flank_left[1] & x_axis <= params$flank_left[2]) |
                      (x_axis >= params$flank_right[1] & x_axis <= params$flank_right[2]))

  if (!length(idx_center) || !length(idx_flank)) return(setNames(NA_real_, nm))

  Tthr <- stats::quantile(mu[idx_flank], 0.25, na.rm = TRUE)
  deficit <- pmax(0, Tthr - mu[idx_center])

  qchasm <- sum(deficit, na.rm = TRUE) / max(1, sum(is.finite(mu[idx_center])))

  setNames(qchasm, nm)
}


compute_qchasm_outward <- function(wins, cov, params, bin = "global",
                                   flank_quantile = 0.25, recover_k = 10, smooth_k = 7) {
  W <- params$center
  width_bp <- 2 * params$Wmax + 1L
  mu <- cov_profile_viewmeans(cov, wins$profile, width_bp)
  nm_pref <- sprintf("qChasmOut_±%dbp_%s", W, bin)

  out <- setNames(rep(NA_real_, 8L), paste0(nm_pref, c(
    "_AUC", "_extent_left", "_extent_right", "_extent_total",
    "_symmetry", "_slope_in", "_slope_out", "_depth_z"
  )))

  if (!is.numeric(mu) || all(!is.finite(mu))) return(out)

  if (!is.na(smooth_k) && smooth_k >= 3 && (smooth_k %% 2 == 1)) {
    mu_sm <- stats::filter(mu, rep(1/smooth_k, smooth_k), sides = 2)
    mu_sm <- zoo::na.locf(mu_sm, na.rm = FALSE)
    mu_sm <- zoo::na.locf(mu_sm, na.rm = FALSE, fromLast = TRUE)
    mu[!is.finite(mu_sm)] <- median(mu, na.rm = TRUE)
  }

  x <- as.integer(names(mu))
  idx_center <- which(abs(x) <= W)
  idx_flankL <- which(x >= params$flank_left[1]  & x <= params$flank_left[2])
  idx_flankR <- which(x >= params$flank_right[1] & x <= params$flank_right[2])
  idx_flank <- sort(unique(c(idx_flankL, idx_flankR)))

  if (!length(idx_flank)) {
    n <- length(mu); pad <- ceiling(0.1 * n)
    idx_flank <- c(1:pad, (n - pad + 1):n)
  }

  mu_f <- mu[idx_flank]
  if (all(!is.finite(mu_f))) mu_f[] <- median(mu, na.rm = TRUE)

  base_level <- quantile(mu_f, flank_quantile, na.rm = TRUE)
  fl_med <- median(mu_f, na.rm = TRUE)
  fl_mad <- mad(mu_f, constant = 1, na.rm = TRUE)

  sweep_side <- function(mu, x, dir = c("left","right"), base_level, recover_k) {
    dir <- match.arg(dir)
    i0 <- which.min(abs(x))
    step <- if (dir == "left") -1L else 1L
    auc <- 0; consec <- 0; used <- integer(0)
    i <- i0
    repeat {
      i <- i + step
      if (i < 1 || i > length(mu)) break
      d <- max(0, base_level - mu[i])
      if (d > 0) {
        auc <- auc + d; consec <- 0; used <- c(used, i)
      } else {
        consec <- consec + 1
        if (consec >= recover_k) break
      }
    }
    extent <- if (length(used)) max(abs(x[used]), na.rm = TRUE) else 0
    list(AUC = auc, extent = extent, idx = used)
  }

  left  <- sweep_side(mu, x, "left",  base_level, recover_k)
  right <- sweep_side(mu, x, "right", base_level, recover_k)

  auc_total <- left$AUC + right$AUC
  extent_tot <- left$extent + right$extent

  symmetry <- if ((left$AUC + right$AUC) > 0) {
    abs(left$AUC - right$AUC) / (left$AUC + right$AUC)
  } else { NA_real_ }

  slope_window <- min(5L, max(length(left$idx), length(right$idx)))

  slope_in  <- if (length(left$idx)  >= slope_window) mean(diff(mu[sort(left$idx)[1:slope_window]])) else NA_real_
  slope_out <- if (length(right$idx) >= slope_window) mean(diff(mu[sort(right$idx)[1:slope_window]])) else NA_real_

  mu_c <- if (length(idx_center)) mu[idx_center] else mu[which.min(abs(x))]
  depth_z <- (fl_med - min(mu_c, na.rm = TRUE)) / (fl_mad + EPSILON)

  out[] <- c(auc_total, left$extent, right$extent, extent_tot,
             symmetry, slope_in, slope_out, depth_z)
  out
}


compute_halo <- function(wins, cov, params,
                         band_center_bp = 170, band_width_bp = 20,
                         annulus_inner = 80, annulus_outer = 500) {

  width_bp <- 2 * params$Wmax + 1L
  mu <- cov_profile_viewmeans(cov, wins$profile, width_bp)

  nm <- c("halo_power_170", "halo_qratio")
  if (!is.numeric(mu) || all(!is.finite(mu))) return(setNames(rep(NA_real_, 2L), nm))

  x <- as.integer(names(mu))

  idx_ann <- which(abs(x) >= annulus_inner & abs(x) <= annulus_outer)
  if (length(idx_ann) < 64) return(setNames(rep(NA_real_, 2L), nm))

  y <- mu[idx_ann] - mean(mu[idx_ann], na.rm = TRUE)
  n <- length(y)
  nfft <- 2^(ceiling(log2(n)))
  y_pad <- c(y, rep(0, nfft - n))
  Y <- stats::fft(y_pad)
  P <- (Mod(Y)^2) / nfft
  freq <- (0:(nfft-1)) / nfft

  f_lo <- 1 / (band_center_bp + band_width_bp)
  f_hi <- 1 / (band_center_bp - band_width_bp)

  main <- (freq >= f_lo) & (freq <= f_hi)
  side <- ((freq >= max(1/(band_center_bp*1.6), 0)) & (freq < f_lo)) |
          ((freq > f_hi) & (freq <= min(1/(band_center_bp*0.6), 0.5)))

  halo_power <- sum(P[main], na.rm = TRUE)
  side_power <- sum(P[side], na.rm = TRUE)

  qratio <- halo_power / (side_power + EPSILON)

  setNames(c(halo_power, qratio), nm)
}


compute_endmotif_features <- function(frag.gr, wins, params, k = 4) {

  center_frags <- frag.gr[overlapsAny(frag.gr, wins$center)]
  base_names <- c("A","C","G","T")

  out_names <- c(
    paste0("endmotif_", base_names, "_5p_freq"),
    paste0("endmotif_", base_names, "_3p_freq"),
    "endmotif_GC_5p","endmotif_GC_3p",
    "endmotif_entropy_5p","endmotif_entropy_3p"
  )

  out <- setNames(rep(NA_real_, length(out_names)), out_names)

  if (length(center_frags) < MIN_FRAGMENTS_FOR_STATS)
    return(out)

  tryCatch({

    # 5' kmers
    k5 <- get_5prime_kmer_intervals(center_frags, k)
    s5 <- getSeq(BSgenome.Hsapiens.UCSC.hg19, k5)

    # 3' kmers
    is_minus <- as.character(strand(center_frags)) == "-"
    s3_start <- ifelse(is_minus,
                       start(center_frags),
                       pmax(start(center_frags), end(center_frags) - k + 1L))
    s3_end   <- ifelse(is_minus,
                       pmax(start(center_frags), start(center_frags) + k - 1L),
                       end(center_frags))

    k3 <- GRanges(
      seqnames = seqnames(center_frags),
      ranges   = IRanges(start = s3_start, end = s3_end),
      strand   = strand(center_frags)
    )

    s3 <- getSeq(BSgenome.Hsapiens.UCSC.hg19, k3)

    w <- if ("weight" %in% names(mcols(center_frags)))
           mcols(center_frags)$weight else rep(1, length(center_frags))

    # 5' nucleotide frequencies
    f5 <- colSums(alphabetFrequency(s5, baseOnly = TRUE) * w)
    f5 <- f5 / sum(f5)

    # 3' nucleotide frequencies
    f3 <- alphabetFrequency(s3, baseOnly = TRUE)
    f3 <- colSums(f3)
    f3 <- f3 / sum(f3)

    # Populate base-level frequencies
    for (b in base_names) {
      if (b %in% names(f5))
        out[paste0("endmotif_", b, "_5p_freq")] <- f5[b]
      if (b %in% names(f3))
        out[paste0("endmotif_", b, "_3p_freq")] <- f3[b]
    }

    # GC content
    out["endmotif_GC_5p"] <- (f5["G"] + f5["C"]) / sum(f5[base_names])
    out["endmotif_GC_3p"] <- (f3["G"] + f3["C"]) / sum(f3[base_names])

    # Entropy (k-mer entropy)
    km5 <- oligonucleotideFrequency(s5, width = k)
    p5  <- prop.table(colSums(km5))
    ent5 <- -sum(p5 * log2(p5 + EPSILON))

    km3 <- oligonucleotideFrequency(s3, width = k)
    p3  <- prop.table(colSums(km3))
    ent3 <- -sum(p3 * log2(p3 + EPSILON))

    out["endmotif_entropy_5p"] <- ent5
    out["endmotif_entropy_3p"] <- ent3

  }, error = function(e) {
    cat("Warning: End motif calc failed:", e$message, "\n")
  })

  out
}


      
# ----------------------------------------------------------------------
# PROCESS ONE CRE FILE
# ----------------------------------------------------------------------
process_cre_file <- function(cre_file, cre_name, frag.gr, cov, params, enabled_modules) {
  cat("\nProcessing CRE:", cre_name, "\n")

  cre.raw <- safe_import_cre(cre_file)
  if (!length(cre.raw)) {
    cat("  ⚠️ CRE import empty:", cre_name, "\n")
    return(NULL)
  }

  # Align CRE seqlevels to matching fragment seqlevels
  cre.gr <- harmonize_to_frag(cre.raw, frag.gr)
  if (!length(cre.gr)) {
    cat("  ⚠️ CRE unusable after harmonization\n")
    return(NULL)
  }

  cat("  Regions:", length(cre.gr), "\n")
  summit_mode <- all(width(cre.gr) == 1L)

  if (summit_mode) {
    cat("  Mode: SUMMIT (all 1bp)\n")
    centers <- cre.gr
    center_win <- safe_center_window(centers, width_bp = 2*params$center + 1L)
    profile_win <- safe_center_window(centers, width_bp = 2*params$Wmax + 1L)
    wins <- list(
      center  = center_win,
      profile = profile_win,
      flankL  = safe_flank_window(centers, params$flank_left),
      flankR  = safe_flank_window(centers, params$flank_right)
    )
  } else {
    cat("  Mode: RANGE file\n")
    w <- width(cre.gr)
    mids <- start(cre.gr) + floor(w/2)
    mid.gr <- GRanges(seqnames = seqnames(cre.gr),
                      ranges   = IRanges(mids, width=1L),
                      strand   = strand(cre.gr))
    profile_win <- safe_center_window(mid.gr, width_bp = 2*params$Wmax+1L)
    wins <- list(
      center  = cre.gr,
      profile = profile_win,
      flankL  = safe_flank_window(mid.gr, params$flank_left),
      flankR  = safe_flank_window(mid.gr, params$flank_right)
    )
  }

  validate_windows(wins, params)
  overlap_cache <- cache_overlaps(frag.gr, wins)

  cat("  Center hits:", sum(overlap_cache$center),
      "| Flank hits:", sum(overlap_cache$flankL) + sum(overlap_cache$flankR), "\n")

  # Compute CRE features
  feats <- list()

  if ("relcov" %in% enabled_modules) {
    feats$relcov <- compute_rel_cov(frag.gr, wins, cov, params, overlap_cache)
  }
  if ("swfse" %in% enabled_modules) {
    feats$swfse <- compute_swfse(frag.gr, wins, params, overlap_cache)
  }
  if ("fld" %in% enabled_modules) {
    feats$fld <- compute_fld(frag.gr, wins, params, overlap_cache)
  }
  if ("qchasm" %in% enabled_modules) {
    feats$qchasm <- compute_qchasm(frag.gr, wins, cov, params)
  }
  if ("qchasm_outward" %in% enabled_modules) {
    feats$qchasm_outward <- compute_qchasm_outward(wins, cov, params)
  }
  if ("halo" %in% enabled_modules) {
    feats$halo <- compute_halo(wins, cov, params)
  }
  if ("endmotif" %in% enabled_modules) {
    feats$endmotif <- compute_endmotif_features(frag.gr, wins, params)
  }

  return(feats)
}

# ======================= FILE DISCOVERY ========================
discover_cre_files <- function(cre_dirs) {
  cre_dirs <- trimws(unlist(strsplit(cre_dirs, "[,;]")))
  cat("\n📂 CRE directories:\n")
  for (d in cre_dirs) cat("   -", d, "\n")
  allf <- unique(unlist(lapply(cre_dirs, function(d) {
    if (!dir.exists(d)) {
      cat("⚠️ dir missing:", d, "\n")
      return(character(0))
    }
    list.files(d,
               pattern    = "\\.(bed|rds)$",
               full.names = TRUE,
               recursive  = TRUE,
               ignore.case = TRUE)
  })))
  if (!length(allf)) stop("No CRE files.")
  base <- tolower(gsub("\\.(bed|rds)$", "", basename(allf)))
  uniq <- tapply(allf, base, function(f) {
    if (length(f) > 1) {
      r <- f[grepl("\\.rds$", f, ignore.case = TRUE)]
      if (length(r)) {
        cat("ℹ️ dup .rds+.bed; using", basename(r[1]), "and ignoring",
            paste(basename(setdiff(f, r[1])), collapse = ", "), "\n")
        r[1]
      } else {
        f[1]
      }
    } else {
      f[1]
    }
  })
  final <- unname(unlist(uniq))
  rel <- function(p) paste0(basename(dirname(p)), "/", basename(p))
  nms <- gsub("\\.(bed|rds)$", "", rel(final), ignore.case = TRUE)
  cat("📄 Found", length(final), "unique CRE files\n")
  list(files = final, names = nms)
}

# ============================ MAIN =============================
cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\ncfDNA CRE Fragmentomics Pipeline (v30.1) [mode=", mode, "]\n",
  paste(rep("=", 70), collapse = ""),
  "\n\n",
  sep = ""
)

# Discover fragment files
frag.paths <- sort(list.files(
  fragDir,
  pattern    = "\\.(rds|starch)$",
  full.names = TRUE,
  ignore.case = TRUE
))

if (length(frag.paths) < index || index < 1) {
  stop("Index out of range for fragDir files. Found ", length(frag.paths),
       ", requested ", index)
}

curr.path <- frag.paths[index]
s <- gsub("\\.starch$|\\.rds$", "", basename(curr.path))

dir.create(file.path(outDir, "tmp"), recursive = TRUE, showWarnings = FALSE)
outfile      <- file.path(outDir, "tmp", paste0(s, "_CRE_metrics.rds"))
outfile_flat <- gsub("_CRE_metrics.rds$", "_CRE_metrics_flat.rds", outfile)
meta_out     <- file.path(outDir, "tmp", paste0(s, "_CRE_metadata_", mode, ".yaml"))
log_out      <- file.path(outDir, "tmp", paste0(s, "_processing_", mode, ".log"))

# Skip if exists
if (file.exists(outfile)) {
  cat("✓ Sample already processed:", outfile, "\n")
  quit(save = "no")
}

sink(log_out, split = TRUE)
cat("Processing started:", as.character(Sys.time()), "\n")
cat("Sample:", s, "\n")
cat("Mode:", mode, "\n")
cat("Fragment file:", curr.path, "\n\n")

# ---------------------- Load fragments -------------------------
frag.gr <- tryCatch(
  load_fragments(
    fp             = curr.path,
    frag_min       = frag_min,
    frag_max       = frag_max,
    use_preprocessed = (mode == "preprocessed")
  ),
  error = function(e) {
    cat("ERROR loading fragments:", e$message, "\n")
    stop(e)
  }
)

# ---------------------- Mode-dependent handling ----------------
gc_validation_summary <- NULL

if (mode == "full") {
  # FULL: blacklist + GC correction (Cristiano)
  cat("\nMode 'full': applying blacklist + GC correction\n")
  target_file <- select_gc_reference(platform)

  frag.gr.original <- frag.gr

  frag.gr <- apply_fragment_gc_correction(
    frag.gr     = frag.gr,
    genome      = HG19_GENOME,
    blacklist   = filters.hg19,
    target_file = target_file,
    by_chrom    = TRUE
  )

  # ------------------ GC validation --------------------
  cat("\nValidating GC correction...\n")
  n_original  <- length(frag.gr.original)
  n_corrected <- length(frag.gr)
  retention_pct <- (n_corrected / n_original) * 100

  cat(sprintf("  Fragments: %s → %s (%.1f%% retained)\n",
              format(n_original, big.mark = ","),
              format(n_corrected, big.mark = ","),
              retention_pct))

  gc_issues   <- character(0)
  gc_warnings <- character(0)

  if (retention_pct < 30) {
    gc_issues <- c(gc_issues, sprintf("Low retention: %.1f%%", retention_pct))
  } else if (retention_pct < 50) {
    gc_warnings <- c(gc_warnings, sprintf("Moderate retention: %.1f%%", retention_pct))
  }

  if (n_corrected < 5000) {
    gc_issues <- c(gc_issues, sprintf("Too few fragments after GC: %s",
                                      format(n_corrected, big.mark = ",")))
  }

  if (!("weight" %in% names(mcols(frag.gr)))) {
    gc_issues <- c(gc_issues, "Weight column missing after GC correction")
  } else {
    weights <- mcols(frag.gr)$weight
    n_na       <- sum(is.na(weights))
    n_inf      <- sum(is.infinite(weights))
    n_zero     <- sum(weights == 0, na.rm = TRUE)
    n_negative <- sum(weights < 0,  na.rm = TRUE)

    weight_mean  <- mean(weights, na.rm = TRUE)
    weight_range <- range(weights, na.rm = TRUE)

    cat(sprintf("  Weights: mean=%.3f, range=[%.3f, %.3f]\n",
                weight_mean, weight_range[1], weight_range[2]))

    if (n_na > 0)       gc_issues   <- c(gc_issues,   sprintf("%d NA weights", n_na))
    if (n_inf > 0)      gc_issues   <- c(gc_issues,   sprintf("%d infinite weights", n_inf))
    if (n_zero > 0)     gc_warnings <- c(gc_warnings, sprintf("%d zero weights", n_zero))
    if (n_negative > 0) gc_issues   <- c(gc_issues,   sprintf("%d negative weights", n_negative))

    if (weight_range[1] > 0 && weight_range[2] > 0) {
      wr <- weight_range[2] / weight_range[1]
      if (wr > 100) {
        gc_warnings <- c(gc_warnings, sprintf("Large weight ratio: %.1f", wr))
      }
    }
  }

  # Check GC bin coverage
  if ("gc_bin" %in% names(mcols(frag.gr))) {
    n_gc_bins <- length(unique(mcols(frag.gr)$gc_bin))
    cat(sprintf("  GC bins covered: %d\n", n_gc_bins))
    if (n_gc_bins < 20) {
      gc_warnings <- c(gc_warnings, sprintf("Few GC bins: %d", n_gc_bins))
    }
  }

  # Sanity check coverage on major chromosome
  if (length(frag.gr) > 0 && "weight" %in% names(mcols(frag.gr))) {
    chr_counts <- table(as.character(seqnames(frag.gr)))
    major_chr <- names(chr_counts)[which.max(chr_counts)]
    idx_chr <- which(seqnames(frag.gr) == major_chr)
    if (length(idx_chr) > 0) {
      idx_sub   <- idx_chr[seq_len(min(1000L, length(idx_chr)))]
      frag_test <- frag.gr[idx_sub]
      w_test    <- mcols(frag.gr)$weight[idx_sub]
      cov_test  <- coverage(frag_test, weight = w_test)
      cov_vals  <- sapply(cov_test, function(x) sum(x, na.rm = TRUE))
      if (all(cov_vals == 0)) {
        gc_issues <- c(gc_issues,
                       sprintf("Zero weighted coverage in test subset on %s", major_chr))
      }
    }
  }

  if (length(gc_issues) > 0) {
    cat("\n❌ GC CORRECTION VALIDATION FAILED:\n")
    for (ii in gc_issues) cat("  ✗", ii, "\n")
    stop("GC correction failed validation.")
  }

  if (length(gc_warnings) > 0) {
    cat("\n⚠️  GC CORRECTION WARNINGS:\n")
    for (w in gc_warnings) cat("  !", w, "\n")
  } else {
    cat("  ✓ GC validation passed with no warnings\n")
  }

  gc_validation_summary <- list(
    enabled       = TRUE,
    passed        = TRUE,
    retention_pct = retention_pct,
    n_original    = n_original,
    n_corrected   = n_corrected,
    n_warnings    = length(gc_warnings),
    warnings      = if (length(gc_warnings)) gc_warnings else NULL
  )

  rm(frag.gr.original)

} else if (mode == "nogc") {
  # NOGC: blacklist only, no GC weighting, weights = 1
  cat("\nMode 'nogc': applying blacklist only, no GC correction (weights=1)\n")
  frag.gr <- apply_blacklist_only(frag.gr, filters.hg19)
  mcols(frag.gr)$weight <- rep(1, length(frag.gr))
  gc_validation_summary <- list(
    enabled = FALSE,
    passed  = TRUE,
    note    = "mode='nogc': blacklist only, no GC; weights=1"
  )

} else if (mode == "preprocessed") {
  cat("\nMode 'preprocessed': skipping blacklist + GC, using precomputed weights\n")

  gc_issues   <- character(0)
  gc_warnings <- character(0)

  if (!("weight" %in% names(mcols(frag.gr)))) {
      gc_issues <- c(gc_issues, "Weight column missing in preprocessed fragments")
  } else {
      weights <- mcols(frag.gr)$weight
      n_na       <- sum(is.na(weights))
      n_inf      <- sum(is.infinite(weights))
      n_zero     <- sum(weights == 0, na.rm = TRUE)
      n_negative <- sum(weights < 0,  na.rm = TRUE)
    
      weight_mean  <- mean(weights, na.rm = TRUE)
      weight_range <- range(weights, na.rm = TRUE)
    
      cat(sprintf("  Weights: mean=%.3f, range=[%.3f, %.3f]\n",
                  weight_mean, weight_range[1], weight_range[2]))
    
      if (n_na > 0)       gc_issues   <- c(gc_issues,   sprintf("%d NA weights", n_na))
      if (n_inf > 0)      gc_issues   <- c(gc_issues,   sprintf("%d infinite weights", n_inf))
      if (n_zero > 0)     gc_warnings <- c(gc_warnings, sprintf("%d zero weights", n_zero))
      if (n_negative > 0) gc_issues   <- c(gc_issues,   sprintf("%d negative weights", n_negative))
    
      if (weight_range[1] > 0 && weight_range[2] > 0) {
        wr <- weight_range[2] / weight_range[1]
        if (wr > 100) {
          gc_warnings <- c(gc_warnings, sprintf("Large weight ratio: %.1f", wr))
        }
      }
    }
}

# ------------------ Coverage (always uses current weights) ----------------
cat("\nCalculating genome coverage (weighted)...\n")
if (!("weight" %in% names(mcols(frag.gr)))) {
  stop("Internal error: 'weight' column missing before coverage() call.")
}
wt  <- mcols(frag.gr)$weight
cov <- coverage(frag.gr, weight = wt)


fl_left  <- parse_window(flank_left)
fl_right <- parse_window(flank_right)

# Ensure coverage profile window spans center + flanks
Wmax <- max(
  c(
    abs(fl_left),
    abs(fl_right),
    abs(center_bp)
  ),
  na.rm = TRUE
)

params <- list(
  center      = center_bp,
  flank_left  = fl_left,
  flank_right = fl_right,
  Wmax        = Wmax,
  frag_min    = frag_min,
  frag_max    = frag_max,
  frag_window = 20,
  frag_step   = 10
)

cat("\nParameters:\n  Center window: ±", center_bp, " bp\n", sep = "")
cat("  Left flank: ", fl_left[1],  " to ", fl_left[2],  "\n", sep = "")
 cat("  Right flank:", fl_right[1], " to ", fl_right[2], "\n", sep = "")
cat("  Frag range: ", frag_min, "-", frag_max, " bp\n", sep = "")
cat("  Mode:", mode, " (platform:", platform, ")\n")

# ------------------ CRE discovery --------------------------
cre_info <- discover_cre_files(creDirs)

feature_results <- setNames(vector("list", length(enabled_modules)), sort(enabled_modules))
processing_log  <- list()

cat("\n", paste(rep("-", 70), collapse=""), "\n",
    "Beginning feature extraction\n",
    paste(rep("-", 70), collapse=""), "\n", sep="")

for (i in seq_along(cre_info$files)) {
  cf <- cre_info$files[i]
  cname <- cre_info$names[i]

  result <- tryCatch({
    feats <- process_cre_file(
      cre_file       = cf,
      cre_name       = cname,
      frag.gr        = frag.gr,
      cov            = cov,
      params         = params,
      enabled_modules= enabled_modules
    )
    processing_log[[cname]] <- list(
      status     = "success",
      n_features = if (length(feats)) length(unlist(feats)) else 0
    )
    feats
  }, error = function(e) {
    cat("  ❌ ERROR in", cname, ":", e$message, "\n")
    processing_log[[cname]] <- list(
      status  = "error",
      message = e$message
    )
    NULL
  })

  if (!is.null(result)) {
    for (fam in names(result)) {
      if (is.null(feature_results[[fam]])) feature_results[[fam]] <- list()
      feature_results[[fam]][[cname]] <- result[[fam]]
    }
  }
}

# ------------------ Combine features ------------------------
cat("\n", paste(rep("-", 70), collapse=""), "\n",
    "Combining results\n",
    paste(rep("-", 70), collapse=""), "\n\n", sep="")

features_out   <- list()
feature_schema <- list()

for (fam in names(feature_results)) {
  fam_list <- feature_results[[fam]]
  if (!is.null(fam_list) && length(fam_list) > 0) {
    cat("Feature family:", fam, "- combining", length(fam_list), "CRE sets\n")
    mat <- tryCatch({
      do.call(rbind, fam_list)
    }, error = function(e) {
      # force numeric & re-bind if types differ
      lst <- lapply(fam_list, function(x) {
        y <- as.numeric(x)
        names(y) <- names(x)
        y
      })
      do.call(rbind, lst)
    })
    if (is.null(dim(mat))) {
      mat <- t(as.matrix(mat))
    }
    rownames(mat) <- names(fam_list)
    features_out[[fam]]   <- mat
    feature_schema[[fam]] <- colnames(mat)
    cat("  →", nrow(mat), "CRE sets ×", ncol(mat), "features\n")
  }
}

# Optionally make a flat table of all features
flat <- NULL
if (length(features_out)) {
  fam_names <- names(features_out)
  prefixed <- lapply(fam_names, function(fn) {
    m <- features_out[[fn]]
    colnames(m) <- paste0(fn, "::", colnames(m))
    as.data.frame(m, check.names = FALSE)
  })
  flat <- Reduce(function(a, b) {
    rn <- intersect(rownames(a), rownames(b))
    a2 <- a[rn, , drop = FALSE]
    b2 <- b[rn, , drop = FALSE]
    cbind(a2, b2)
  }, prefixed)
}

# ------------------ Build final result object ---------------
final_results <- list(
  sample_name    = s,
  features       = features_out,
  cre_names      = if (length(features_out)) rownames(features_out[[names(features_out)[1]]]) else character(0),
  feature_schema = feature_schema,
  processing_log = processing_log,
  parameters     = list(
    analysis_type   = "CRE",
    cre_directories = unlist(strsplit(creDirs, "[,;]")),
    windows         = list(
      center      = center_bp,
      flank_left  = fl_left,
      flank_right = fl_right,
      Wmax        = Wmax
    ),
    fragments       = list(
      min = frag_min,
      max = frag_max
    ),
    enabled_modules = enabled_modules,
    mode            = mode,
    gc_correction   = (mode == "full"),
    platform        = platform,
    gc_validation   = gc_validation_summary,
    processing_info = list(
      timestamp        = as.character(Sys.time()),
      R_version        = as.character(R.version$version.string),
      total_fragments  = length(frag.gr),
      total_cre_files  = length(cre_info$files),
      pipeline_version = "WGS-CRE-v30.1-mode"
    ),
    notes = list(
      mode_logic    = "Per-file: SUMMIT (1bp) vs RANGE (>1bp); profiles from wins$profile; overlaps from wins$center",
      dedup_inputs  = "RDS preferred over BED when basenames match",
      gc_modes      = "mode='full': blacklist+GC; 'nogc': blacklist only, weights=1; 'preprocessed': uses pre-weighted fragments, no GC/blacklist",
      optimizations = "Cached overlaps, weighted coverage"
    )
  )
)

# ------------------ Save outputs ----------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n",
    "Saving results\n",
    paste(rep("=", 70), collapse=""), "\n\n", sep="")

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(final_results, outfile)
cat("✓ Features saved:", outfile, "\n")

if (!is.null(flat)) {
  saveRDS(flat, outfile_flat)
  cat("✓ Flat table saved:", outfile_flat, "\n")
}

# ------------------ Metadata YAML ---------------------------
metadata <- list(
  Sample            = s,
  Timestamp         = as.character(Sys.time()),
  CRE_directories   = unlist(strsplit(creDirs, "[,;]")),
  Total_CRE_files   = length(cre_info$files),
  Enabled_modules   = enabled_modules,
  Feature_families  = lapply(features_out, ncol),
  Total_features    = sum(sapply(features_out, ncol)),
  Fragment_count    = length(frag.gr),
  Mode              = mode,
  GC_correction     = (mode == "full"),
  Platform          = platform,
  GC_validation     = gc_validation_summary,
  Processing_status = lapply(processing_log, function(x) x$status),
  Pipeline_version  = "WGS-CRE-v30.1-mode"
)

write_yaml(metadata, meta_out)
cat("✓ Metadata saved:", meta_out, "\n")

# ------------------ Finish log & summary ---------------------
sink()
cat("✓ Log saved:", log_out, "\n")

cat("\n", paste(rep("=", 70), collapse=""), "\n",
    "ANALYSIS COMPLETE\n",
    paste(rep("=", 70), collapse=""), "\n", sep="")
cat("Sample:", s, "\n")
cat("Total fragments analyzed:", length(frag.gr), "\n")
cat("Total CRE sets processed:", length(cre_info$files), "\n")
cat("Successful:", sum(sapply(processing_log, function(x) x$status == "success")), "\n")
cat("Failed:",     sum(sapply(processing_log, function(x) x$status == "error")), "\n")
cat("Total feature families:", length(features_out), "\n")
cat("Total features extracted:", sum(sapply(features_out, ncol)), "\n")
cat("\nOutput files:\n  - Results:", outfile,
    "\n  - Flat:", outfile_flat,
    "\n  - Metadata:", meta_out,
    "\n  - Log:", log_out, "\n\n")
