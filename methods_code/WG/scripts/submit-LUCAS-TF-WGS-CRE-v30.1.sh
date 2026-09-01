#!/bin/bash
#SBATCH --job-name=wgs_cre_v30.1
#SBATCH --array=1-386
#SBATCH --cpus-per-task=1
#SBATCH --mem=24GB
#SBATCH --time=7-0:00:00
#SBATCH --output=logs/slurm-%x-%J-%a.out.txt
#SBATCH --error=logs/slurm-%x-%J-%a.err.txt


# ===================================================================
# Environment setup
# ===================================================================
module load conda
conda activate reproduce-lucas-wflow
export OMP_NUM_THREADS=1

# ===================================================================
# Input / Output directories
# ===================================================================
# RAW OR PREPROCESSED FRAGMENT DIRECTORY
fragDir=/dcs04/scharpf/data/annapragada/Data/granges/LUCAS_Hi

# Multiple CRE directories separated by commas or semicolons
creDirs="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/ref/TFBS/indv"

# Output directory
outDir=/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/results/corrected-LUCAS

# Path to MONO V30 script
script=/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/scripts/WGS-CRE-v30.1-Mono.R

mkdir -p "$outDir"
mkdir -p logs

# ===================================================================
# Parameters
# ===================================================================
center_bp=100
flank_left="-3000:-2750"
flank_right="2750:3000"
frag_min=100
frag_max=220

# mode = full | nogc | preprocessed
#   full         → blacklist + GC correction
#   nogc         → blacklist only, weights=1
#   preprocessed → skip blacklist+GC, use preweighted files
mode="preprocessed"  # <-- CHANGE THIS AS NEEDED

platform="hiseq"

# Only compute selected feature families
enabled_modules="relcov,swfse,fld"

# ===================================================================
# Run WGS-CRE Pipeline
# ===================================================================
echo "[$(date)] Starting WGS-CRE v30.1 Mono (mode=$mode)"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Fragment directory: $fragDir"
echo "CRE directories: $creDirs"
echo "Output directory: $outDir"

time Rscript "$script" \
  "$SLURM_ARRAY_TASK_ID" \
  "$fragDir" \
  "$creDirs" \
  "$outDir" \
  "$center_bp" \
  "$flank_left" \
  "$flank_right" \
  "$frag_min" \
  "$frag_max" \
  "$mode" \
  "$platform" \
  "$enabled_modules"

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID completed."
