#!/bin/bash
#SBATCH --job-name=tf_access
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
fragDir="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/results/corrected-LUCAS"

tfbsDir="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/ref/TFBS/indv"

outDir="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/results/access-LUCAS"

script="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/scripts/utils_accessibillity-precompute-module-v2.R"

mkdir -p "$outDir"
mkdir -p logs

# ===================================================================
# Parameters
# ===================================================================
window=3000
bin_size=25
frag_min=1
frag_max=2000
platform="hiseq"

# ✔ NEW PARAMETER: save_mode = combined | per_tf
save_mode="per_tf"

# ===================================================================
# Status Info
# ===================================================================
echo "[$(date)] Starting TF Accessibility Computation"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Fragment directory: $fragDir"
echo "TFBS directory: $tfbsDir"
echo "Output directory: $outDir"
echo "Save mode: $save_mode"

# ===================================================================
# Run the updated R script
# ===================================================================
time Rscript "$script" \
  "$SLURM_ARRAY_TASK_ID" \
  "$fragDir" \
  "$tfbsDir" \
  "$outDir" \
  "$window" \
  "$bin_size" \
  "$frag_min" \
  "$frag_max" \
  "$platform" \
  "$save_mode"

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID completed."
