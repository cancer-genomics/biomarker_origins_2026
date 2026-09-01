#!/bin/bash
#SBATCH --job-name=split_tf_cache
#SBATCH --partition=cancergen
#SBATCH --array=1-384
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=3-0:00:00
#SBATCH --output=logs/slurm-%x-%J-%a.out.txt
#SBATCH --error=logs/slurm-%x-%J-%a.err.txt

# ===================================================================
module load conda
conda activate reproduce-lucas-wflow
export OMP_NUM_THREADS=1

# ===================================================================
inputDir="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/results/access-LUCAS"

outputDir="$inputDir"   # per-sample folders created inside same directory

script="/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/scripts/utils-split_cached_tf_matrices.R"

mkdir -p logs

echo "[$(date)] Starting TF cache split job"
echo "Task index: $SLURM_ARRAY_TASK_ID"

time Rscript "$script" \
  "$SLURM_ARRAY_TASK_ID" \
  "$inputDir" \
  "$outputDir"

echo "[$(date)] Completed index $SLURM_ARRAY_TASK_ID."
