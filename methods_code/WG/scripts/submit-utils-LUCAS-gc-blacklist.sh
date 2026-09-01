#!/bin/bash
#SBATCH                   #$ -cwd
#SBATCH                   #$  -j n
#SBATCH                   #$ -l h_fsize=100G
#SBATCH --mem=24GB        #$ -l mem_free=1G
#SBATCH                   #$ -l h_vmem=5G
#SBATCH --time=7-00:00:00 #$ -l h_rt=730:00:00
#SBATCH --array=1-386
#SBATCH --job-name=utils_gc-blacklist-module
#SBATCH --cpus-per-task=1
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
fragDir=/dcs04/scharpf/data/annapragada/Data/granges/LUCAS_Hi
outDir=/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/results/corrected-LUCAS

script=/dcs07/scharpf/data/jzavras/LUCAS_Olink/Lucas-Cancer-Screening/scripts/lung-cancer-screening-paper/methods_code/WG/scripts/utils_gc-blacklist-module.R

mkdir -p "$outDir"
mkdir -p logs

# ===================================================================
# Parameters
# ===================================================================
platform="hiseq"
frag_min=0
frag_max=0

# ===================================================================
# Run GC Correction, QC, and Blacklist Filtering
# ===================================================================
echo "[$(date)] Starting Utils GC and Blacklist Corrections"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

time Rscript "$script" \
  "$SLURM_ARRAY_TASK_ID" \
  "$fragDir" \
  "$outDir" \
  "$platform" \
  "$frag_min" \
  "$frag_max"

echo "[$(date)] Task $SLURM_ARRAY_TASK_ID completed."
