#!/bin/bash
#SBATCH --job-name=nmd_infer_all
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/infer_all_%j.out
#SBATCH --error=logs/infer_all_%j.err

# Run inference on every isoform in the model's H5 (split="all") and write
# results to results_4ct/predictions_all_atg500_stop500.tsv. Then scp that
# file to your local ~/claude_projects/nmd/code/nmd_predictor_comparison/.

set -e
mkdir -p logs

source ~/miniconda3/etc/profile.d/conda.sh   # or wherever conda init lives
conda activate nmd_model

cd ~/cc/nmd_orf_model_v5_4ct

python run_infer_all.py \
    --config config.yaml \
    --atg-window 500 \
    --stop-window 500 \
    --split all

echo "DONE. Output at: ~/cc/nmd_orf_model_v5_4ct/results_4ct/predictions_all_atg500_stop500.tsv"
