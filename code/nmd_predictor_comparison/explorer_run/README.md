# Full-cohort inference on Explorer

We need the trained NMD predictor's score on every isoform in the model's H5 (`split="all"`), not just the held-out test set. The current `evaluate.py` writes only the test-set predictions; this runner replicates the same model-loading + inference loop but with `NMDDataset(split="all")`, so every isoform in `nmd_orf_data.h5` gets a row in the output.

## Run

```bash
# On Explorer
ssh p.castaldi@explorer.northeastern.edu
cd ~/cc/nmd_orf_model_v5_4ct

# Pull the runner into the model repo (this script depends on model.py/utils.py)
cp ~/claude_projects/nmd/code/nmd_predictor_comparison/explorer_run/run_infer_all.py .

# Submit (≤30 min on a single GPU)
sbatch ~/claude_projects/nmd/code/nmd_predictor_comparison/explorer_run/slurm_infer_all.sh

# OR run interactively in a GPU session
srun --partition=gpu --gres=gpu:1 --time=01:00:00 --pty bash
conda activate nmd_model
python run_infer_all.py --config config.yaml --atg-window 500 --stop-window 500
```

The job writes:

```
~/cc/nmd_orf_model_v5_4ct/results_4ct/predictions_all_atg500_stop500.tsv
```

with columns `isoform_id, chr, h5_split, label, logit, prob`.

## Pull predictions back to the local repo

```bash
# On laptop (any directory)
scp p.castaldi@explorer.northeastern.edu:~/cc/nmd_orf_model_v5_4ct/results_4ct/predictions_all_atg500_stop500.tsv \
    ~/claude_projects/nmd/code/nmd_predictor_comparison/
```

## Then re-run the comparison

Edit `01_extract_our_isoforms.R` to source predictions from
`predictions_all_atg500_stop500.tsv` (instead of the test-only
`predictions_atg500_stop500.tsv`) and re-run 01 → 04 + the figure.

The expected improvement: head-to-head intersection grows from 561 to
the full size of the model's H5 overlap with our 2,332 cohort (likely
close to the full 2,332 once the train + val splits are included).
