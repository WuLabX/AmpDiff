# AmpDiff

AmpDiff is a motif-guided diffusion framework for antimicrobial peptide (AMP)
sequence optimization. It supports AMP pretraining, activity-guided finetuning,
motif-constrained sampling, and downstream evaluation with AMP activity and
hemolysis predictors.

The project is organized around two generation modes:

- `de`: fixes active AMP motifs and optimizes the remaining positions.
- `inp`: fixes both active motifs and hemolysis-related motifs, with optional
  HemoPI2 hemolysis guidance.

## Repository Layout

- `amp_scripts/`: training, finetuning, sampling, evaluation, and experiment
  scripts for AmpDiff.
- `configs/`: YAML configuration files for AMP pretraining and finetuning.
- `dataset/`: dataset preparation utilities and AMP dataset loaders.
- `data/`: prepared AmpDiff pretrain, finetune, test, and motif files.
- `checkpoints/`: canonical pretrained and finetuned model checkpoints.
- `evaluation/`: generated variants and downstream analysis outputs.
- `doc/`: project notes, dataset plans, and guidance documentation.

## Environment

Create the conda environment from the provided file:

```bash
conda env create -f environment.yaml
```

Activate the created environment before running the commands below. If your
local scorer paths differ from the defaults, pass the corresponding
`--pepnet_root`, `--hemopi2_root`, `--amppred_root`, `--iamp_root`, or
`--unidl_root` arguments.

## Data

Prepared AmpDiff datasets are included under `data/`:

- `data/pretrain/ampdiff_pretrain_all.csv`
- `data/pretrain/ampdiff_pretrain_motif.csv`
- `data/finetune/ampdiff_finetune_de.csv`
- `data/finetune/ampdiff_finetune_inp.csv`
- `data/test/ampdiff_test_de.csv`
- `data/test/ampdiff_test_inp.csv`

FASTA versions are provided for the pretrain and test splits where needed. The
dataset split summary is available in `data/dataset_report.md`.

To rebuild the prepared datasets from the merged raw AMP files:

```bash
python dataset/prepare_ampdiff_datasets.py
```

The preparation script reads from `raw dataset/` by default and writes the
processed outputs to `data/`.

## Pretraining

AmpDiff uses a two-stage pretraining workflow.

Stage 1 trains on the broad AMP pretraining set:

```bash
bash amp_scripts/run_pretrain_stage1_all.sh
```

Stage 2 initializes from the stage 1 checkpoint and trains on the active-motif
subset:

```bash
bash amp_scripts/run_pretrain_stage2_motif.sh
```

The canonical checkpoints are written to:

- `checkpoints/pretrain/stage1_all.pt`
- `checkpoints/pretrain/stage2_motif.pt`

You can also run `amp_scripts/amp_train.py` directly:

```bash
python amp_scripts/amp_train.py \
  --data_path data/pretrain/ampdiff_pretrain_motif.csv \
  --config_path configs/amp_pretrain_stage2_motif.yml \
  --log_path logs/amp_pretrain \
  --train_loss fr \
  --motif_type prosite,regular,merci \
  --mode de \
  --scorer none \
  --require_active_motif True \
  --min_optimizable 1 \
  --init_checkpoint checkpoints/pretrain/stage1_all.pt \
  --ckpt_subdir pretrain \
  --ckpt_name stage2_motif.pt
```

## Finetuning

Finetuning uses supervised AMP activity data and one selectable activity
guidance predictor. Supported activity guidance predictors are:

- `pepnet`
- `amppred_mfa`
- `iamp_attenpred`
- `unidl4biopep`

Finetune the `de` model:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/ampdiff_finetune_de.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/stage2_motif.pt \
  --guidance_scorer pepnet \
  --mode de
```

Finetune the `inp` model with optional hemolysis guidance:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/ampdiff_finetune_inp.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/stage2_motif.pt \
  --guidance_scorer pepnet \
  --mode inp \
  --hemolysis-guidance
```

Canonical finetuned checkpoints are written to:

- `checkpoints/de/de.pt`
- `checkpoints/inp/inp.pt`

## Sampling

Use `amp_scripts/amp_sample.py` to optimize seed AMP sequences from a FASTA
file.

Example for `de` mode:

```bash
python amp_scripts/amp_sample.py \
  --ckpt checkpoints/de/de.pt \
  --input_fasta data/test/ampdiff_test_de.fasta \
  --mode de \
  --guidance_scorer pepnet \
  --independent_scorers auto \
  --num_samples 10 \
  --batch_size 10 \
  --temperature 1.0 \
  --output outputs/de_sampling
```

Example for `inp` mode with hemolysis guidance and hemolysis evaluation:

```bash
python amp_scripts/amp_sample.py \
  --ckpt checkpoints/inp/inp.pt \
  --input_fasta data/test/ampdiff_test_inp.fasta \
  --mode inp \
  --guidance_scorer pepnet \
  --hemolysis-guidance \
  --eval-hemolysis \
  --independent_scorers auto \
  --num_samples 10 \
  --batch_size 10 \
  --temperature 1.0 \
  --output outputs/inp_sampling
```

Sampling writes `sample_amp_result.csv`, `generated_amps.fasta`, and, when
post-generation evaluators are enabled, `scores.tsv`.

For the main generation experiments across checkpoints, modes, temperatures,
and seeds:

```bash
bash amp_scripts/run_main_generation.sh
```

## Evaluation

Evaluate generated AMP FASTA files with any registered scorers:

```bash
python amp_scripts/amp_eval.py \
  --fasta outputs/de_sampling/Seed2024_Modede_Temp1.0/generated_amps.fasta \
  --scorer pepnet,amppred_mfa,iamp_attenpred,unidl4biopep \
  --output_dir outputs/de_sampling_eval
```

Add `hemopi2` to the scorer list when hemolysis scoring is required:

```bash
python amp_scripts/amp_eval.py \
  --fasta outputs/inp_sampling/Seed2024_Modeinp_Temp1.0/generated_amps.fasta \
  --scorer pepnet,amppred_mfa,iamp_attenpred,unidl4biopep,hemopi2 \
  --output_dir outputs/inp_sampling_eval
```

Additional property and substitution analysis is available through:

```bash
python amp_scripts/amp_property_analysis.py
```

## Motifs and Guidance

Active and hemolysis-related motif files are stored under `data/motif/`:

- `data/motif/AMP_act_Motifs.csv`
- `data/motif/AMP_hom_negMotifs.csv`

Motif sources can be selected with:

```bash
--motif_type prosite,regular,merci
```

Use `--motif_type none` to disable motif matching where supported. See
`doc/activity_guidance.md` for details about guidance scorer selection and
independent post-generation evaluation.

## Citation

If you use AmpDiff in your research, please cite this repository or the
corresponding AmpDiff manuscript when available.
