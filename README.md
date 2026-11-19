# AmpDiff

AmpDiff is a motif-guided diffusion framework for antimicrobial peptide (AMP)
sequence optimization. It supports AMP pretraining, activity-guided finetuning,
motif-constrained sampling, and downstream evaluation with AMP activity and
hemolysis predictors.

The project is organized around two generation modes:

- `act`: fixes active AMP motifs and optimizes the remaining positions.
- `bal`: fixes both active motifs and hemolysis-related motifs, with optional
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
- `data/finetune/ampdiff_finetune_act.csv`
- `data/finetune/ampdiff_finetune_bal.csv`
- `data/test/ampdiff_test_act.csv`
- `data/test/ampdiff_test_bal.csv`

FASTA versions are provided for the pretrain and test splits where needed. The
dataset split summary is available in `data/dataset_report.md`.

To rebuild the prepared datasets from the merged raw AMP files:

```bash
python dataset/prepare_ampdiff_datasets.py
```

The preparation script reads from `raw dataset/` by default and writes the
processed outputs to `data/`.

## Pretraining

AmpDiff uses a two-stage pretraining (`stage 1`/`stage 2`) workflow and provides two independent usage paths (`A`/`B`).

### Path A: Raw-data retraining

Use this path when you want to rebuild the processed AMP datasets and train
both pretraining stages from raw AMP records.

1. Download the [raw data](https://doi.org/10.5281/zenodo.24358304) and place the files under `raw dataset/`.
2. Build the processed datasets:

```bash
python dataset/prepare_ampdiff_datasets.py \
  --raw-dir "raw dataset" \
  --out-dir data \
  --work-dir dataset/work \
  --motif-dir data/motif
```

3. Run `stage 1` on the broad AMP pretraining set:

```bash
bash amp_scripts/run_pretrain_stage1_all.sh
```

4. Run `stage 2`, which initializes from the stage 1 checkpoint and trains on
   the active-motif subset:

```bash
bash amp_scripts/run_pretrain_stage2_motif.sh
```

The two stages write their canonical checkpoints to:

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
  --mode act \
  --scorer none \
  --require_active_motif True \
  --min_optimizable 1 \
  --init_checkpoint checkpoints/pretrain/stage1_all.pt \
  --ckpt_subdir pretrain \
  --ckpt_name stage2_motif.pt
```

### Path B: Warm-start from checkpoints

Use this path when you want to skip pretraining and start from released
weights.

1. Download the [pretraining checkpoints](https://doi.org/10.5281/zenodo.24358304).
2. Place the checkpoint files at:

```text
checkpoints/pretrain/stage1_all.pt
checkpoints/pretrain/stage2_motif.pt
```

After completing either Path A or Path B, continue to finetuning with
`checkpoints/pretrain/stage2_motif.pt`.

## Finetuning

Finetuning uses supervised AMP activity data and one selectable activity
guidance predictor. Supported activity guidance predictors are:

- [PepNet](https://github.com/Harkool/PepNet)
- [AMPpred-MFA](https://github.com/Jiangle525/AMPpred-MFA)
- [iAMP-Attenpred](https://github.com/xingwxzz/iAMP-Attenpred)
- [UniDL4BioPep](https://github.com/dzjxzyd/UniDL4BioPep)

Optional hemolysis guidance is provided by:

- [HemoPI2](https://github.com/raghavagps/hemopi2)

### Path A: Raw-data retraining

Use this path when you want to build the supervised finetuning datasets and
train the `act` and `bal` models. Complete either pretraining path first so
that `checkpoints/pretrain/stage2_motif.pt` is available.

1. Download the [raw data](https://doi.org/10.5281/zenodo.24358304) and place the files under `raw dataset/`.
2. Build the processed datasets if they have not already been prepared:

```bash
python dataset/prepare_ampdiff_datasets.py \
  --raw-dir "raw dataset" \
  --out-dir data \
  --work-dir dataset/work \
  --motif-dir data/motif
```

3. Finetune the `act` model:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/ampdiff_finetune_act.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/stage2_motif.pt \
  --guidance_scorer pepnet \
  --mode act
```

4. Finetune the `bal` model with optional hemolysis guidance:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/ampdiff_finetune_bal.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/stage2_motif.pt \
  --guidance_scorer pepnet \
  --mode bal \
  --hemolysis-guidance
```

The canonical finetuned checkpoints are written to:

- `checkpoints/act/act.pt`
- `checkpoints/bal/bal.pt`

### Path B: Warm-start from checkpoints

Use this path when you want to skip finetuning and directly run sampling or
evaluation from released weights.

1. Download the [finetuned checkpoints](https://doi.org/10.5281/zenodo.24358304).
2. Place the checkpoint files at:

```text
checkpoints/act/act.pt
checkpoints/bal/bal.pt
```

After completing either Path A or Path B, continue to the sampling workflow
below.

## Sampling

Use `amp_scripts/amp_sample.py` to optimize seed AMP sequences from a FASTA
file.

Example for `act` mode:

```bash
python amp_scripts/amp_sample.py \
  --ckpt checkpoints/act/act.pt \
  --input_fasta data/test/ampdiff_test_act.fasta \
  --mode act \
  --guidance_scorer pepnet \
  --independent_scorers auto \
  --num_samples 10 \
  --batch_size 10 \
  --temperature 1.0 \
  --output outputs/act_sampling
```

Example for `bal` mode with hemolysis guidance and hemolysis evaluation:

```bash
python amp_scripts/amp_sample.py \
  --ckpt checkpoints/bal/bal.pt \
  --input_fasta data/test/ampdiff_test_bal.fasta \
  --mode bal \
  --guidance_scorer pepnet \
  --hemolysis-guidance \
  --eval-hemolysis \
  --independent_scorers auto \
  --num_samples 10 \
  --batch_size 10 \
  --temperature 1.0 \
  --output outputs/bal_sampling
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
  --fasta outputs/act_sampling/Seed2024_Modeact_Temp1.0/generated_amps.fasta \
  --scorer pepnet,amppred_mfa,iamp_attenpred,unidl4biopep \
  --output_dir outputs/act_sampling_eval
```

Add `hemopi2` to the scorer list when hemolysis scoring is required:

```bash
python amp_scripts/amp_eval.py \
  --fasta outputs/bal_sampling/Seed2024_Modebal_Temp1.0/generated_amps.fasta \
  --scorer pepnet,amppred_mfa,iamp_attenpred,unidl4biopep,hemopi2 \
  --output_dir outputs/bal_sampling_eval
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

```bibtex
@article{HAmpDiff,
  title = {Function-aware autoregressive diffusion enables safe and efficient antimicrobial peptide optimization},
  author = {Wu, Chuya},
  journal = {Brief. Bioinform.},
  year = {2026}
}
```

## Contact

If you have any questions or suggestions regarding this work, please feel free to contact us:

- Zhenyu Yue: [zhenyuyue@ahau.edu.cn](mailto:zhenyuyue@ahau.edu.cn)
- Chuya Wu: [wuchuya@stu.ahau.edu.cn](mailto:wuchuya@stu.ahau.edu.cn)
