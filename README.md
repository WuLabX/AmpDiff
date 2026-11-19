# AmpDiff

AmpDiff: Function-aware autoregressive diffusion enables safe and efficient
antimicrobial peptide optimization. It supports AMP pretraining,
activity-guided finetuning, motif-constrained sampling, and downstream
evaluation with AMP activity and hemolysis predictors.

The project is organized around two generation modes:

- `act`: fixes active AMP motifs and optimizes the remaining positions, with
  activity guidance.
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

## Experimental Setup

All experiments were run on Linux servers with NVIDIA A100 GPUs. The sampling
scripts use the following grid unless otherwise specified:

- Random seeds: `2024`, `2025`, `2026`
- Sampling temperatures: `0.6`, `1.0`, `1.4`

## Environment

The codebase was developed with Python 3.9, PyTorch 1.13.0, and CUDA 11.6.
Create the conda environment with:

```bash
conda env create -f environment.yaml -n AmpDiff
conda activate AmpDiff
```

The environment file installs the main conda and pip dependencies.

Activate the created environment before running the commands below. If your
local scorer paths differ from the defaults, pass the corresponding
`--pepnet_root`, `--hemopi2_root`, `--amppred_root`, `--iamp_root`, or
`--unidl_root` arguments.

Install optional external predictors in the top-level directory (sibling to
the `AmpDiff` folder) when using predictor-guided AMP evaluation:

```text
parent-directory/
├── AmpDiff/
├── PepNet/
├── AMPpred-MFA/
├── iAMP-Attenpred/
├── UniDL4BioPep/
└── hemopi2/
```

Predictor repositories:

- [PepNet](https://github.com/Harkool/PepNet)
- [AMPpred-MFA](https://github.com/Jiangle525/AMPpred-MFA)
- [iAMP-Attenpred](https://github.com/xingwxzz/iAMP-Attenpred)
- [UniDL4BioPep](https://github.com/dzjxzyd/UniDL4BioPep)
- [HemoPI2](https://github.com/raghavagps/hemopi2)

## Data Processing

1. Download the [raw dataset](https://doi.org/10.5281/zenodo.24358304) and
   place the extracted files under `raw dataset/`.
2. Process the raw AMP records:

```bash
python dataset/prepare_ampdiff_datasets.py \
  --raw-dir "raw dataset" \
  --out-dir data \
  --work-dir dataset/work \
  --motif-dir data/motif
```

The processed datasets are written to:

- `data/pretrain/`: datasets for the two pretraining stages
- `data/finetune/`: supervised datasets for `act` and `bal` finetuning
- `data/test/`: held-out datasets for sampling and evaluation

FASTA versions are provided for the pretrain and test splits where needed. The
dataset split summary is available in `data/dataset_report.md`.

## Pretraining

AmpDiff uses a two-stage pretraining (`stage 1`/`stage 2`) workflow and provides two independent usage paths (`A`/`B`).

### Path A: Retraining

Use this path to retrain both pretraining stages from the processed datasets
under `data/pretrain/`. Complete the Data Processing workflow above first.

1. Run `stage 1` on the broad AMP pretraining set:

```bash
bash amp_scripts/run_pretrain_stage1_all.sh
```

2. Run `stage 2`, which initializes from the stage 1 checkpoint and trains on
   the active-motif subset:

```bash
bash amp_scripts/run_pretrain_stage2_motif.sh
```

The two stages write their canonical checkpoints to:

- `checkpoints/pretrain/pretrain_stage1.pt`
- `checkpoints/pretrain/pretrain_stage2.pt`

### Path B: Warm-start from checkpoints

Use this path when you want to skip pretraining and start from released
weights.

1. Download the [pretraining checkpoints](https://doi.org/10.5281/zenodo.24358304).
2. Place the checkpoint files at:

```text
checkpoints/pretrain/pretrain_stage1.pt
checkpoints/pretrain/pretrain_stage2.pt
```

After completing either Path A or Path B, continue to finetuning with
`checkpoints/pretrain/pretrain_stage2.pt`.

## Finetuning

Finetuning uses supervised AMP activity data and selectable activity
guidance predictors.

### Path A: Re-finetuning

Use this path to finetune the `act` and `bal` models again with the processed
datasets under `data/finetune/`. Complete either pretraining path first so that
`checkpoints/pretrain/pretrain_stage2.pt` is available.

1. Re-finetune the `act` model:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/ampdiff_finetune_act.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/pretrain_stage2.pt \
  --guidance_scorer pepnet \
  --mode act
```

2. Re-finetune the `bal` model with optional hemolysis guidance:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/ampdiff_finetune_bal.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/pretrain_stage2.pt \
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
