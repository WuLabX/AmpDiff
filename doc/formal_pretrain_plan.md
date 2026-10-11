# AmpDiff formal pretraining plan

This plan reflects the current AmpDiff objective:

- Stage 1 learns the broad AMP-like sequence distribution.
- Stage 2 adapts the pretrained model to the motif-conditioned `de` task.
- PepNet / HemoPI2 guidance is reserved for downstream `de` / `inp` fine-tuning, not base pretraining.

## 1. Rebuild datasets

Run from `/mnt/wucy/WUCHUYA/AmpDiff`:

```bash
python dataset/prepare_ampdiff_datasets.py
```

Important outputs:

- `data/pretrain/ampdiff_pretrain_all.csv`
- `data/pretrain/ampdiff_pretrain_motif.csv`
- `data/finetune/ampdiff_finetune_de.csv, ampdiff_finetune_inp.csv`
- `data/test/ampdiff_test_de.csv, ampdiff_test_inp.csv`
- `data/dataset_report.md`

Current expected split after CD-HIT guarded rebuild:

- `pretrain_all`: 52,689 sequences.
- `pretrain_motif`: 18,493 sequences, all with active motif hits and at least one optimizable position.
- `finetune`: 7,415 sequences.
- `test`: 565 sequences.
- Test cluster overlap with pretrain / pretrain_motif / finetune: 0 at CD-HIT identity 0.8.

## 2. Stage 1 pretraining

Use all AMP-like pretraining sequences:

```bash
bash amp_scripts/run_pretrain_stage1_all.sh
```

Canonical best checkpoint:

```text
checkpoints/pretrain/stage1_all.pt
```

TensorBoard:

```bash
tensorboard --logdir logs/amp_pretrain --host 0.0.0.0 --port 6006
```

Primary curves to watch:

- `val/loss`
- `val/nll`
- `val/acc`
- `val/top3_acc`
- `val/top5_acc`
- `val/auc`
- `train/masked_per_seq`
- `train/motif_fixed_fraction`
- `train/no_motif_fraction`
- `train/optimizer_updates`

## 3. Stage 2 motif-conditioned pretraining

Start from the best Stage 1 checkpoint:

```bash
bash amp_scripts/run_pretrain_stage2_motif.sh checkpoints/pretrain/stage1_all.pt
```

Canonical best checkpoint:

```text
checkpoints/pretrain/stage2_motif.pt
```

For Stage 2, `train/no_motif_fraction` should stay at 0.0 because the dataset is explicitly motif-conditioned.

## 4. Why `train_loss=fr`

The fixed active motif is the AMP analogue of HuDiff-Nb CDR. The pretraining objective should mainly learn to reconstruct optimizable non-fixed positions around fixed motifs.

Using `train_loss=all` adds fixed motif reconstruction loss. That can be useful as a diagnostic, but it is not the cleanest objective for the formal motif-conditioned pretraining task.

## 5. Smoke test

A tiny config exists only for code-path verification:

```bash
python amp_scripts/amp_train.py \
  --data_path data/pretrain/ampdiff_pretrain_motif.csv \
  --config_path configs/amp_pretrain_smoke.yml \
  --log_path logs/amp_pretrain_smoke \
  --train_loss fr \
  --motif_type prosite,regular,merci \
  --mode de \
  --scorer none \
  --require_active_motif True \
  --min_optimizable 1 \
  --ckpt_subdir smoke \
  --ckpt_name smoke.pt
```

This should not be used for scientific results.

## 6. Next decision gate

After Stage 1 and Stage 2 complete, choose the formal pretrain checkpoint by validation loss first, then inspect:

- whether `val/top3_acc` and `val/top5_acc` improve even if `val/acc` saturates;
- whether `val/auc` remains stable or improves;
- whether `train/masked_per_seq` is nonzero and stable;
- whether Stage 2 improves motif-conditioned generation metrics compared with Stage 1.

Only after this checkpoint is selected should `de` fine-tuning with PepNet guidance and `inp` fine-tuning with PepNet + HemoPI2 guidance begin.
