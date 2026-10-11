# Activity guidance predictor selection

AmpDiff exposes four interchangeable activity predictors:

- `pepnet`
- `amppred_mfa`
- `iamp_attenpred`
- `unidl4biopep`

Exactly one is selected with `--guidance_scorer`. During sampling,
`--independent_scorers auto` selects the remaining three for post-generation
evaluation. Those evaluator models are held in a separate scorer instance and
are never called by `combined_guidance_loss`.

Example:

```bash
python amp_scripts/amp_sample.py \
  --ckpt checkpoints/de/de.pt \
  --input_fasta /path/to/input.fasta \
  --guidance_scorer amppred_mfa \
  --independent_scorers auto \
  --mode de
```

The resulting `scores.tsv` contains scores from `pepnet`,
`iamp_attenpred`, and `unidl4biopep`. Add `--eval-hemolysis` when HemoPI2
is also required. It records `amppred_mfa` in the
`guidance_activity_scorer` column but does not reuse its score as independent
evidence.

Fine-tuning uses the same selection parameter:

```bash
python amp_scripts/amp_finetune.py \
  --data_path data/finetune/finetune.csv \
  --config_path configs/amp_finetune.yml \
  --pretrain_ckpt checkpoints/pretrain/stage1_all.pt \
  --guidance_scorer iamp_attenpred \
  --mode de
```

HemoPI2 remains a separate optional hemolysis objective in `inp` mode. Enable
it with `--hemolysis-guidance` in an environment that provides
`transformers`. It is not one of the four activity predictors.
