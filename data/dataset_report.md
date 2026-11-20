# AmpDiff Dataset Report

This report audits the datasets currently published under `data/`.

## Processing policy

- Valid peptide length: 5-100 aa; standard 20 amino acids only.
- Initial train/test split uses CD-HIT clusters at 0.8 identity.
- Training sets use CD-HIT at 0.9 identity for internal deduplication.
- Training sets are filtered against the test union with CD-HIT-2D at 0.6 identity.

## Published files

| Dataset | Rows | CSV | FASTA |
| --- | ---: | --- | --- |
| `pretrain_stage1` | 26057 | `data/pretrain/ampdiff_pretrain_stage1.csv` | `data/pretrain/ampdiff_pretrain_stage1.fasta` |
| `pretrain_stage2` | 8398 | `data/pretrain/ampdiff_pretrain_stage2.csv` | `data/pretrain/ampdiff_pretrain_stage2.fasta` |
| `finetune_act` | 1294 | `data/finetune/ampdiff_finetune_act.csv` | `data/finetune/ampdiff_finetune_act.fasta` |
| `finetune_bal` | 64 | `data/finetune/ampdiff_finetune_bal.csv` | `data/finetune/ampdiff_finetune_bal.fasta` |
| `test_act` | 247 | `data/test/ampdiff_test_act.csv` | `data/test/ampdiff_test_act.fasta` |
| `test_bal` | 28 | `data/test/ampdiff_test_bal.csv` | `data/test/ampdiff_test_bal.fasta` |

All CSV/FASTA pairs contain the same sequences in the same order. No published dataset contains duplicate sequences.

## Sequence relationships

```json
{
  "pretrain_stage1_vs_pretrain_stage2_exact_overlap": 8200,
  "finetune_act_vs_finetune_bal_exact_overlap": 61,
  "test_act_vs_test_bal_exact_overlap": 6,
  "pretrain_stage1_vs_test_union_exact_overlap": 0,
  "pretrain_stage2_vs_test_union_exact_overlap": 0,
  "finetune_act_vs_test_union_exact_overlap": 0,
  "finetune_bal_vs_test_union_exact_overlap": 0
}
```

## Split summaries

```json
{
  "pretrain_stage1": {
    "length": {
      "n": 26057,
      "min": 11,
      "p25": 18,
      "median": 21,
      "p75": 63,
      "p95": 91,
      "max": 100
    },
    "motif": {
      "n": 26057,
      "seq_with_active_motif": 8267,
      "seq_without_active_motif": 17790,
      "seq_with_hom_neg_motif": 412,
      "active_positions_total": 54994,
      "hom_neg_positions_total": 2940,
      "bal_fixed_positions_total": 57147,
      "active_fixed_fraction": {
        "p25": 0.0,
        "median": 0.0,
        "p75": 0.05263157894736842,
        "p95": 0.3333333333333333,
        "max": 1.0
      },
      "bal_fixed_fraction": {
        "p25": 0.0,
        "median": 0.0,
        "p75": 0.05714285714285714,
        "p95": 0.35714285714285715,
        "max": 1.0
      },
      "optimizable_act_positions": {
        "min": 0,
        "p25": 17,
        "median": 20,
        "p75": 58,
        "p95": 88,
        "max": 100
      },
      "optimizable_bal_positions": {
        "min": 0,
        "p25": 16,
        "median": 20,
        "p75": 58,
        "p95": 88,
        "max": 100
      }
    },
    "hemolysis_summary": {},
    "combo_class": {}
  },
  "pretrain_stage2": {
    "length": {
      "n": 8398,
      "min": 11,
      "p25": 20,
      "median": 49,
      "p75": 76,
      "p95": 95,
      "max": 100
    },
    "motif": {
      "n": 8398,
      "seq_with_active_motif": 8398,
      "seq_without_active_motif": 0,
      "seq_with_hom_neg_motif": 231,
      "active_positions_total": 55196,
      "hom_neg_positions_total": 1607,
      "bal_fixed_positions_total": 55997,
      "active_fixed_fraction": {
        "p25": 0.05970149253731343,
        "median": 0.14035087719298245,
        "p75": 0.25,
        "p95": 0.5333333333333333,
        "max": 0.972972972972973
      },
      "bal_fixed_fraction": {
        "p25": 0.06,
        "median": 0.14285714285714285,
        "p75": 0.25,
        "p95": 0.546875,
        "max": 1.0
      },
      "optimizable_act_positions": {
        "min": 1,
        "p25": 16,
        "median": 39,
        "p75": 70,
        "p95": 90,
        "max": 97
      },
      "optimizable_bal_positions": {
        "min": 0,
        "p25": 16,
        "median": 39,
        "p75": 70,
        "p95": 90,
        "max": 97
      }
    },
    "hemolysis_summary": {},
    "combo_class": {}
  },
  "finetune_act": {
    "length": {
      "n": 1294,
      "min": 11,
      "p25": 15,
      "median": 20,
      "p75": 30,
      "p95": 49,
      "max": 100
    },
    "motif": {
      "n": 1294,
      "seq_with_active_motif": 1294,
      "seq_without_active_motif": 0,
      "seq_with_hom_neg_motif": 76,
      "active_positions_total": 9704,
      "hom_neg_positions_total": 563,
      "bal_fixed_positions_total": 9980,
      "active_fixed_fraction": {
        "p25": 0.17391304347826086,
        "median": 0.2857142857142857,
        "p75": 0.4666666666666667,
        "p95": 0.8181818181818182,
        "max": 0.9791666666666666
      },
      "bal_fixed_fraction": {
        "p25": 0.18181818181818182,
        "median": 0.3,
        "p75": 0.5,
        "p95": 0.8333333333333334,
        "max": 0.9791666666666666
      },
      "optimizable_act_positions": {
        "min": 1,
        "p25": 9,
        "median": 13,
        "p75": 21,
        "p95": 43,
        "max": 95
      },
      "optimizable_bal_positions": {
        "min": 1,
        "p25": 8,
        "median": 13,
        "p75": 21,
        "p95": 43,
        "max": 95
      }
    },
    "hemolysis_summary": {
      "unknown": 392,
      "hemolytic": 449,
      "mixed": 453
    },
    "combo_class": {
      "continuous__unknown": 392,
      "continuous__continuous": 675,
      "continuous__hemolytic": 18,
      "continuous__mixed": 209
    }
  },
  "finetune_bal": {
    "length": {
      "n": 64,
      "min": 11,
      "p25": 16,
      "median": 20,
      "p75": 26,
      "p95": 31,
      "max": 42
    },
    "motif": {
      "n": 64,
      "seq_with_active_motif": 64,
      "seq_without_active_motif": 0,
      "seq_with_hom_neg_motif": 64,
      "active_positions_total": 447,
      "hom_neg_positions_total": 494,
      "bal_fixed_positions_total": 689,
      "active_fixed_fraction": {
        "p25": 0.1875,
        "median": 0.25,
        "p75": 0.4375,
        "p95": 0.76,
        "max": 0.9259259259259259
      },
      "bal_fixed_fraction": {
        "p25": 0.3793103448275862,
        "median": 0.47368421052631576,
        "p75": 0.6875,
        "p95": 0.8461538461538461,
        "max": 0.9473684210526315
      },
      "optimizable_act_positions": {
        "min": 2,
        "p25": 9,
        "median": 12,
        "p75": 20,
        "p95": 25,
        "max": 38
      },
      "optimizable_bal_positions": {
        "min": 1,
        "p25": 5,
        "median": 9,
        "p75": 15,
        "p95": 19,
        "max": 37
      }
    },
    "hemolysis_summary": {
      "hemolytic": 21,
      "mixed": 43
    },
    "combo_class": {
      "continuous__hemolytic": 1,
      "continuous__continuous": 55,
      "continuous__mixed": 8
    }
  },
  "test_act": {
    "length": {
      "n": 247,
      "min": 6,
      "p25": 13,
      "median": 18,
      "p75": 25,
      "p95": 42,
      "max": 84
    },
    "motif": {
      "n": 247,
      "seq_with_active_motif": 247,
      "seq_without_active_motif": 0,
      "seq_with_hom_neg_motif": 41,
      "active_positions_total": 1918,
      "hom_neg_positions_total": 373,
      "bal_fixed_positions_total": 2126,
      "active_fixed_fraction": {
        "p25": 0.20833333333333334,
        "median": 0.35294117647058826,
        "p75": 0.6,
        "p95": 0.8571428571428571,
        "max": 0.9285714285714286
      },
      "bal_fixed_fraction": {
        "p25": 0.23529411764705882,
        "median": 0.42105263157894735,
        "p75": 0.65625,
        "p95": 0.8888888888888888,
        "max": 0.9565217391304348
      },
      "optimizable_act_positions": {
        "min": 1,
        "p25": 6,
        "median": 10,
        "p75": 17,
        "p95": 36,
        "max": 80
      },
      "optimizable_bal_positions": {
        "min": 1,
        "p25": 5,
        "median": 9,
        "p75": 16,
        "p95": 36,
        "max": 80
      }
    },
    "hemolysis_summary": {
      "mixed": 140,
      "hemolytic": 107
    },
    "combo_class": {
      "continuous__continuous": 200,
      "continuous__mixed": 44,
      "continuous__hemolytic": 3
    }
  },
  "test_bal": {
    "length": {
      "n": 28,
      "min": 15,
      "p25": 18,
      "median": 24,
      "p75": 24,
      "p95": 24,
      "max": 24
    },
    "motif": {
      "n": 28,
      "seq_with_active_motif": 28,
      "seq_without_active_motif": 0,
      "seq_with_hom_neg_motif": 28,
      "active_positions_total": 159,
      "hom_neg_positions_total": 467,
      "bal_fixed_positions_total": 500,
      "active_fixed_fraction": {
        "p25": 0.16666666666666666,
        "median": 0.16666666666666666,
        "p75": 0.2222222222222222,
        "p95": 0.5625,
        "max": 0.9166666666666666
      },
      "bal_fixed_fraction": {
        "p25": 0.625,
        "median": 0.8888888888888888,
        "p75": 0.9583333333333334,
        "p95": 0.9583333333333334,
        "max": 0.9583333333333334
      },
      "optimizable_act_positions": {
        "min": 2,
        "p25": 11,
        "median": 20,
        "p75": 20,
        "p95": 20,
        "max": 20
      },
      "optimizable_bal_positions": {
        "min": 1,
        "p25": 1,
        "median": 2,
        "p75": 6,
        "p95": 9,
        "max": 18
      }
    },
    "hemolysis_summary": {
      "hemolytic": 13,
      "mixed": 15
    },
    "combo_class": {
      "continuous__continuous": 28
    }
  }
}
```
