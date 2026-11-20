#!/usr/bin/env python3
"""Validate the published AmpDiff datasets and rebuild their reports."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Sequence, Tuple

from prepare_ampdiff_datasets import load_motif_matchers, summarize_task


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "data"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

DATASETS = {
    "pretrain_stage1": "pretrain/ampdiff_pretrain_stage1.csv",
    "pretrain_stage2": "pretrain/ampdiff_pretrain_stage2.csv",
    "finetune_act": "finetune/ampdiff_finetune_act.csv",
    "finetune_bal": "finetune/ampdiff_finetune_bal.csv",
    "test_act": "test/ampdiff_test_act.csv",
    "test_bal": "test/ampdiff_test_bal.csv",
}

POLICY = {
    "pretrain_stage1": "broad AMP-like sequence pool for stage 1 pretraining",
    "pretrain_stage2": "active-motif sequence pool for stage 2 pretraining",
    "finetune_act": "activity-supervised act-mode finetuning set",
    "finetune_bal": "activity- and hemolysis-supervised bal-mode finetuning set",
    "test_act": "held-out act-mode test set",
    "test_bal": "held-out bal-mode test set",
}


def read_csv(path: Path) -> Tuple[List[str], List[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if "sequence" not in fields:
        raise ValueError(f"Missing sequence column: {path}")
    return fields, rows


def read_fasta(path: Path) -> List[str]:
    sequences: List[str] = []
    chunks: List[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if chunks:
                    sequences.append("".join(chunks))
                    chunks = []
            else:
                chunks.append(line)
    if chunks:
        sequences.append("".join(chunks))
    return sequences


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def generate_reports(data_dir: Path = DEFAULT_DATA_DIR) -> dict:
    motif_dir = data_dir / "motif"
    active_exact, active_prosite, hom_neg_motifs = load_motif_matchers(
        motif_dir, {"prosite", "regular", "merci"}
    )

    dataset_sequences: Dict[str, List[str]] = {}
    split_summaries: Dict[str, object] = {}
    output_files: Dict[str, dict[str, str]] = {}
    validation: Dict[str, dict[str, object]] = {}

    for name, relative_csv in DATASETS.items():
        csv_path = data_dir / relative_csv
        fasta_path = csv_path.with_suffix(".fasta")
        fields, rows = read_csv(csv_path)
        sequences = [row["sequence"].strip() for row in rows]
        fasta_sequences = read_fasta(fasta_path)
        if sequences != fasta_sequences:
            raise ValueError(f"CSV/FASTA sequence mismatch: {csv_path} vs {fasta_path}")
        duplicate_count = len(sequences) - len(set(sequences))
        if duplicate_count:
            raise ValueError(f"Duplicate sequences in {csv_path}: {duplicate_count}")
        invalid_sequences = [
            seq for seq in sequences
            if not 5 <= len(seq) <= 100 or not set(seq).issubset(VALID_AA)
        ]
        if invalid_sequences:
            raise ValueError(f"Invalid sequences in {csv_path}: {len(invalid_sequences)}")

        activity = {
            row["sequence"].strip(): SimpleNamespace(
                hemolysis_summary=row.get("hemolysis_summary", ""),
                combo_class=row.get("combo_class", ""),
            )
            for row in rows
            if "hemolysis_summary" in fields or "combo_class" in fields
        }
        dataset_sequences[name] = sequences
        split_summaries[name] = summarize_task(
            sequences, activity, active_exact, active_prosite, hom_neg_motifs
        )
        output_files[name] = {
            "csv": relative(csv_path),
            "fasta": relative(fasta_path),
        }
        validation[name] = {
            "rows": len(sequences),
            "unique_sequences": len(set(sequences)),
            "duplicate_sequences": duplicate_count,
            "csv_fasta_identical_order": True,
            "fields": fields,
        }

    sequence_sets = {name: set(seqs) for name, seqs in dataset_sequences.items()}
    test_union = sequence_sets["test_act"] | sequence_sets["test_bal"]
    relationships = {
        "pretrain_stage1_vs_pretrain_stage2_exact_overlap": len(
            sequence_sets["pretrain_stage1"] & sequence_sets["pretrain_stage2"]
        ),
        "finetune_act_vs_finetune_bal_exact_overlap": len(
            sequence_sets["finetune_act"] & sequence_sets["finetune_bal"]
        ),
        "test_act_vs_test_bal_exact_overlap": len(
            sequence_sets["test_act"] & sequence_sets["test_bal"]
        ),
        "pretrain_stage1_vs_test_union_exact_overlap": len(
            sequence_sets["pretrain_stage1"] & test_union
        ),
        "pretrain_stage2_vs_test_union_exact_overlap": len(
            sequence_sets["pretrain_stage2"] & test_union
        ),
        "finetune_act_vs_test_union_exact_overlap": len(
            sequence_sets["finetune_act"] & test_union
        ),
        "finetune_bal_vs_test_union_exact_overlap": len(
            sequence_sets["finetune_bal"] & test_union
        ),
    }

    report = {
        "parameters": {
            "data_dir": relative(data_dir),
            "motif_dir": relative(motif_dir),
            "valid_length": [5, 100],
            "split_cluster_identity": 0.8,
            "training_internal_dedup_identity": 0.9,
            "training_vs_test_filter_identity": 0.6,
            "seed": 2026,
            "motif_types": ["merci", "prosite", "regular"],
        },
        "policy": POLICY,
        "motifs": {
            "active_exact": len(active_exact),
            "active_prosite": len(active_prosite),
            "hom_neg": len(hom_neg_motifs),
        },
        "outputs": output_files,
        "validation": validation,
        "splits": split_summaries,
        "sequence_relationships": relationships,
    }

    json_path = data_dir / "dataset_report.json"
    json_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, default=dict) + "\n",
        encoding="utf-8",
    )

    md = [
        "# AmpDiff Dataset Report",
        "",
        "This report audits the datasets currently published under `data/`.",
        "",
        "## Processing policy",
        "",
        "- Valid peptide length: 5-100 aa; standard 20 amino acids only.",
        "- Initial train/test split uses CD-HIT clusters at 0.8 identity.",
        "- Training sets use CD-HIT at 0.9 identity for internal deduplication.",
        "- Training sets are filtered against the test union with CD-HIT-2D at 0.6 identity.",
        "",
        "## Published files",
        "",
        "| Dataset | Rows | CSV | FASTA |",
        "| --- | ---: | --- | --- |",
    ]
    for name in DATASETS:
        info = validation[name]
        files = output_files[name]
        md.append(f"| `{name}` | {info['rows']} | `{files['csv']}` | `{files['fasta']}` |")
    md.extend([
        "",
        "All CSV/FASTA pairs contain the same sequences in the same order. No published dataset contains duplicate sequences.",
        "",
        "## Sequence relationships",
        "",
        "```json",
        json.dumps(relationships, indent=2),
        "```",
        "",
        "## Split summaries",
        "",
        "```json",
        json.dumps(split_summaries, indent=2, default=dict),
        "```",
    ])
    (data_dir / "dataset_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    generated = generate_reports()
    counts = {name: item["rows"] for name, item in generated["validation"].items()}
    print(json.dumps(counts, indent=2))
