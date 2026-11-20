#!/usr/bin/env python3
"""Build the published AmpDiff train sets with CD-HIT filtering."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

from audit_dataset_release import generate_reports


ROOT = Path(__file__).resolve().parents[1]
TRAIN_FILES = (
    "pretrain/ampdiff_pretrain_stage1.csv",
    "pretrain/ampdiff_pretrain_stage2.csv",
    "finetune/ampdiff_finetune_act.csv",
    "finetune/ampdiff_finetune_bal.csv",
)
TEST_FILES = (
    "test/ampdiff_test_act.csv",
    "test/ampdiff_test_bal.csv",
)


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_rows(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def write_fasta(seqs: Iterable[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for index, seq in enumerate(seqs):
            handle.write(f">seq_{index}|len={len(seq)}\n{seq}\n")


def read_fasta_sequences(path: Path) -> list[str]:
    sequences: list[str] = []
    chunks: list[str] = []
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


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)


def cdhit_word_size(identity: float) -> str:
    if identity >= 0.7:
        return "5"
    if identity >= 0.6:
        return "4"
    if identity >= 0.5:
        return "3"
    return "2"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=ROOT / "dataset" / "work" / "prepared")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--work-dir", type=Path, default=ROOT / "dataset" / "work" / "cdhit60")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    work_dir = args.work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    test_sequences: list[str] = []
    for relative_path in TEST_FILES:
        source_csv = input_dir / relative_path
        fields, rows = read_rows(source_csv)
        sequences = [row["sequence"] for row in rows]
        test_sequences.extend(sequences)
        output_csv = output_dir / relative_path
        write_rows(output_csv, fields, rows)
        write_fasta(sequences, output_csv.with_suffix(".fasta"))

    test_sequences = sorted(set(test_sequences))
    test_fasta = work_dir / "ampdiff_test_union.fasta"
    write_fasta(test_sequences, test_fasta)

    summary: list[dict[str, object]] = []
    for relative_path in TRAIN_FILES:
        source_csv = input_dir / relative_path
        fieldnames, rows = read_rows(source_csv)
        seq_to_row = {row["sequence"]: row for row in rows}

        raw_fasta = work_dir / f"{source_csv.stem}.fasta"
        internal_fasta = work_dir / f"{source_csv.stem}.cdhit90.fasta"
        final_prefix = work_dir / f"{source_csv.stem}.cdhit90_vs_test60"
        write_fasta(seq_to_row, raw_fasta)

        run([
            shutil.which("cd-hit") or "cd-hit",
            "-i", str(raw_fasta),
            "-o", str(internal_fasta),
            "-c", "0.9",
            "-n", cdhit_word_size(0.9),
            "-d", "0",
            "-M", "0",
            "-T", "0",
        ])
        run([
            shutil.which("cd-hit-2d") or "cd-hit-2d",
            "-i", str(test_fasta),
            "-i2", str(internal_fasta),
            "-o", str(final_prefix),
            "-c", "0.6",
            "-n", cdhit_word_size(0.6),
            "-d", "0",
            "-M", "0",
            "-T", "0",
        ])

        kept_sequences = read_fasta_sequences(final_prefix)
        kept_rows = [seq_to_row[seq] for seq in kept_sequences]
        output_csv = output_dir / relative_path
        written = write_rows(output_csv, fieldnames, kept_rows)
        write_fasta(kept_sequences, output_csv.with_suffix(".fasta"))
        summary.append({
            "input": str(source_csv),
            "output_csv": str(output_csv),
            "output_fasta": str(output_csv.with_suffix(".fasta")),
            "input_rows": len(rows),
            "after_internal_90": len(read_fasta_sequences(internal_fasta)),
            "after_vs_test_60": written,
        })

    generate_reports(output_dir)
    for item in summary:
        print(item)


if __name__ == "__main__":
    main()
