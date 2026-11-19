#!/usr/bin/env python3
"""
Prepare AmpDiff pretrain / finetune / test datasets from the merged raw AMP data.

Design notes:
- Pretraining follows HuDiff-Nb's unsupervised denoising setup: sequence only.
- Finetuning uses high-confidence experimental activity records aggregated per sequence.
- Test sequences are held out at a CD-HIT cluster level to reduce similarity leakage.
- Outputs are written next to existing project data without overwriting legacy files.
"""

from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import random
import shutil
import statistics
import subprocess
import tempfile
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

AA_MASS = {
    "A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13, "F": 165.19,
    "G": 75.07, "H": 155.16, "I": 131.17, "K": 146.19, "L": 131.17,
    "M": 149.21, "N": 132.12, "P": 115.13, "Q": 146.15, "R": 174.20,
    "S": 105.09, "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19,
}

ACTIVITY_METRICS = {"MIC", "MBC", "IC50", "EC50"}
ALLOWED_ACTIVITY_COMPARATORS = {"", "=", "~", "<", "<="}
PROSITE_KEYS = ("PROSITE", "Cysteine Regex")
MERCI_KEYS = ("MERCI",)


@dataclass
class MasterRecord:
    sequence: str
    sequence_length: int
    evidence_groups: str
    databases: str
    activity_summary: str
    hemolysis_summary: str
    combo_class: str
    evidence_rows: int


@dataclass
class ActivityRecord:
    sequence: str
    value: float
    activity_uM_min: float
    activity_uM_median: float
    activity_records: int
    activity_metrics: str
    activity_units: str
    databases: str
    hemolysis_summary: str
    combo_class: str
    sequence_length: int


def normalize_sequence(seq: str) -> str:
    return (seq or "").strip().upper()


def is_valid_sequence(seq: str, min_len: int, max_len: int) -> bool:
    return min_len <= len(seq) <= max_len and set(seq).issubset(VALID_AA)


def stable_sort_key(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()


def peptide_mw(seq: str) -> float:
    """Approximate peptide molecular weight in Da from average residue masses."""
    if not seq:
        return 0.0
    return sum(AA_MASS[aa] for aa in seq) - (len(seq) - 1) * 18.015


def parse_float(raw: str) -> Optional[float]:
    raw = (raw or "").strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def normalize_unit(unit: str) -> str:
    u = (unit or "").strip().lower()
    replacements = {
        "μ": "u",
        "µ": "u",
        "micro": "u",
        "碌": "u",
    }
    for old, new in replacements.items():
        u = u.replace(old, new)
    u = u.replace(" ", "")
    return u


def concentration_to_uM(value: float, unit: str, seq: str) -> Optional[float]:
    unit = normalize_unit(unit)
    if value <= 0:
        return None
    if unit in {"um", "umol/l", "umolar"}:
        return value
    if unit == "nm":
        return value / 1000.0
    if unit == "mm":
        return value * 1000.0
    if unit in {"ug/ml", "ug/ml.", "ug/ml)", "ug/ml,", "ug/mL".lower(), "mg/l", "mg/liter"}:
        mw = peptide_mw(seq)
        return value * 1000.0 / mw if mw > 0 else None
    if unit in {"mg/ml"}:
        mw = peptide_mw(seq)
        return value * 1_000_000.0 / mw if mw > 0 else None
    return None


def read_master(raw_dir: Path, min_len: int, max_len: int) -> Dict[str, MasterRecord]:
    records: Dict[str, MasterRecord] = {}
    path = raw_dir / "amp_peptide_master.csv"
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            seq = normalize_sequence(row.get("representative_sequence") or row.get("sequence_key"))
            if row.get("is_standard_aa_sequence") != "yes":
                continue
            if not is_valid_sequence(seq, min_len, max_len):
                continue
            records[seq] = MasterRecord(
                sequence=seq,
                sequence_length=len(seq),
                evidence_groups=row.get("evidence_groups", ""),
                databases=row.get("databases", ""),
                activity_summary=row.get("activity_binary_summary", ""),
                hemolysis_summary=row.get("hemolysis_binary_summary", ""),
                combo_class=row.get("combo_class", ""),
                evidence_rows=int(row.get("evidence_rows") or 0),
            )
    return records


def build_activity_records(raw_dir: Path, master: Dict[str, MasterRecord]) -> Dict[str, ActivityRecord]:
    values: Dict[str, List[Tuple[float, str, str, str]]] = defaultdict(list)
    path = raw_dir / "amp_confirmed_continuous.csv"
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            seq = normalize_sequence(row.get("sequence_key") or row.get("sequence_original"))
            if seq not in master:
                continue
            metric = (row.get("activity_metric_raw") or "").strip().upper()
            if metric not in ACTIVITY_METRICS:
                continue
            comparator = (row.get("activity_comparator_raw") or "").strip()
            if comparator not in ALLOWED_ACTIVITY_COMPARATORS:
                continue
            raw_value = parse_float(row.get("activity_continuous_value_raw", ""))
            if raw_value is None:
                continue
            uM = concentration_to_uM(raw_value, row.get("activity_unit_raw", ""), seq)
            if uM is None or not math.isfinite(uM) or uM <= 0:
                continue
            values[seq].append((uM, metric, row.get("activity_unit_raw", ""), row.get("database", "")))

    records: Dict[str, ActivityRecord] = {}
    for seq, seq_values in values.items():
        uM_values = sorted(v[0] for v in seq_values)
        min_uM = uM_values[0]
        median_uM = statistics.median(uM_values)
        rec = master[seq]
        records[seq] = ActivityRecord(
            sequence=seq,
            value=-math.log10(min_uM),
            activity_uM_min=min_uM,
            activity_uM_median=median_uM,
            activity_records=len(seq_values),
            activity_metrics="|".join(sorted({v[1] for v in seq_values})),
            activity_units="|".join(sorted({v[2] for v in seq_values if v[2]})),
            databases="|".join(sorted({v[3] for v in seq_values if v[3]})),
            hemolysis_summary=rec.hemolysis_summary,
            combo_class=rec.combo_class,
            sequence_length=rec.sequence_length,
        )
    return records


def write_fasta(seqs: Sequence[str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for i, seq in enumerate(seqs):
            handle.write(f">seq_{i}|len={len(seq)}\n{seq}\n")


def run_cdhit(seqs: Sequence[str], identity: float, work_dir: Path) -> Optional[Path]:
    cdhit = shutil.which("cd-hit")
    inp = work_dir / "all_valid_sequences.fasta"
    out = work_dir / "all_valid_sequences.cdhit"
    existing_clstr = out.with_suffix(out.suffix + ".clstr")
    if not cdhit:
        if existing_clstr.exists():
            return existing_clstr
        return None
    write_fasta(seqs, inp)
    word_size = 5 if identity >= 0.7 else 4
    cmd = [
        cdhit,
        "-i", str(inp),
        "-o", str(out),
        "-c", str(identity),
        "-n", str(word_size),
        "-d", "0",
        "-M", "0",
        "-T", "0",
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return out.with_suffix(out.suffix + ".clstr")


def parse_cdhit_clusters(clstr_path: Optional[Path], seqs: Sequence[str]) -> Dict[str, int]:
    if clstr_path is None or not clstr_path.exists():
        return {seq: i for i, seq in enumerate(seqs)}
    id_to_seq = {f"seq_{i}": seq for i, seq in enumerate(seqs)}
    seq_to_cluster: Dict[str, int] = {}
    current = -1
    with clstr_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">Cluster"):
                current = int(line.split()[-1])
                continue
            marker = ">"
            end = "..."
            if marker in line and end in line:
                seq_id = line.split(marker, 1)[1].split(end, 1)[0]
                seq_id = seq_id.split("|", 1)[0]
                if seq_id in id_to_seq:
                    seq_to_cluster[id_to_seq[seq_id]] = current
    for i, seq in enumerate(seqs):
        seq_to_cluster.setdefault(seq, max(seq_to_cluster.values(), default=-1) + 1 + i)
    return seq_to_cluster


def choose_test_clusters(
    eligible_test_sequences: Set[str],
    seq_to_cluster: Dict[str, int],
    test_frac: float,
    seed: int,
) -> Set[int]:
    clusters: Dict[int, List[str]] = defaultdict(list)
    for seq in eligible_test_sequences:
        clusters[seq_to_cluster[seq]].append(seq)
    target = max(1, int(round(len(eligible_test_sequences) * test_frac)))
    cluster_ids = sorted(clusters, key=lambda cid: (stable_sort_key(str(cid)), cid))
    rng = random.Random(seed)
    rng.shuffle(cluster_ids)

    selected: Set[int] = set()
    selected_count = 0
    for cid in cluster_ids:
        selected.add(cid)
        selected_count += len(clusters[cid])
        if selected_count >= target:
            break
    return selected


def summarize_lengths(seqs: Sequence[str]) -> Dict[str, Optional[float]]:
    if not seqs:
        return {k: None for k in ["n", "min", "p25", "median", "p75", "p95", "max"]}
    lengths = sorted(len(s) for s in seqs)
    def q(p: float) -> int:
        return lengths[int((len(lengths) - 1) * p)]
    return {
        "n": len(lengths),
        "min": q(0),
        "p25": q(0.25),
        "median": q(0.5),
        "p75": q(0.75),
        "p95": q(0.95),
        "max": q(1),
    }


def write_csv(path: Path, rows: Iterable[dict], fieldnames: Sequence[str]) -> int:
    count = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def prosite_to_regex(pattern: str) -> Optional[re.Pattern]:
    parts = pattern.split("-")
    regex_parts = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^[xX]\((\d+)(?:,(\d+))?\)$", part)
        if m:
            lo, hi = m.group(1), m.group(2)
            regex_parts.append(f".{{{lo},{hi}}}" if hi else f".{{{lo}}}")
            continue
        if part.lower() == "x":
            regex_parts.append(".")
            continue
        if part.startswith("[") and "]" in part:
            regex_parts.append(part)
            continue
        regex_parts.append(re.escape(part))
    try:
        return re.compile("".join(regex_parts))
    except re.error:
        return None


def load_motif_matchers(motif_dir: Path, enabled_types: Set[str]) -> Tuple[List[str], List[re.Pattern], List[str]]:
    active_exact: List[str] = []
    active_prosite: List[re.Pattern] = []
    active_path = motif_dir / "AMP_act_Motifs.csv"
    if active_path.exists():
        with active_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                motif = (row.get("pattern") or "").strip()
                type_family = (row.get("type") or "").strip()
                if not motif or not type_family:
                    continue
                is_prosite = any(k in type_family for k in PROSITE_KEYS)
                is_merci = any(k in type_family for k in MERCI_KEYS)
                if is_prosite:
                    if "prosite" in enabled_types:
                        pattern = prosite_to_regex(motif)
                        if pattern is not None:
                            active_prosite.append(pattern)
                elif is_merci:
                    if "merci" in enabled_types:
                        active_exact.append(motif)
                elif "regular" in enabled_types:
                    active_exact.append(motif)

    hemolytic: List[str] = []
    hem_path = motif_dir / "AMP_hom_negMotifs.csv"
    if hem_path.exists():
        with hem_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                motif = (row.get("motif") or "").strip()
                if motif:
                    hemolytic.append(motif)
    return active_exact, active_prosite, hemolytic


def find_positions(seq: str, exact_motifs: Sequence[str], patterns: Sequence[re.Pattern]) -> Set[int]:
    covered: Set[int] = set()
    for motif in exact_motifs:
        start = 0
        while True:
            idx = seq.find(motif, start)
            if idx == -1:
                break
            covered.update(range(idx, idx + len(motif)))
            start = idx + 1
    for pattern in patterns:
        for match in pattern.finditer(seq):
            covered.update(range(match.start(), match.end()))
    return covered



def known_hemolysis(summary: str) -> bool:
    return summary not in {"", "unknown"}


def motif_counts(
    seq: str,
    active_exact: Sequence[str],
    active_prosite: Sequence[re.Pattern],
    hom_neg_motifs: Sequence[str],
) -> Dict[str, float]:
    active_pos = find_positions(seq, active_exact, active_prosite)
    hom_neg_pos = find_positions(seq, hom_neg_motifs, [])
    bal_fixed = active_pos | hom_neg_pos
    return {
        "length": len(seq),
        "active_positions": len(active_pos),
        "hom_neg_positions": len(hom_neg_pos),
        "bal_fixed_positions": len(bal_fixed),
        "optimizable_act_positions": len(seq) - len(active_pos),
        "optimizable_bal_positions": len(seq) - len(bal_fixed),
        "active_fixed_fraction": len(active_pos) / len(seq) if seq else 0.0,
        "bal_fixed_fraction": len(bal_fixed) / len(seq) if seq else 0.0,
        "has_active_motif": int(bool(active_pos)),
        "has_hom_neg_motif": int(bool(hom_neg_pos)),
    }


def profile_flags(
    seq: str,
    active_exact: Sequence[str],
    active_prosite: Sequence[re.Pattern],
    hom_neg_motifs: Sequence[str],
) -> Dict[str, bool]:
    profile = motif_counts(seq, active_exact, active_prosite, hom_neg_motifs)
    return {
        "active": bool(profile["has_active_motif"]),
        "hom_neg": bool(profile["has_hom_neg_motif"]),
        "optimizable_act": profile["optimizable_act_positions"] >= 1,
        "optimizable_bal": profile["optimizable_bal_positions"] >= 1,
    }


def split_task(
    candidates: Sequence[str],
    seq_to_cluster: Dict[str, int],
    test_frac: float,
    seed: int,
) -> Tuple[List[str], List[str], Set[int]]:
    selected_clusters = choose_test_clusters(set(candidates), seq_to_cluster, test_frac, seed)
    test = sorted(
        (seq for seq in candidates if seq_to_cluster[seq] in selected_clusters),
        key=lambda s: (seq_to_cluster[s], len(s), stable_sort_key(s)),
    )
    train = sorted(
        (seq for seq in candidates if seq_to_cluster[seq] not in selected_clusters),
        key=lambda s: (len(s), stable_sort_key(s)),
    )
    return train, test, selected_clusters


def summarize_task(
    seqs: Sequence[str],
    activity: Dict[str, ActivityRecord],
    active_exact: Sequence[str],
    active_prosite: Sequence[re.Pattern],
    hom_neg_motifs: Sequence[str],
) -> Dict[str, object]:
    profiles = [motif_counts(seq, active_exact, active_prosite, hom_neg_motifs) for seq in seqs]

    def q(values: Sequence[float], p: float) -> float:
        if not values:
            return 0.0
        ordered = sorted(values)
        return ordered[int((len(ordered) - 1) * p)]

    return {
        "length": summarize_lengths(seqs),
        "motif": {
            "n": len(profiles),
            "seq_with_active_motif": sum(1 for p in profiles if p["has_active_motif"]),
            "seq_without_active_motif": sum(1 for p in profiles if not p["has_active_motif"]),
            "seq_with_hom_neg_motif": sum(1 for p in profiles if p["has_hom_neg_motif"]),
            "active_positions_total": int(sum(p["active_positions"] for p in profiles)),
            "hom_neg_positions_total": int(sum(p["hom_neg_positions"] for p in profiles)),
            "bal_fixed_positions_total": int(sum(p["bal_fixed_positions"] for p in profiles)),
            "active_fixed_fraction": {
                "p25": q([p["active_fixed_fraction"] for p in profiles], 0.25),
                "median": q([p["active_fixed_fraction"] for p in profiles], 0.5),
                "p75": q([p["active_fixed_fraction"] for p in profiles], 0.75),
                "p95": q([p["active_fixed_fraction"] for p in profiles], 0.95),
                "max": q([p["active_fixed_fraction"] for p in profiles], 1.0),
            },
            "bal_fixed_fraction": {
                "p25": q([p["bal_fixed_fraction"] for p in profiles], 0.25),
                "median": q([p["bal_fixed_fraction"] for p in profiles], 0.5),
                "p75": q([p["bal_fixed_fraction"] for p in profiles], 0.75),
                "p95": q([p["bal_fixed_fraction"] for p in profiles], 0.95),
                "max": q([p["bal_fixed_fraction"] for p in profiles], 1.0),
            },
            "optimizable_act_positions": {
                "min": q([p["optimizable_act_positions"] for p in profiles], 0.0),
                "p25": q([p["optimizable_act_positions"] for p in profiles], 0.25),
                "median": q([p["optimizable_act_positions"] for p in profiles], 0.5),
                "p75": q([p["optimizable_act_positions"] for p in profiles], 0.75),
                "p95": q([p["optimizable_act_positions"] for p in profiles], 0.95),
                "max": q([p["optimizable_act_positions"] for p in profiles], 1.0),
            },
            "optimizable_bal_positions": {
                "min": q([p["optimizable_bal_positions"] for p in profiles], 0.0),
                "p25": q([p["optimizable_bal_positions"] for p in profiles], 0.25),
                "median": q([p["optimizable_bal_positions"] for p in profiles], 0.5),
                "p75": q([p["optimizable_bal_positions"] for p in profiles], 0.75),
                "p95": q([p["optimizable_bal_positions"] for p in profiles], 0.95),
                "max": q([p["optimizable_bal_positions"] for p in profiles], 1.0),
            },
        },
        "hemolysis_summary": Counter(activity[seq].hemolysis_summary for seq in seqs if seq in activity),
        "combo_class": Counter(activity[seq].combo_class for seq in seqs if seq in activity),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="/mnt/wucy/WUCHUYA/AmpDiff/raw dataset")
    parser.add_argument("--out-dir", default="/mnt/wucy/WUCHUYA/AmpDiff/data")
    parser.add_argument("--work-dir", default="/mnt/wucy/WUCHUYA/AmpDiff/dataset/work")
    parser.add_argument("--min-len", type=int, default=5)
    parser.add_argument("--max-len", type=int, default=100)
    parser.add_argument("--identity", type=float, default=0.8)
    parser.add_argument("--test-frac", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--motif-dir", default="/mnt/wucy/WUCHUYA/AmpDiff/data/motif")
    parser.add_argument("--motif-type", default="prosite,regular,merci")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    work_dir = Path(args.work_dir)
    motif_dir = Path(args.motif_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    for sub in ["pretrain", "finetune", "test"]:
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    enabled_motif_types = {
        item.strip().lower()
        for item in args.motif_type.split(",")
        if item.strip() and item.strip().lower() != "none"
    }
    master = read_master(raw_dir, args.min_len, args.max_len)
    activity = build_activity_records(raw_dir, master)
    active_exact, active_prosite, hom_neg_motifs = load_motif_matchers(
        motif_dir, enabled_motif_types
    )

    all_valid = sorted(master, key=lambda s: (len(s), stable_sort_key(s)))
    clstr = run_cdhit(all_valid, args.identity, work_dir)
    seq_to_cluster = parse_cdhit_clusters(clstr, all_valid)

    task_flags = {
        seq: profile_flags(seq, active_exact, active_prosite, hom_neg_motifs)
        for seq in activity
    }
    act_candidates = [
        seq for seq in activity
        if task_flags[seq]["active"] and task_flags[seq]["optimizable_act"]
    ]
    act_test_candidates = [
        seq for seq in act_candidates
        if known_hemolysis(activity[seq].hemolysis_summary)
    ]
    bal_candidates = [
        seq for seq, rec in activity.items()
        if task_flags[seq]["active"]
        and task_flags[seq]["hom_neg"]
        and task_flags[seq]["optimizable_bal"]
        and known_hemolysis(rec.hemolysis_summary)
    ]

    _, test_act, test_act_clusters = split_task(
        act_test_candidates, seq_to_cluster, args.test_frac, args.seed
    )
    finetune_act = sorted(
        (seq for seq in act_candidates if seq_to_cluster[seq] not in test_act_clusters),
        key=lambda s: (len(s), stable_sort_key(s)),
    )
    finetune_bal, test_bal, test_bal_clusters = split_task(
        bal_candidates, seq_to_cluster, args.test_frac, args.seed + 1
    )
    heldout_clusters = test_act_clusters | test_bal_clusters

    finetune_act_set = set(finetune_act)
    finetune_bal_set = set(finetune_bal)
    test_act_set = set(test_act)
    test_bal_set = set(test_bal)
    pretrain_all = sorted(
        (
            seq for seq, rec in master.items()
            if rec.activity_summary in {"positive", "continuous", "mixed"}
            and seq_to_cluster[seq] not in heldout_clusters
            and seq not in finetune_act_set
            and seq not in finetune_bal_set
            and seq not in test_act_set
            and seq not in test_bal_set
        ),
        key=lambda s: (len(s), stable_sort_key(s)),
    )
    pretrain_motif = [
        seq for seq in pretrain_all
        if (
            (profile := motif_counts(seq, active_exact, active_prosite, hom_neg_motifs))["has_active_motif"]
            and profile["optimizable_act_positions"] >= 1
        )
    ]

    finetune_fields = list(asdict(next(iter(activity.values()))).keys()) if activity else ["sequence", "value"]
    outputs = {
        "pretrain_all": out_dir / "pretrain" / "ampdiff_pretrain_all.csv",
        "pretrain": out_dir / "pretrain" / "ampdiff_pretrain.csv",
        "pretrain_fasta": out_dir / "pretrain" / "ampdiff_pretrain.fasta",
        "pretrain_motif": out_dir / "pretrain" / "ampdiff_pretrain_motif.csv",
        "pretrain_motif_fasta": out_dir / "pretrain" / "ampdiff_pretrain_motif.fasta",
        "finetune_act": out_dir / "finetune" / "ampdiff_finetune_act.csv",
        "finetune_bal": out_dir / "finetune" / "ampdiff_finetune_bal.csv",
        "test_act": out_dir / "test" / "ampdiff_test_act.csv",
        "test_act_fasta": out_dir / "test" / "ampdiff_test_act.fasta",
        "test_bal": out_dir / "test" / "ampdiff_test_bal.csv",
        "test_bal_fasta": out_dir / "test" / "ampdiff_test_bal.fasta",
        "report_json": out_dir / "dataset_report.json",
        "report_md": out_dir / "dataset_report.md",
    }

    write_csv(outputs["pretrain_all"], ({"sequence": seq} for seq in pretrain_all), ["sequence"])
    write_csv(outputs["pretrain"], ({"sequence": seq} for seq in pretrain_all), ["sequence"])
    write_fasta(pretrain_all, outputs["pretrain_fasta"])
    write_csv(outputs["pretrain_motif"], ({"sequence": seq} for seq in pretrain_motif), ["sequence"])
    write_fasta(pretrain_motif, outputs["pretrain_motif_fasta"])
    write_csv(outputs["finetune_act"], (asdict(activity[seq]) for seq in finetune_act), finetune_fields)
    write_csv(outputs["finetune_bal"], (asdict(activity[seq]) for seq in finetune_bal), finetune_fields)
    write_csv(outputs["test_act"], (asdict(activity[seq]) for seq in test_act), finetune_fields)
    write_fasta(test_act, outputs["test_act_fasta"])
    write_csv(outputs["test_bal"], (asdict(activity[seq]) for seq in test_bal), finetune_fields)
    write_fasta(test_bal, outputs["test_bal_fasta"])

    split_clusters = {
        "pretrain_all": {seq_to_cluster[s] for s in pretrain_all},
        "pretrain_motif": {seq_to_cluster[s] for s in pretrain_motif},
        "finetune_act": {seq_to_cluster[s] for s in finetune_act},
        "finetune_bal": {seq_to_cluster[s] for s in finetune_bal},
        "test_act": {seq_to_cluster[s] for s in test_act},
        "test_bal": {seq_to_cluster[s] for s in test_bal},
    }
    report = {
        "parameters": vars(args),
        "policy": {
            "pretrain": "broad AMP-like sequence-only pool; no active-motif requirement",
            "pretrain_motif": "active motif hit and >=1 optimizable residue",
            "finetune_act": "activity-supervised; active motif hit and >=1 optimizable residue",
            "finetune_bal": "activity-supervised; active motif + hom_neg motif hits; known hemolysis; >=1 optimizable residue",
            "test_act": "held-out act-valid clusters",
            "test_bal": "held-out bal-valid clusters",
        },
        "motifs": {
            "enabled_types": sorted(enabled_motif_types),
            "active_exact": len(active_exact),
            "active_prosite": len(active_prosite),
            "hom_neg": len(hom_neg_motifs),
        },
        "raw_valid_master_sequences": len(master),
        "activity_supervised_sequences": len(activity),
        "candidate_counts": {
            "act": len(act_candidates),
            "act_test_known_hemolysis": len(act_test_candidates),
            "bal": len(bal_candidates),
        },
        "outputs": {key: str(value) for key, value in outputs.items()},
        "splits": {
            "pretrain_all": summarize_task(pretrain_all, activity, active_exact, active_prosite, hom_neg_motifs),
            "pretrain_motif": summarize_task(pretrain_motif, activity, active_exact, active_prosite, hom_neg_motifs),
            "finetune_act": summarize_task(finetune_act, activity, active_exact, active_prosite, hom_neg_motifs),
            "finetune_bal": summarize_task(finetune_bal, activity, active_exact, active_prosite, hom_neg_motifs),
            "test_act": summarize_task(test_act, activity, active_exact, active_prosite, hom_neg_motifs),
            "test_bal": summarize_task(test_bal, activity, active_exact, active_prosite, hom_neg_motifs),
        },
        "cluster_summary": {
            "identity": args.identity,
            "test_act_vs_finetune_act_overlap": len(split_clusters["test_act"] & split_clusters["finetune_act"]),
            "test_bal_vs_finetune_bal_overlap": len(split_clusters["test_bal"] & split_clusters["finetune_bal"]),
            "test_act_vs_pretrain_overlap": len(split_clusters["test_act"] & split_clusters["pretrain_all"]),
            "test_bal_vs_pretrain_overlap": len(split_clusters["test_bal"] & split_clusters["pretrain_all"]),
            "test_act_clusters": len(test_act_clusters),
            "test_bal_clusters": len(test_bal_clusters),
            "heldout_cluster_union": len(heldout_clusters),
        },
    }
    outputs["report_json"].write_text(
        json.dumps(report, ensure_ascii=False, indent=2, default=lambda x: dict(x)),
        encoding="utf-8",
    )

    md = [
        "# AmpDiff Dataset Report",
        "",
        "## Split policy",
        f"- Valid peptide length: {args.min_len}-{args.max_len} aa; standard 20 amino acids only.",
        "- Pretrain keeps broad AMP-like sequences and does not require active motif hits.",
        "- Pretrain motif stage requires active motif hits and at least one optimizable residue.",
        "- finetune_act/test_act require active motif hits and at least one optimizable residue.",
        "- finetune_bal/test_bal require active motif + hom_neg motif hits, known hemolysis labels, and at least one optimizable residue.",
        f"- CD-HIT cluster holdout identity: {args.identity}.",
        "",
        "## Counts",
        f"- Raw valid master sequences: {len(master)}",
        f"- Activity-supervised sequences: {len(activity)}",
        f"- act candidates: {len(act_candidates)}",
        f"- act test candidates with known hemolysis: {len(act_test_candidates)}",
        f"- bal candidates: {len(bal_candidates)}",
        f"- pretrain_all: {len(pretrain_all)}",
        f"- pretrain_motif: {len(pretrain_motif)}",
        f"- finetune_act: {len(finetune_act)}",
        f"- finetune_bal: {len(finetune_bal)}",
        f"- test_act: {len(test_act)}",
        f"- test_bal: {len(test_bal)}",
        "",
        "## Split summaries",
        "```json",
        json.dumps(report["splits"], indent=2, default=lambda x: dict(x)),
        "```",
        "",
        "## Similarity summary",
        "```json",
        json.dumps(report["cluster_summary"], indent=2),
        "```",
    ]
    outputs["report_md"].write_text("\n".join(md) + "\n", encoding="utf-8")

    print(json.dumps({
        "pretrain_all": len(pretrain_all),
        "pretrain_motif": len(pretrain_motif),
        "finetune_act": len(finetune_act),
        "finetune_bal": len(finetune_bal),
        "test_act": len(test_act),
        "test_bal": len(test_bal),
        "report": str(outputs["report_md"]),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
