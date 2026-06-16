#!/usr/bin/env python3
"""Orchestrator: turn a biological dataset into a leakage-aware train/val/test split.

This is the entrypoint the skill drives. It selects a modality strategy, produces
cluster labels, assigns whole clusters to splits via the shared core, then reports
the honest diagnostics every run must surface: realized split fractions, the
train<->test nearest-neighbour similarity distribution, and a provenance manifest.

Wired modalities: `sequence` (MMseqs2 with a k-mer fallback) and `metadata` (the
group / metadata axis alone, for histology-style tile/slide/patient data). Other
modalities (structure, small_molecule, protein_ligand, temporal) plug in at
`_cluster_for_modality`. The group axis composes with ANY modality via `group_keys`.

CLI:
    # sequence split
    python -m scripts.split --modality sequence --fasta data.fasta \
        --train 0.8 --val 0.1 --test 0.1 --min-seq-id 0.3 --out splits/

    # group / metadata split (propose-then-confirm)
    python -m scripts.split --modality metadata --metadata histo.csv \
        --id-col tile_id --auto-detect-groups            # prints proposal, exits
    python -m scripts.split --modality metadata --metadata histo.csv \
        --id-col tile_id --group-col patient_id --out splits/   # splits
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Sequence

from . import cluster_sequence, group_detect, leakage_metrics
from .core import assign_clusters_to_splits, clusters_from_pairs


def _hash_items(items: Sequence) -> str:
    h = hashlib.sha256()
    for it in items:
        h.update(repr(it).encode())
        h.update(b"\x00")
    return h.hexdigest()[:16]


def read_fasta(path: str) -> tuple[list[str], list[str]]:
    ids, seqs = [], []
    cur_id, cur = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    ids.append(cur_id)
                    seqs.append("".join(cur))
                cur_id, cur = line[1:].split()[0], []
            elif line:
                cur.append(line)
    if cur_id is not None:
        ids.append(cur_id)
        seqs.append("".join(cur))
    return ids, seqs


def _tool_backend(name: str, actual_version, validated_version, extra: dict) -> dict:
    """Backend manifest entry recording the ACTUAL tool version, warning on drift from
    the validated one (clustering shifts across releases, so the pin is not assumed)."""
    backend = {"backend": name, "version": actual_version,
               "version_validated": validated_version, **extra}
    if actual_version and validated_version and actual_version != validated_version:
        backend["warning"] = (
            f"installed {name} {actual_version} differs from the validated "
            f"{validated_version}; clustering can shift across releases - re-validate or "
            "install the validated version (recorded so the split stays auditable)")
    return backend


def _cluster_for_modality(modality: str, ids, items, *, threshold: float, prefer_tool: bool) -> tuple[list[int], dict]:
    """Return (cluster_labels, backend_info) for a clustering-based modality.

    Handles `sequence`, `small_molecule`, and `structure`. Temporal and
    protein_ligand do not go through this path (see make_split).
    """
    if modality == "sequence":
        if prefer_tool and cluster_sequence.mmseqs_available():
            labels = cluster_sequence.mmseqs_cluster(ids, items, min_seq_id=threshold)
            return labels, _tool_backend(
                "mmseqs2", cluster_sequence.mmseqs_version(), cluster_sequence.PINNED_MMSEQS_VERSION,
                {"min_seq_id": threshold, "cluster_mode": 1})
        pairs = cluster_sequence.kmer_similar_pairs(items, threshold=threshold)
        labels = clusters_from_pairs(len(items), pairs)
        return labels, {"backend": "kmer_jaccard_fallback", "threshold": threshold, "k": 3,
                        "warning": "MMseqs2 not available; used k-mer Jaccard approximation. "
                                   "Install MMseqs2 for production splits."}

    if modality == "small_molecule":
        from . import split_small_molecule as sm
        if not sm.rdkit_available():
            raise RuntimeError("RDKit required for small_molecule modality; `pip install rdkit`")
        # Scaffold split is the baseline; Butina at the configured Tanimoto is harder.
        labels = sm.scaffold_labels(items)
        return labels, {"backend": "rdkit_bemis_murcko_scaffold", "note": sm.PINNED_RDKIT_NOTE}

    if modality == "structure":
        from . import cluster_structure as cs
        if not cs.foldseek_available():
            raise RuntimeError(
                "Foldseek required for structure modality and no pure-Python fallback exists "
                "(needs 3D coordinates). Install Foldseek, or supply precomputed cluster labels, "
                "or fall back to a sequence split and disclose that structural leakage is unchecked."
            )
        labels = cs.foldseek_cluster(items, ids=ids, tmscore_threshold=threshold)
        return labels, _tool_backend(
            "foldseek", cs.foldseek_version(), cs.PINNED_FOLDSEEK_VERSION,
            {"tmscore_threshold": threshold})

    raise NotImplementedError(f"modality '{modality}' not a clustering modality")


def _sequence_similarity(threshold_k: int = 3):
    return lambda a, b: cluster_sequence.kmer_jaccard(a, b, k=threshold_k)


def _similarity_for_modality(modality: str):
    """Return a similarity callable for scoring cross-split leakage, or None if not scorable here."""
    if modality == "sequence":
        return _sequence_similarity()
    if modality == "small_molecule":
        from . import split_small_molecule as sm
        if sm.rdkit_available():
            return sm.tanimoto
    return None


_METADATA_MODALITIES = ("metadata", "grouped")


def _grouping_block(group_keys, final_labels, modality, group_provenance) -> dict:
    """The honest 'group axis' report: sizes, nulls, and whether similarity glued
    distinct groups together (clusters_merged_across_groups)."""
    dist = group_detect.group_size_distribution(group_keys)
    keys_in_cluster: dict[int, set] = {}
    for lab, k in zip(final_labels, group_keys):
        if k is not None:
            keys_in_cluster.setdefault(lab, set()).add(k)
    merged = sum(1 for ks in keys_in_cluster.values() if len(ks) > 1)
    block = {
        "applied": True,
        "n_groups": dist["n_groups"],
        "group_size_distribution": {k: dist[k] for k in ("min", "max", "median", "p95", "n_singletons")},
        "n_null_group_keys": dist["n_null_group_keys"],
        "composed_with_similarity": modality not in _METADATA_MODALITIES,
        "clusters_merged_across_groups": merged,
    }
    if group_provenance:
        # Carry user-facing context (which column, who confirmed it, the proposal).
        for k in ("group_col", "id_col", "confirmed_by", "proposal"):
            if k in group_provenance:
                block[k] = group_provenance[k]
    return block


def make_split(
    ids: Sequence[str],
    items: Sequence,
    modality: str = "sequence",
    fractions: dict[str, float] | None = None,
    threshold: float = 0.3,
    seed: int = 0,
    stratify_labels: Sequence | None = None,
    prefer_tool: bool = True,
    group_keys: Sequence | None = None,
    group_provenance: dict | None = None,
) -> dict:
    """Produce a leakage-aware split and a full diagnostic report.

    Returns a dict with: per-item assignment, realized fractions, invariant
    reports, leakage stats (nearest-neighbour distribution), and a provenance
    manifest. The orchestrator never silently hides a failed invariant; callers
    inspect `report["invariants"]`.

    Two leakage axes, composed in a *single* union-find pass:

    * **similarity** - the per-modality clustering (sequence identity, scaffold, ...).
    * **group / metadata** - when ``group_keys`` is given (item-aligned), replicate
      units that share a source entity (patient, site, batch) are forced into one
      split. ``modality="metadata"`` runs the group axis ALONE (no similarity tool,
      no similarity leakage scored) for cases like histology where no biological
      similarity is available.
    """
    if fractions is None:
        fractions = {"train": 0.8, "val": 0.1, "test": 0.1}

    sim = None
    if modality in _METADATA_MODALITIES:
        if group_keys is None:
            raise ValueError("modality 'metadata' requires group_keys (the confirmed grouping column)")
        cluster_labels = clusters_from_pairs(len(items), group_detect.pairs_from_groups(group_keys))
        backend = {"backend": "group_metadata_only",
                   "note": "group axis only; no biological-similarity tool was run and no "
                           "similarity leakage is scored (there is no sequence/structure/molecule "
                           "to compare). Group integrity is the guarantee here."}
        if group_provenance:
            backend["group_provenance"] = {k: v for k, v in group_provenance.items() if k != "proposal"}
    else:
        cluster_labels, backend = _cluster_for_modality(
            modality, ids, items, threshold=threshold, prefer_tool=prefer_tool
        )
        if group_keys is not None:
            # Compose: lift similarity labels back to edges, add group edges, re-cluster
            # ONCE. Whole-component assignment then respects both axes simultaneously.
            sim_pairs = group_detect.pairs_from_labels(cluster_labels)
            grp_pairs = group_detect.pairs_from_groups(group_keys)
            cluster_labels = clusters_from_pairs(len(items), list(sim_pairs) + grp_pairs)
            backend = {**backend, "composed_with_groups": True}
        sim = _similarity_for_modality(modality)

    result = assign_clusters_to_splits(cluster_labels, fractions, seed=seed, stratify_labels=stratify_labels)

    invariants = [
        leakage_metrics.check_partition(result.assignment, len(items)),
        leakage_metrics.check_cluster_integrity(result.assignment, cluster_labels),
        leakage_metrics.check_ratio(result.realized_fractions, fractions),
        leakage_metrics.check_no_empty_splits(result.assignment, fractions),
    ]
    if group_keys is not None:
        invariants.append(leakage_metrics.check_group_integrity(result.assignment, group_keys))

    # Scaffold disjointness is a meaningful extra invariant for molecules.
    if modality == "small_molecule":
        from . import split_small_molecule as sm
        if sm.rdkit_available():
            scaffolds = [sm.bemis_murcko_scaffold(s) for s in items]
            invariants.append(leakage_metrics.check_scaffold_disjoint(scaffolds, result.assignment))

    # Leakage stats for any modality we can score with a similarity callable.
    stats = None
    if sim is not None:
        stats = leakage_metrics.compute_leakage_stats(items, result.assignment, sim, threshold)

    # Channel for advisories that are NOT hard invariant failures (empty splits are now
    # the no_empty_splits invariant; this stays for future soft notices, e.g. large
    # discard fractions). Kept so the report shape is stable for callers.
    warnings_out: list[str] = []

    manifest = {
        "modality": modality,
        "n_items": len(items),
        "n_clusters": result.n_clusters,
        "fractions_requested": fractions,
        "fractions_realized": result.realized_fractions,
        "threshold": threshold,
        "seed": seed,
        "backend": backend,
        "input_hash": _hash_items(items),
        "split_hash": _hash_items(result.assignment),
    }

    report = {
        "assignment": result.assignment,
        "ids": list(ids),
        "cluster_labels": cluster_labels,
        "realized_fractions": result.realized_fractions,
        "invariants": [asdict(r) for r in invariants],
        "all_invariants_pass": all(bool(r) for r in invariants),
        "leakage": None if stats is None else {
            "max_cross_similarity": stats.max_cross_similarity,
            "threshold": stats.threshold,
            "n_violations": stats.n_violations,
            "passes_threshold": stats.passes_threshold,
            "nn_similarity": stats.nn_similarity,
        },
        "provenance": manifest,
        "warnings": warnings_out,
    }
    if group_keys is not None:
        report["grouping"] = _grouping_block(group_keys, cluster_labels, modality, group_provenance)
    return report


def _print_proposal(proposal) -> None:
    """Render a ranked grouping proposal for the user to confirm (propose-then-confirm)."""
    print("Grouping-key proposal (propose-then-confirm: nothing is split yet)\n")
    rec = proposal.recommended
    if rec is not None:
        print(f"  RECOMMENDED: --group-col {rec.column}  "
              f"({rec.n_groups} groups, sizes {rec.group_size_min}-{rec.group_size_max}, "
              f"median {rec.group_size_median:g}; canonical='{rec.canonical_key}')")
    print("\n  Ranked candidates:")
    for c in proposal.candidates:
        tag = "DISQUALIFIED" if c.disqualified else f"score={c.rank_score:.2f}"
        sat = "satisfiable" if c.satisfiable_for else "NOT satisfiable for requested ratio"
        print(f"   - {c.column:<20} {tag:<14} groups={c.n_groups:<6} {sat}")
        for w in c.warnings:
            print(f"       ! {w}")
    for note in proposal.notes:
        print(f"\n  note: {note}")
    print("\nRe-run with --group-col <column> once you have confirmed the unit of generalization.")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Leakage-aware biological train/val/test split")
    p.add_argument("--modality", default="sequence",
                   choices=["sequence", "structure", "small_molecule", "protein_ligand",
                            "temporal", "metadata"])
    p.add_argument("--fasta", help="FASTA input for sequence modality")
    p.add_argument("--train", type=float, default=0.8)
    p.add_argument("--val", type=float, default=0.1)
    p.add_argument("--test", type=float, default=0.1)
    p.add_argument("--min-seq-id", type=float, default=0.3, dest="threshold",
                   help="identity/similarity threshold (sequence: MMseqs2 --min-seq-id)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="splits", help="output directory")
    # Group / metadata axis (composes with any modality; metadata = group axis only).
    p.add_argument("--metadata", help="CSV/TSV of per-item metadata for group-aware splitting")
    p.add_argument("--id-col", help="metadata column matching item ids (auto-guessed if omitted)")
    p.add_argument("--group-col", help="confirmed grouping column; no group key spans two splits")
    p.add_argument("--stratify-col", help="metadata column to stratify on (also disqualified as a group key)")
    p.add_argument("--auto-detect-groups", action="store_true",
                   help="print the ranked grouping proposal and EXIT without splitting "
                        "unless --group-col is also given (propose-then-confirm)")
    p.add_argument("--drop-null-groups", action="store_true",
                   help="drop items whose group key is null instead of treating them as singletons")
    p.add_argument("--no-tools", action="store_false", dest="prefer_tool",
                   help="force the hermetic k-mer / pure-Python fallback instead of MMseqs2/Foldseek "
                        "(deterministic; for reproducible tests and demos)")
    args = p.parse_args(argv)

    fractions = {k: v for k, v in {"train": args.train, "val": args.val, "test": args.test}.items() if v > 0}

    # ----- load items per modality -----
    if args.modality == "metadata":
        if not args.metadata:
            p.error("--metadata is required for modality metadata")
        table = group_detect.load_metadata(args.metadata)
        id_col = args.id_col or group_detect.guess_id_column(table)
        if id_col is None:
            p.error("could not guess an id column; pass --id-col")
        ids = [str(v) for v in table.column(id_col)]
        if len(set(ids)) != len(ids):
            p.error(f"id column '{id_col}' has duplicate values; a duplicate id is itself a finding")
        items = list(ids)  # no biological content in metadata-only mode
        aligned = table
    elif args.modality == "sequence":
        if not args.fasta:
            p.error("--fasta is required for sequence modality")
        ids, items = read_fasta(args.fasta)
        aligned = None
        if args.metadata:
            table = group_detect.load_metadata(args.metadata)
            id_col = args.id_col or group_detect.guess_id_column(table)
            if id_col is None:
                p.error("could not guess a metadata id column; pass --id-col")
            aligned = table.align_to(ids, id_col)
    else:
        p.error(f"modality '{args.modality}' not yet wired in the CLI")
        return 2

    # ----- group axis: propose-then-confirm -----
    group_keys = None
    group_provenance = None
    if args.modality == "metadata" or args.metadata:
        if args.auto_detect_groups and not args.group_col:
            proposal = group_detect.propose_grouping(
                aligned, ids=ids, fractions=fractions, stratify_col=args.stratify_col)
            _print_proposal(proposal)
            return 0
        if not args.group_col:
            if args.modality == "metadata":
                p.error("modality metadata needs a grouping column: run with --auto-detect-groups "
                        "to see the proposal, then re-run with --group-col <column>")
            # sequence + metadata but no group column and no detection requested: ignore metadata.
        else:
            keys = list(aligned.column(args.group_col))
            if args.drop_null_groups:
                keep = [i for i, k in enumerate(keys) if k is not None]
                ids = [ids[i] for i in keep]
                items = [items[i] for i in keep]
                keys = [keys[i] for i in keep]
            group_keys = keys
            group_provenance = {"group_col": args.group_col,
                                "id_col": args.id_col or group_detect.guess_id_column(aligned),
                                "confirmed_by": "cli_flag"}

    report = make_split(ids, items, modality=args.modality, fractions=fractions,
                        threshold=args.threshold, seed=args.seed, prefer_tool=args.prefer_tool,
                        group_keys=group_keys, group_provenance=group_provenance)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for split in set(report["assignment"]):
        sel = [ids[i] for i, s in enumerate(report["assignment"]) if s == split]
        (out / f"{split}.ids.txt").write_text("\n".join(sel) + "\n")
    (out / "report.json").write_text(json.dumps(report, indent=2))

    pr = report["provenance"]
    print(f"modality={pr['modality']} n={pr['n_items']} clusters={pr['n_clusters']} "
          f"backend={pr['backend'].get('backend')}")
    print(f"realized fractions: {report['realized_fractions']}")
    print(f"all invariants pass: {report['all_invariants_pass']}")
    if report["leakage"]:
        lk = report["leakage"]
        print(f"max train<->test similarity: {lk['max_cross_similarity']:.3f} "
              f"(threshold {lk['threshold']}, violations {lk['n_violations']})")
    if "grouping" in report:
        g = report["grouping"]
        print(f"grouping: col={g.get('group_col')} n_groups={g['n_groups']} "
              f"sizes(min/med/max)={g['group_size_distribution']['min']}/"
              f"{g['group_size_distribution']['median']:g}/{g['group_size_distribution']['max']} "
              f"composed_with_similarity={g['composed_with_similarity']} "
              f"clusters_merged_across_groups={g['clusters_merged_across_groups']}")
    if isinstance(pr["backend"], dict) and pr["backend"].get("warning"):
        print(f"WARNING: {pr['backend']['warning']}", file=sys.stderr)
    for w in report.get("warnings", []):
        print(f"WARNING: {w}", file=sys.stderr)
    print(f"wrote {out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
