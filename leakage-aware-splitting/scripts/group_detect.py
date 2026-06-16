#!/usr/bin/env python3
"""Group / metadata leakage detection for leakage-aware splitting.

Biological splits have two orthogonal leakage axes:

1. **Similarity** - near-duplicate items (homologous sequences, shared scaffolds)
   handled by the per-modality clustering in cluster_*.py / split_*.py.
2. **Group / metadata** - replicate units that share a *source entity* must travel
   to the same split. The canonical footgun is histology: tiles come from slides,
   slides from patients, patients from sites/scanners. If you split by tile, the
   same patient lands on both sides and the model is scored on memorizing patients,
   not pathology. The same shape covers subject/donor/case (clinical),
   batch/plate/well/site/scanner (assay batch effects), sample/library/replicate
   (sequencing), study/source/dataset (pooled data), and organism/strain.

This module profiles a dataset's metadata, proposes a grouping key anchored on the
question *"what is the unit of generalization?"*, and emits the intra-group edges
that the shared union-find core (`core.clusters_from_pairs`) turns into whole-group
clusters. A group key is just another edge source, so group-aware splitting composes
with similarity clustering in a single union-find pass and inherits every existing
guarantee (whole-component assignment, cluster integrity, reproducibility).

**Detection is propose-then-confirm, never silent.** Automatic group-column
detection is not a solved problem - no off-the-shelf tool does it reliably - so
`propose_grouping` returns a *ranked proposal and applies nothing*. A human (or the
agent, after asking) picks the unit of generalization. The CLI enforces this: it
prints the ranked proposal and exits unless an explicit `--group-col` is given.

Loading uses the stdlib `csv` module by default (no pandas dependency), mirroring
the `mmseqs_available()`-with-fallback idiom elsewhere in the skill, with optional
pandas when importable and explicitly opted into.
"""

from __future__ import annotations

import csv
import difflib
import re
import statistics
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

# --------------------------------------------------------------------------- #
# Grouping-key vocabulary                                                      #
# --------------------------------------------------------------------------- #
# Families of column-name keywords that name a grouping unit. Weighted toward
# families that are genuine *units of generalization* (a model should generalize to
# a new patient / site / study) over merely technical sub-units (a new well, a new
# image tile). The weight feeds the proposal score; the hierarchy bonus (coarseness)
# is added on top, so a coarse, generalization-unit column wins on both counts.
GROUPING_KEY_DICTIONARY: dict[str, list[str]] = {
    "subject": ["patient", "subject", "donor", "case", "participant", "individual",
                "animal", "mouse", "rat", "person"],
    "study": ["study", "source", "dataset", "cohort", "collection", "project"],
    "site": ["site", "center", "centre", "scanner", "device", "institution",
             "hospital", "clinic", "lab", "machine", "instrument"],
    "organism": ["organism", "strain", "species", "genotype", "isolate"],
    "sample": ["sample", "specimen", "library", "replicate", "aliquot", "biopsy", "subjectsample"],
    "batch": ["batch", "plate", "well", "run", "lane", "flowcell", "experiment", "assay"],
    "image": ["slide", "wsi", "image", "tile", "frame", "patch", "field", "roi"],
}

FAMILY_WEIGHTS: dict[str, float] = {
    "subject": 1.0,    # patient / donor - the canonical generalization unit
    "study": 0.9,      # source / dataset / cohort - pooled-data leakage
    "site": 0.85,      # scanner / center - batch / site effects
    "organism": 0.8,
    "sample": 0.6,     # specimen / replicate
    "batch": 0.55,     # plate / well / run
    "image": 0.4,      # slide / tile / frame - rarely the generalization unit alone
}

# Column names that look like a prediction target, not a grouping unit. Grouping by
# the label would force all of one class into a single split (catastrophic).
LABEL_HINTS = {"label", "labels", "class", "classes", "target", "targets", "y",
               "outcome", "response", "activity", "active", "inactive", "category",
               "pic50", "pki", "ic50", "ec50", "affinity", "binding", "value", "score"}

# Null tokens normalized to None on load (case-insensitive, stripped).
NULL_TOKENS = {"", "na", "nan", "none", "null", "n/a", "."}

# role_hint thresholds (uniqueness_ratio = n_unique / n_rows).
_GROUPING_UNIQUENESS_MAX = 0.95   # above this a column is too unique to group usefully
_CONTINUOUS_UNIQUENESS_MIN = 0.5  # numeric AND this unique -> looks like a measurement
_FREEFORM_UNIQUENESS_MIN = 0.9    # near-unique non-numeric text
_LABEL_MAX_CARDINALITY = 20       # label-ish columns have few distinct values
_FUZZY_RATIO = 0.85               # difflib ratio threshold for a fuzzy name match


# --------------------------------------------------------------------------- #
# Metadata table (column-major, stdlib-only)                                   #
# --------------------------------------------------------------------------- #

@dataclass
class MetadataTable:
    """Column-major metadata: columns in file order, data[col] aligned by row.

    Column-major keeps each column a flat list, which is exactly what profiling and
    edge-emission consume, and lets ``align_to`` permute rows to item order before
    any column becomes a per-item array (so the orchestrator stays array-in/array-out).
    """

    columns: list[str]
    data: dict[str, list]  # values are str or None
    n_rows: int

    def column(self, name: str) -> list:
        if name not in self.data:
            raise KeyError(f"column '{name}' not in metadata (have: {self.columns})")
        return self.data[name]

    def align_to(self, ids: Sequence[str], id_col: str) -> "MetadataTable":
        """Return a copy whose rows are reordered/subset to match ``ids``.

        Uses ``id_col`` as the join key. Errors on a duplicate id in the table (a
        duplicate id is itself a leakage finding - two rows claiming the same item)
        and on any requested id missing from the table. Alignment happens before
        columns are handed to the splitter, so downstream code never has to reason
        about row order.
        """
        key = self.column(id_col)
        pos: dict[str, int] = {}
        dups: list[str] = []
        for i, v in enumerate(key):
            sv = None if v is None else str(v)
            if sv in pos:
                dups.append(sv)
            else:
                pos[sv] = i
        if dups:
            sample = ", ".join(sorted(set(dups))[:5])
            raise ValueError(
                f"id column '{id_col}' has {len(set(dups))} duplicated id(s) (e.g. {sample}); "
                "duplicate ids are themselves a finding - resolve before splitting"
            )
        missing = [str(x) for x in ids if str(x) not in pos]
        if missing:
            sample = ", ".join(missing[:5])
            raise ValueError(
                f"{len(missing)} item id(s) not found in metadata id column '{id_col}' "
                f"(e.g. {sample})"
            )
        order = [pos[str(x)] for x in ids]
        new_data = {c: [self.data[c][i] for i in order] for c in self.columns}
        return MetadataTable(columns=list(self.columns), data=new_data, n_rows=len(order))


def _normalize_cell(value: str | None) -> str | None:
    if value is None:
        return None
    v = value.strip()
    if v.lower() in NULL_TOKENS:
        return None
    return v


def _sniff_delimiter(path: Path, sample: str) -> str:
    ext = path.suffix.lower()
    if ext in (".tsv", ".tab"):
        return "\t"
    if ext == ".csv":
        return ","
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter
    except csv.Error:
        return ","


def load_metadata(path: str, *, prefer_pandas: bool = False) -> MetadataTable:
    """Load a CSV/TSV metadata file into a MetadataTable.

    Stdlib ``csv`` by default (delimiter from extension or ``csv.Sniffer``; the null
    tokens in ``NULL_TOKENS`` become ``None``). When ``prefer_pandas`` is set and
    pandas is importable, reads via pandas and coerces to the same MetadataTable, so
    behaviour is identical to callers regardless of backend.
    """
    p = Path(path)
    if prefer_pandas:
        try:
            import pandas as pd  # type: ignore
        except ImportError:
            pass
        else:
            sep = _sniff_delimiter(p, p.read_text()[:4096])
            df = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)
            cols = [str(c) for c in df.columns]
            data = {c: [_normalize_cell(x) for x in df[c].tolist()] for c in cols}
            return MetadataTable(columns=cols, data=data, n_rows=len(df))

    text = p.read_text()
    delimiter = _sniff_delimiter(p, text[:4096])
    reader = csv.reader(text.splitlines(), delimiter=delimiter)
    rows = list(reader)
    if not rows:
        return MetadataTable(columns=[], data={}, n_rows=0)
    header = [h.strip() for h in rows[0]]
    data: dict[str, list] = {c: [] for c in header}
    for row in rows[1:]:
        for j, col in enumerate(header):
            data[col].append(_normalize_cell(row[j]) if j < len(row) else None)
    return MetadataTable(columns=header, data=data, n_rows=len(rows) - 1)


# --------------------------------------------------------------------------- #
# Column profiling                                                             #
# --------------------------------------------------------------------------- #

@dataclass
class ColumnProfile:
    name: str
    n_unique: int
    n_null: int
    uniqueness_ratio: float
    largest_group_frac: float
    name_match: str | None          # matched GROUPING_KEY_DICTIONARY family, or None
    match_kind: str                 # exact | substring | fuzzy | none
    examples: list[str]
    role_hint: str                  # grouping_candidate|id_like|constant|label_like|continuous_like|freeform


def _match_family(colname: str) -> tuple[str | None, str, float]:
    """Best (family, kind, score) match of a column name to the grouping dictionary."""
    norm = re.sub(r"[^a-z0-9]+", " ", colname.lower()).strip()
    flat = norm.replace(" ", "")
    tokens = norm.split()
    best: tuple[str | None, str, float] = (None, "none", 0.0)
    for family, kws in GROUPING_KEY_DICTIONARY.items():
        for kw in kws:
            if kw in tokens:
                return family, "exact", 1.0  # exact token match cannot be beaten
            if kw in flat and len(kw) >= 4:
                if 0.8 > best[2]:
                    best = (family, "substring", 0.8)
            else:
                for tok in tokens:
                    r = difflib.SequenceMatcher(None, tok, kw).ratio()
                    if r >= _FUZZY_RATIO and r > best[2]:
                        best = (family, "fuzzy", r)
    return best


def _fraction_numeric(values: list[str]) -> float:
    if not values:
        return 0.0
    n = 0
    for v in values:
        try:
            float(v)
            n += 1
        except (TypeError, ValueError):
            pass
    return n / len(values)


def _looks_like_label(colname: str) -> bool:
    norm = re.sub(r"[^a-z0-9]+", " ", colname.lower()).strip()
    return any(tok in LABEL_HINTS for tok in norm.split())


def _infer_role(name: str, n_rows: int, n_unique: int, nonnull: list[str],
                name_match: str | None) -> str:
    if n_unique <= 1:
        return "constant"
    uniqueness = n_unique / n_rows if n_rows else 0.0
    if uniqueness >= 1.0:
        return "id_like"
    numeric = _fraction_numeric(nonnull) >= 0.95
    groups_repeat = uniqueness < _GROUPING_UNIQUENESS_MAX
    if _looks_like_label(name) and n_unique <= _LABEL_MAX_CARDINALITY:
        return "label_like"
    if name_match is not None and groups_repeat:
        return "grouping_candidate"
    if numeric and uniqueness >= _CONTINUOUS_UNIQUENESS_MIN:
        return "continuous_like"
    if uniqueness >= _FREEFORM_UNIQUENESS_MIN:
        return "freeform"
    if groups_repeat:
        return "grouping_candidate"
    return "freeform"


def profile_columns(table: MetadataTable) -> list[ColumnProfile]:
    """Per-column statistics + a role_hint that gates false alarms."""
    profiles = []
    for col in table.columns:
        values = table.data[col]
        nonnull = [v for v in values if v is not None]
        n_null = len(values) - len(nonnull)
        counts: dict[str, int] = defaultdict(int)
        examples: list[str] = []
        for v in nonnull:
            counts[v] += 1
            if v not in examples and len(examples) < 5:
                examples.append(v)
        n_unique = len(counts)
        uniqueness = n_unique / table.n_rows if table.n_rows else 0.0
        largest = max(counts.values()) / table.n_rows if counts and table.n_rows else 0.0
        family, kind, _score = _match_family(col)
        role = _infer_role(col, table.n_rows, n_unique, nonnull, family)
        profiles.append(ColumnProfile(
            name=col, n_unique=n_unique, n_null=n_null,
            uniqueness_ratio=uniqueness, largest_group_frac=largest,
            name_match=family, match_kind=kind, examples=examples, role_hint=role,
        ))
    return profiles


# --------------------------------------------------------------------------- #
# Functional dependencies / hierarchy                                          #
# --------------------------------------------------------------------------- #

@dataclass
class Hierarchy:
    fd_edges: set                   # (A, B): A functionally determines B
    strict_finer: set               # (A, B): A strictly finer than B (A->B, not B->A)
    coarseness: dict                # col -> #columns it is strictly coarser than
    chains: list                    # human-readable finest->coarsest orderings
    excluded_null_pairs: list       # (A, B, n_excluded) pairs with rows dropped for nulls


def _determines(a_vals: list, b_vals: list) -> tuple[bool, int]:
    """Does A determine B? Returns (holds, n_rows_excluded_for_nulls)."""
    seen: dict = {}
    excluded = 0
    for x, y in zip(a_vals, b_vals):
        if x is None or y is None:
            excluded += 1
            continue
        if x in seen:
            if seen[x] != y:
                return False, excluded
        else:
            seen[x] = y
    return True, excluded


def functional_dependencies(table: MetadataTable, cols: Sequence[str]) -> set:
    """Set of (A, B) where every value of A co-occurs with exactly one value of B.

    O(k^2 * n_rows) for k columns (k is small - only surviving candidates). Rows null
    in either column are excluded from that pair's test.
    """
    edges = set()
    for a in cols:
        for b in cols:
            if a == b:
                continue
            holds, _ = _determines(table.data[a], table.data[b])
            if holds:
                edges.add((a, b))
    return edges


def infer_hierarchy(table: MetadataTable, cols: Sequence[str]) -> Hierarchy:
    """Reduce functional dependencies among ``cols`` to a finest->coarsest picture.

    tile->slide->patient: each tile is on one slide, each slide on one patient, so
    tile determines slide determines patient. The *coarsest* key (patient) gives the
    most protection, so it gets the largest coarseness score; the whole chain is
    surfaced so the user can pick the unit of generalization.
    """
    fd = functional_dependencies(table, cols)
    excluded = []
    for (a, b) in fd:
        _, n_ex = _determines(table.data[a], table.data[b])
        if n_ex:
            excluded.append((a, b, n_ex))
    strict = {(a, b) for (a, b) in fd if (b, a) not in fd}
    coarseness = {c: 0 for c in cols}
    for (a, b) in strict:
        # a is strictly finer than b -> b is coarser; credit b.
        coarseness[b] += 1
    # Build human-readable chains by ordering columns within a connected strict DAG
    # by coarseness (finest first). Robust to non-total orders (just lists by rank).
    related = sorted({c for e in strict for c in e}, key=lambda c: (coarseness[c], c))
    chains = []
    if related:
        chains.append(related)  # finest -> coarsest
    return Hierarchy(fd_edges=fd, strict_finer=strict, coarseness=coarseness,
                     chains=chains, excluded_null_pairs=excluded)


# --------------------------------------------------------------------------- #
# Proposal (pure - ranks, applies nothing)                                     #
# --------------------------------------------------------------------------- #

@dataclass
class GroupingCandidate:
    column: str
    canonical_key: str
    n_groups: int
    group_size_min: int
    group_size_max: int
    group_size_median: float
    uniqueness_ratio: float
    satisfiable_for: bool
    rank_score: float
    rationale: dict
    warnings: list = field(default_factory=list)
    disqualified: bool = False


@dataclass
class GroupingProposal:
    candidates: list           # ranked best-first; disqualified sink to the bottom
    recommended: GroupingCandidate | None
    tied: list                 # >1 entry => ambiguous, the agent MUST ask
    profiles: list
    hierarchy: Hierarchy
    notes: list = field(default_factory=list)


def _group_sizes(values: Sequence) -> dict:
    sizes: dict = defaultdict(int)
    for v in values:
        if v is not None:
            sizes[v] += 1
    return sizes


def _satisfiability(n_groups: int, largest_frac: float, fractions: dict, tol: float = 0.15) -> tuple[bool, str]:
    """Can the requested fractions plausibly be hit with whole groups? (ok, reason).

    Three ways it fails: fewer groups than splits (a split must be empty); one group
    bigger than the largest split target (overshoots no matter where it goes - the
    'clean split impossible' case); or so few groups that the smallest split expects
    less than one whole group (it is liable to come out empty). Each is the kind of
    thing that otherwise silently produces a misleading, e.g. empty, test set.
    """
    nz = [v for v in fractions.values() if v > 0]
    n_splits = len(nz)
    if n_groups < n_splits:
        return False, f"only {n_groups} group(s) for {n_splits} splits - some split would be empty"
    if largest_frac > max(fractions.values()) + tol:
        return False, (f"one group holds {largest_frac:.0%} of rows, over the largest split "
                       "target - a clean split at this ratio is impossible")
    if min(nz) * n_groups < 0.75:
        return False, (f"{n_groups} groups is too few for the smallest requested split "
                       f"(fraction {min(nz):g} expects <1 group) - a split may come out empty; "
                       "use a coarser ratio or gather more groups")
    return True, ""


def propose_grouping(
    table: MetadataTable,
    *,
    ids: Sequence[str] | None = None,
    fractions: dict | None = None,
    stratify_col: str | None = None,
    hints: Sequence[str] | None = None,
) -> GroupingProposal:
    """Rank grouping-key candidates. PURE: returns a proposal, applies nothing.

    Score = dictionary-family weight + hierarchy coarseness + satisfiability. Hard
    disqualifications (emitted as warnings, never silently dropped): unique-per-row
    (grouping is a no-op), constant (one group = no split possible), and any column
    equal to ``stratify_col`` (grouping by the label leaks the target). Ties at the
    top are surfaced in ``tied`` so the caller asks rather than guessing.
    """
    if fractions is None:
        fractions = {"train": 0.8, "val": 0.1, "test": 0.1}
    profiles = profile_columns(table)
    by_name = {p.name: p for p in profiles}

    # Hierarchy is computed over columns that could plausibly be groups.
    candidate_cols = [p.name for p in profiles
                      if p.role_hint in ("grouping_candidate", "id_like")]
    hierarchy = infer_hierarchy(table, candidate_cols)

    hint_set = {h.lower() for h in (hints or [])}
    candidates: list[GroupingCandidate] = []
    for p in profiles:
        warnings: list[str] = []
        disqualified = False
        if p.role_hint == "constant" or p.n_unique <= 1:
            disqualified = True
            warnings.append("constant column (a single group cannot be split)")
        if p.uniqueness_ratio >= 1.0:
            disqualified = True
            warnings.append("unique per row - grouping by this column is a no-op (it is the row id)")
        if p.largest_group_frac >= 0.9 and p.uniqueness_ratio < 1.0 and p.n_unique >= 2:
            disqualified = True
            warnings.append(f"near-constant: one group holds {p.largest_group_frac:.0%} of rows; "
                            "grouping on it forces ~all data into a single split (refused - it "
                            "would emit a misleading split rather than a useful one)")
        if stratify_col is not None and p.name == stratify_col:
            disqualified = True
            warnings.append(f"equals the stratification/label column '{stratify_col}' - "
                            "grouping by the label would leak the target")
        if p.role_hint in ("label_like",) and not disqualified:
            warnings.append("name looks like a prediction target, not a grouping unit")

        # Only score plausible group columns; skip pure freeform/continuous noise
        # unless the name matched the dictionary or the user hinted at it.
        plausible = (p.role_hint == "grouping_candidate"
                     or p.name_match is not None
                     or p.name.lower() in hint_set)
        if not plausible and not disqualified:
            continue

        sizes = _group_sizes(table.data[p.name])
        size_vals = sorted(sizes.values())
        n_groups = len(size_vals)
        gmin = size_vals[0] if size_vals else 0
        gmax = size_vals[-1] if size_vals else 0
        gmed = statistics.median(size_vals) if size_vals else 0.0
        sat_ok, sat_reason = _satisfiability(n_groups, p.largest_group_frac, fractions)
        satisfiable = (not disqualified) and sat_ok
        if not disqualified and not sat_ok:
            warnings.append(sat_reason)

        family_weight = FAMILY_WEIGHTS.get(p.name_match, 0.0) if p.name_match else 0.0
        coarseness = hierarchy.coarseness.get(p.name, 0)
        hint_bonus = 0.5 if p.name.lower() in hint_set else 0.0
        score = family_weight + 0.5 * coarseness + (0.5 if satisfiable else 0.0) + hint_bonus
        if disqualified:
            score = float("-inf")

        candidates.append(GroupingCandidate(
            column=p.name,
            canonical_key=p.name_match or p.name,
            n_groups=n_groups,
            group_size_min=gmin, group_size_max=gmax, group_size_median=gmed,
            uniqueness_ratio=p.uniqueness_ratio,
            satisfiable_for=satisfiable,
            rank_score=score,
            rationale={
                "family": p.name_match, "family_weight": family_weight,
                "match_kind": p.match_kind, "coarseness": coarseness,
                "satisfiable_for": satisfiable, "hint_bonus": hint_bonus,
                "role_hint": p.role_hint,
            },
            warnings=warnings, disqualified=disqualified,
        ))

    # Rank: non-disqualified by score desc (ties broken by column name); disqualified last.
    candidates.sort(key=lambda c: (c.disqualified, -c.rank_score, c.column))

    live = [c for c in candidates if not c.disqualified]
    recommended = live[0] if live else None
    tied = []
    if len(live) >= 2 and abs(live[0].rank_score - live[1].rank_score) < 1e-9:
        tied = [c for c in live if abs(c.rank_score - live[0].rank_score) < 1e-9]

    notes: list[str] = []
    if recommended is None:
        notes.append("No viable grouping column found - every candidate was disqualified "
                     "(unique-per-row, constant, or the label). Confirm there is a real "
                     "grouping unit, or split on similarity only.")
    if tied:
        cols = ", ".join(c.column for c in tied)
        notes.append(f"Ambiguous: {cols} tie for best. Ask the user which is the unit of "
                     "generalization before splitting.")
    if hierarchy.strict_finer:
        nest = "; ".join(f"{a} ⊂ {b}" for (a, b) in sorted(hierarchy.strict_finer))
        notes.append(f"Nested columns detected ({nest}). The coarsest gives the most "
                     "protection; the chain is surfaced so you can pick the right level.")
    if hierarchy.excluded_null_pairs:
        notes.append("Some functional-dependency tests excluded rows with null keys; "
                     "see hierarchy.excluded_null_pairs.")

    return GroupingProposal(candidates=candidates, recommended=recommended, tied=tied,
                            profiles=profiles, hierarchy=hierarchy, notes=notes)


# --------------------------------------------------------------------------- #
# Edge emission (the bridge to core.clusters_from_pairs)                       #
# --------------------------------------------------------------------------- #

def _spanning_pairs(keys: Sequence, *, skip_null: bool) -> list:
    """Spanning chain (n-1 edges) per distinct key. NOT the O(n^2) clique.

    Connecting each member to the previous member of the same key is enough for
    union-find to merge the whole group; it costs one edge per item instead of one
    per pair. With ``skip_null``, None keys emit no edges (they become singletons).
    """
    last: dict = {}
    pairs: list = []
    for i, k in enumerate(keys):
        if skip_null and k is None:
            continue
        if k in last:
            pairs.append((last[k], i))
        last[k] = i
    return pairs


def pairs_from_groups(group_keys: Sequence) -> list:
    """Intra-group edges for union-find. Null keys -> singletons (no edge)."""
    return _spanning_pairs(group_keys, skip_null=True)


def pairs_from_labels(labels: Sequence[int]) -> list:
    """Lift cluster *labels* back to spanning edges.

    A label-producing backend (e.g. MMseqs2) returns per-item cluster ids, not pairs.
    To re-merge those clusters with group edges in a single union-find pass, turn the
    labels back into spanning edges. Labels are never null.
    """
    return _spanning_pairs(labels, skip_null=False)


def guess_id_column(table: MetadataTable) -> str | None:
    """Best-effort id column: a unique-per-row column, preferring id-named ones.

    Deterministic. Returns None if no column is unique per row (caller should then
    require an explicit --id-col).
    """
    profiles = {p.name: p for p in profile_columns(table)}
    unique_cols = [c for c in table.columns if profiles[c].uniqueness_ratio >= 1.0]
    if not unique_cols:
        return None
    id_named = [c for c in unique_cols
                if re.search(r"(^|[_\W])id($|[_\W])", c.lower()) or c.lower() == "id"]
    return (id_named or unique_cols)[0]


def group_size_distribution(group_keys: Sequence) -> dict:
    """Summary of group sizes for the report (min/max/median/p95/singletons/nulls)."""
    sizes = sorted(_group_sizes(group_keys).values())
    n_null = sum(1 for k in group_keys if k is None)
    if not sizes:
        return {"n_groups": 0, "min": 0, "max": 0, "median": 0.0, "p95": 0,
                "n_singletons": 0, "n_null_group_keys": n_null}
    p95_idx = max(0, min(len(sizes) - 1, int(round(0.95 * (len(sizes) - 1)))))
    return {
        "n_groups": len(sizes),
        "min": sizes[0], "max": sizes[-1],
        "median": statistics.median(sizes),
        "p95": sizes[p95_idx],
        "n_singletons": sum(1 for s in sizes if s == 1),
        "n_null_group_keys": n_null,
    }
