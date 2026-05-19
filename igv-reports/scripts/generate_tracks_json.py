#!/usr/bin/env python3
"""generate_tracks_json.py — build an igv-reports tracks.json from a YAML spec.

Author: Samuel Ahuno
Purpose:
  ONT methylation viewers need named, colored, y-axis-locked tracks that
  the positional `create_report --tracks` API cannot express. The path is
  `--track-config <json>`, but hand-writing that JSON for 4-8 samples
  with 5mC + 5hmC bedGraph pairs each is tedious and error-prone.

  This helper consumes a small YAML spec (see
  examples/methylation_ont/tracks_spec.example.yaml) and emits the JSON
  with the right defaults baked in:

    * BAM tracks  -> colorBy=basemod2, showSoftClips=false, displayMode=COLLAPSED
    * bedGraph    -> type=wig, min=0, max=100 (methylation percent)
    * Annotation  -> displayMode honored, color honored
    * Group color -> reads from `group_colors:` map keyed by sample.group

Usage:
  python generate_tracks_json.py \
      --spec examples/methylation_ont/tracks_spec.example.yaml \
      --run-dir examples/methylation_ont \
      --out examples/methylation_ont/tracks.json

  --run-dir is prepended to any relative `url:` path in the spec, so the
  emitted JSON has absolute paths that create_report can resolve from any
  working directory.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not available. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(2)


BAM_DEFAULTS = {
    "format": "bam",
    "type": "alignment",
    "colorBy": "basemod2",
    "showSoftClips": False,
    "displayMode": "COLLAPSED",
}

BEDGRAPH_DEFAULTS = {
    "format": "bedgraph",
    "type": "wig",
    "min": 0,
    "max": 100,
}


# YAML shortcut keys (annotation: - default: <KEY>) map to the
# databases_config.yaml field for each genome plus display metadata.
# Colors are Okabe-Ito where chosen — colorblind-safe. format/displayMode
# match what build_igvreports.py emits on the non-track-config path.
ANNOTATION_DEFAULTS = {
    "cgi": {
        "display_name": "CpG islands",
        "yaml_key": "CpGIslands",
        "format": "bed",
        "displayMode": "EXPANDED",
        "color": "rgb(0,158,115)",       # Okabe-Ito green
    },
    "gencode": {
        "display_name": "Gencode",
        "yaml_key": "gtf",
        "format": "gff",                  # works for .gtf.gz and .gff3.gz
        "displayMode": "EXPANDED",
        "color": None,                    # IGV.js renders its own gene-track palette
    },
    "repmasker": {
        "display_name": "RepeatMasker",
        "yaml_key": "repMaskerBed",
        "format": "bed",
        "displayMode": "COLLAPSED",
        "color": None,
    },
    "epdnew_coding": {
        "display_name": "EPDnew (coding)",
        "yaml_key": "EPDnewCoding",
        "format": "bed",
        "displayMode": "EXPANDED",
        "color": "rgb(213,94,0)",         # Okabe-Ito vermillion
    },
    "epdnew_noncoding": {
        "display_name": "EPDnew (non-coding)",
        "yaml_key": "EPDnewNonCoding",
        "format": "bed",
        "displayMode": "EXPANDED",
        "color": "rgb(86,180,233)",       # Okabe-Ito sky blue
    },
}


def load_db_config(path: Path) -> dict:
    """Load databases_config.yaml; return {} on miss. Same semantics as the
    twin function in build_igvreports.py so the two stay aligned."""
    if not path.exists():
        sys.stderr.write(
            f"[generate_tracks_json] WARNING: db-config not found at {path}\n"
            "  Annotation entries using `default:` shortcuts will fail to resolve.\n"
            "  Use explicit `url:` paths, or set $IGV_REPORTS_DB_CONFIG.\n"
        )
        return {}
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def resolve_annotation_default(default_key: str, genome: str, cfg: dict) -> dict:
    """Look up a built-in annotation by short key (`cgi`, `gencode`, ...) for
    the given genome in the databases YAML. Returns a partial track dict with
    `display_name` / `url` / `indexURL` / `format` / `displayMode` / `color`
    populated; caller merges with name-overrides from the YAML.

    Raises SystemExit if the key is unknown, the genome is absent, or the
    resolved path doesn't exist on disk."""
    if default_key not in ANNOTATION_DEFAULTS:
        valid = ", ".join(sorted(ANNOTATION_DEFAULTS))
        raise SystemExit(
            f"ERROR: unknown annotation default '{default_key}'. Valid: {valid}"
        )
    meta = ANNOTATION_DEFAULTS[default_key]
    g = cfg.get("reference_genomes", {}).get("local", {}).get(genome, {})
    if not g:
        raise SystemExit(
            f"ERROR: db-config has no entry for genome '{genome}' "
            f"(needed to resolve `default: {default_key}`)."
        )
    yaml_key = meta["yaml_key"]
    raw = g.get(yaml_key)
    if not raw:
        raise SystemExit(
            f"ERROR: db-config has no '{yaml_key}' for genome '{genome}' "
            f"(needed to resolve `default: {default_key}`)."
        )
    # For hg38 gencode, prefer the bgzip+tabix .gff3.gz sibling if present
    # (mirrors build_igvreports.py:resolve_default_tracks gencode handling).
    url = raw
    if default_key == "gencode" and genome == "hg38":
        sibling = Path(raw).parent / "gencode.v47.annotation.gff3.gz"
        if sibling.exists() and (sibling.parent / (sibling.name + ".tbi")).exists():
            url = str(sibling)
    if not Path(url).exists():
        raise SystemExit(
            f"ERROR: resolved path missing on disk for `default: {default_key}` "
            f"({genome}): {url}"
        )
    # indexURL: include only if it actually exists. tabix .tbi is the standard
    # sibling for bgzipped tracks; igv.js falls back gracefully when absent.
    index_url = None
    for cand in (url + ".tbi", url + ".csi"):
        if Path(cand).exists():
            index_url = cand
            break

    track: dict = {
        "display_name": meta["display_name"],
        "url": url,
        "format": meta["format"],
        "displayMode": meta["displayMode"],
    }
    if index_url is not None:
        track["indexURL"] = index_url
    if meta["color"] is not None:
        track["color"] = meta["color"]
    return track


def abspath_relative_to(p: str, run_dir: Path) -> str:
    """Resolve `p` to an absolute path. If `p` is already absolute, return as-is."""
    pp = Path(p)
    if pp.is_absolute():
        return str(pp)
    return str((run_dir / pp).resolve())


def build_annotation_tracks(spec: dict, run_dir: Path, cfg: dict | None = None) -> list[dict]:
    """Build the annotation-track list. Each entry in `spec["annotation"]`
    is either:

      Explicit (existing behavior):
        - name: "Gencode v47"
          url: /abs/or/relative/path.gff3.gz
          indexURL: /abs/or/relative/path.gff3.gz.tbi  (optional)
          format: gff                                  (optional, default bed)
          displayMode: EXPANDED                        (optional)
          color: "rgb(...)"                            (optional)

      Shortcut (NEW — needs top-level `genome:` in spec and a loaded `cfg`):
        - default: gencode    # one of: cgi, gencode, repmasker,
                              #         epdnew_coding, epdnew_noncoding
          name: "Gencode v47"  # OPTIONAL override of the canned display name
          color: "rgb(...)"    # OPTIONAL override of the canned color
          displayMode: COLLAPSED  # OPTIONAL override

    Shortcut entries are resolved through resolve_annotation_default() against
    the databases YAML keyed by the spec's top-level `genome:`."""
    out: list[dict] = []
    genome = spec.get("genome")
    for a in spec.get("annotation", []):
        if "default" in a:
            if not genome:
                raise SystemExit(
                    "ERROR: annotation entry uses `default:` but spec is missing "
                    "top-level `genome:` — add e.g. `genome: hg38` to the YAML."
                )
            resolved = resolve_annotation_default(a["default"], genome, cfg or {})
            track = {
                "name": a.get("name", resolved["display_name"]),
                "url": resolved["url"],
                "format": a.get("format", resolved["format"]),
                "type": "annotation",
                "displayMode": a.get("displayMode", resolved["displayMode"]),
            }
            if "indexURL" in resolved:
                track["indexURL"] = resolved["indexURL"]
            if a.get("color") or resolved.get("color"):
                track["color"] = a.get("color", resolved.get("color"))
            out.append(track)
            continue
        # Explicit-path entry — preserves the prior behavior verbatim.
        track = {
            "name": a["name"],
            "url": abspath_relative_to(a["url"], run_dir),
            "format": a.get("format", "bed"),
            "type": "annotation",
            "displayMode": a.get("displayMode", "EXPANDED"),
        }
        if a.get("indexURL"):
            track["indexURL"] = abspath_relative_to(a["indexURL"], run_dir)
        if a.get("color"):
            track["color"] = a["color"]
        out.append(track)
    return out


def build_sample_tracks(spec: dict, run_dir: Path) -> list[dict]:
    group_colors = spec.get("group_colors", {})
    out: list[dict] = []
    for s in spec.get("samples", []):
        name = s["name"]
        group = s.get("group", "default")
        gc = group_colors.get(group, {})

        # BAM (per-read basemod2 view).
        if s.get("bam"):
            bam_abs = abspath_relative_to(s["bam"], run_dir)
            track = {"name": name, "url": bam_abs, "indexURL": bam_abs + ".bai"}
            track.update(BAM_DEFAULTS)
            out.append(track)

        # 5mC bedGraph.
        if s.get("bedgraph_5mC"):
            track = {
                "name": f"{name} 5mC",
                "url": abspath_relative_to(s["bedgraph_5mC"], run_dir),
            }
            track.update(BEDGRAPH_DEFAULTS)
            if gc.get("5mC"):
                track["color"] = gc["5mC"]
            out.append(track)

        # 5hmC bedGraph.
        if s.get("bedgraph_5hmC"):
            track = {
                "name": f"{name} 5hmC",
                "url": abspath_relative_to(s["bedgraph_5hmC"], run_dir),
            }
            track.update(BEDGRAPH_DEFAULTS)
            if gc.get("5hmC"):
                track["color"] = gc["5hmC"]
            out.append(track)

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", required=True, help="YAML spec (see tracks_spec.example.yaml)")
    ap.add_argument("--run-dir", required=True, help="dir that relative urls in spec are resolved against")
    ap.add_argument("--out", required=True, help="output tracks.json path")
    ap.add_argument("--db-config", default=os.environ.get("IGV_REPORTS_DB_CONFIG"), help=(
        "Databases YAML used to resolve `annotation: - default: <key>` shortcuts "
        "(cgi/gencode/repmasker/epdnew_coding/epdnew_noncoding) for the spec's "
        "`genome:`. Defaults to $IGV_REPORTS_DB_CONFIG. "
        "Not loaded if no shortcut entries appear. The YAML schema is "
        "`reference_genomes.local.<genome>.{CpGIslands,gtf,repMaskerBed,"
        "EPDnewCoding,EPDnewNonCoding}` — see references/databases_config_paths.md."
    ))
    ap.add_argument("--force", action="store_true",
                    help="overwrite --out if it already exists (default: refuse and exit 2 so hand-edits aren't clobbered)")
    args = ap.parse_args()

    spec_path = Path(args.spec)
    if not spec_path.exists():
        raise SystemExit(f"ERROR: spec not found: {spec_path}")
    run_dir = Path(args.run_dir).resolve()
    if not run_dir.exists():
        raise SystemExit(f"ERROR: run-dir not found: {run_dir}")

    with spec_path.open() as fh:
        spec = yaml.safe_load(fh)

    # Only load the db-config if any annotation entry uses the shortcut form;
    # specs that hand-paste paths remain self-contained.
    needs_cfg = any("default" in a for a in spec.get("annotation", []))
    if needs_cfg:
        if not args.db_config:
            raise SystemExit(
                "ERROR: spec has `default:` annotation shortcuts but --db-config "
                "is not set and $IGV_REPORTS_DB_CONFIG is empty.\n"
                "       Pass --db-config /path/to/databases.yaml, or convert the "
                "shortcuts to explicit `url:` entries."
            )
        cfg = load_db_config(Path(args.db_config))
    else:
        cfg = {}

    tracks = build_annotation_tracks(spec, run_dir, cfg) + build_sample_tracks(spec, run_dir)

    out_path = Path(args.out)
    if out_path.exists() and not args.force:
        raise SystemExit(
            f"ERROR: {out_path} already exists. A user may have hand-edited it after generation.\n"
            "       Pass --force to overwrite, or move the existing file aside and rerun."
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(tracks, fh, indent=2)
        fh.write("\n")

    print(f"Wrote {len(tracks)} tracks to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
