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


def abspath_relative_to(p: str, run_dir: Path) -> str:
    """Resolve `p` to an absolute path. If `p` is already absolute, return as-is."""
    pp = Path(p)
    if pp.is_absolute():
        return str(pp)
    return str((run_dir / pp).resolve())


def build_annotation_tracks(spec: dict, run_dir: Path) -> list[dict]:
    out: list[dict] = []
    for a in spec.get("annotation", []):
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
    args = ap.parse_args()

    spec_path = Path(args.spec)
    if not spec_path.exists():
        raise SystemExit(f"ERROR: spec not found: {spec_path}")
    run_dir = Path(args.run_dir).resolve()
    if not run_dir.exists():
        raise SystemExit(f"ERROR: run-dir not found: {run_dir}")

    with spec_path.open() as fh:
        spec = yaml.safe_load(fh)

    tracks = build_annotation_tracks(spec, run_dir) + build_sample_tracks(spec, run_dir)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(tracks, fh, indent=2)
        fh.write("\n")

    print(f"Wrote {len(tracks)} tracks to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
