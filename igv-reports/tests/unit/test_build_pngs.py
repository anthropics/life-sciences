"""Unit tests for the --also-png plumbing in build_igvreports.py.

Author: Samuel Ahuno
Purpose:
  Exercises the helpers that bridge the HTML build to igver: sites BED
  parsing + UID assignment, flanked regions BED writer, input.txt writer,
  igver-cmd resolution, and manifest writing.

  We don't actually invoke igver here — the manifest writer reconstructs
  filenames from the same convention igver uses (validated against
  igver's _parse_bed_file source: `<chrom>-<start>-<end>.<name>.<ext>`).
  Cross-artifact consistency depends on this filename contract; if igver
  ever changes it, this test plus verify_cohort will catch the drift.

Run:
  pytest tests/unit/test_build_pngs.py -v
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import build_igvreports as b  # noqa: E402


def _write_bed(path: Path, rows: list[tuple]) -> None:
    """Helper — `rows` is a list of (chrom, start, end[, name][, ...])."""
    with path.open("w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


# ----- _read_sites_bed_rows -----


def test_read_sites_bed_assigns_uids_when_name_missing(tmp_path):
    bed = tmp_path / "sites.hg38.bed"
    _write_bed(bed, [("chr1", 100, 200), ("chr2", 300, 400)])
    rows = b._read_sites_bed_rows(bed)
    assert len(rows) == 2
    # Auto-UIDs are zero-padded to 3 digits so directory listings sort right
    # and `region_010` doesn't sort before `region_2`.
    assert rows[0]["name"] == "region_001"
    assert rows[1]["name"] == "region_002"
    assert rows[0]["bed_row_idx"] == 1
    assert rows[1]["bed_row_idx"] == 2


def test_read_sites_bed_preserves_existing_names(tmp_path):
    bed = tmp_path / "sites.hg38.bed"
    _write_bed(bed, [("chr2", 100, 200, "DNMT3A_full_gene"), ("chr7", 300, 400, "TP53")])
    rows = b._read_sites_bed_rows(bed)
    assert rows[0]["name"] == "DNMT3A_full_gene"
    assert rows[1]["name"] == "TP53"


def test_read_sites_bed_skips_comment_and_track_lines(tmp_path):
    bed = tmp_path / "sites.hg38.bed"
    bed.write_text(
        "#chrom\tstart\tend\tname\n"
        "track name=foo\n"
        "browser dense\n"
        "chr1\t100\t200\treal_row\n"
    )
    rows = b._read_sites_bed_rows(bed)
    assert len(rows) == 1
    assert rows[0]["name"] == "real_row"


def test_read_sites_bed_handles_mixed_named_and_unnamed(tmp_path):
    # If some rows have names and others don't, unnamed ones still get
    # auto-UIDs based on file position so manifests stay deterministic.
    bed = tmp_path / "sites.hg38.bed"
    _write_bed(bed, [
        ("chr1", 100, 200, "named_first"),
        ("chr2", 300, 400),
        ("chr3", 500, 600, "named_third"),
    ])
    rows = b._read_sites_bed_rows(bed)
    assert [r["name"] for r in rows] == ["named_first", "region_002", "named_third"]


# ----- _write_igver_regions_bed -----


def test_write_igver_regions_bed_applies_flanking(tmp_path):
    rows = [
        {"chrom": "chr1", "start": 100, "end": 200, "name": "A", "bed_row_idx": 1},
        {"chrom": "chr2", "start": 50,  "end": 150, "name": "B", "bed_row_idx": 2},
    ]
    out = tmp_path / "igver_regions.bed"
    b._write_igver_regions_bed(rows, flanking=300, out=out)
    lines = out.read_text().splitlines()
    # Row 1: 100-300=−200, clamped to 0; end 200+300=500.
    assert lines[0] == "chr1\t0\t500\tA"
    # Row 2: 50−300=−250, clamped to 0; end 150+300=450.
    assert lines[1] == "chr2\t0\t450\tB"


def test_write_igver_regions_bed_zero_flanking_passes_rows_verbatim(tmp_path):
    rows = [{"chrom": "chrX", "start": 1000, "end": 2000, "name": "promoter", "bed_row_idx": 1}]
    out = tmp_path / "igver_regions.bed"
    b._write_igver_regions_bed(rows, flanking=0, out=out)
    assert out.read_text().strip() == "chrX\t1000\t2000\tpromoter"


# ----- _write_igver_input_list -----


def test_write_igver_input_list_one_path_per_line(tmp_path):
    tracks = ["/path/to/tumor.bam", "/path/to/normal.bam", "/path/to/calls.vcf"]
    out = tmp_path / "igver_input.txt"
    b._write_igver_input_list(tracks, out)
    assert out.read_text().splitlines() == tracks


# ----- _resolve_igver_cmd -----


def test_resolve_igver_cmd_explicit_override_wins():
    override = "apptainer exec /path/to/igver.sif igver"
    assert b._resolve_igver_cmd(override) == override.split()


def test_resolve_igver_cmd_env_var_falls_back(monkeypatch):
    monkeypatch.setenv("IGVER_CMD", "/usr/local/bin/igver --debug")
    # which() must not find igver for this branch to fire; mock it to None.
    with patch.object(b.shutil, "which", return_value=None), \
         patch.object(b.Path, "exists", return_value=False):
        assert b._resolve_igver_cmd(None) == ["/usr/local/bin/igver", "--debug"]


def test_resolve_igver_cmd_path_lookup(monkeypatch):
    monkeypatch.delenv("IGVER_CMD", raising=False)
    with patch.object(b.shutil, "which", return_value="/usr/bin/igver"):
        assert b._resolve_igver_cmd(None) == ["/usr/bin/igver"]


def test_resolve_igver_cmd_raises_when_not_found(monkeypatch):
    monkeypatch.delenv("IGVER_CMD", raising=False)
    with patch.object(b.shutil, "which", return_value=None), \
         patch.object(b.Path, "exists", return_value=False):
        with pytest.raises(SystemExit, match="igver not found"):
            b._resolve_igver_cmd(None)


# ----- build_pngs_with_igver — mocked subprocess -----


def test_build_pngs_with_igver_writes_manifest_and_inputs(tmp_path, monkeypatch):
    # Set up a synthetic sites BED + tracks list + mock igver.
    bed = tmp_path / "sites.hg38.bed"
    _write_bed(bed, [
        ("chr1", 100, 200, "alpha"),
        ("chr2", 300, 400, "beta"),
    ])
    tracks = ["/data/sample.bam", "/data/calls.vcf"]
    html_path = tmp_path / "sample.hg38.html"
    html_path.write_text("<html/>")
    out_dir = tmp_path / "png_sample.hg38"
    log = logging.getLogger("test")

    # Force the resolver to return a stub so the test doesn't need igver.
    monkeypatch.setenv("IGVER_CMD", "/usr/bin/true")
    with patch.object(b.shutil, "which", return_value="/usr/bin/true"):
        manifest = b.build_pngs_with_igver(
            sites=bed,
            tracks=tracks,
            genome="hg38",
            flanking=300,
            out_dir=out_dir,
            log=log,
            html_path=html_path,
            igver_cmd=None,
            dpi=300,
            display_mode="collapse",
        )

    # 1. Intermediate files exist with the expected content.
    regions_bed = out_dir / "igver_regions.bed"
    input_txt = out_dir / "igver_input.txt"
    assert regions_bed.exists()
    assert input_txt.exists()
    assert regions_bed.read_text() == "chr1\t0\t500\talpha\nchr2\t0\t700\tbeta\n"
    assert input_txt.read_text() == "/data/sample.bam\n/data/calls.vcf\n"

    # 2. Manifest has one row per region with the right schema and the
    #    expected PNG-filename convention (validated against igver source).
    lines = manifest.read_text().splitlines()
    assert lines[0].startswith("#bed_row_idx\tuid\tchrom\t")
    data_rows = lines[1:]
    assert len(data_rows) == 2

    cols0 = data_rows[0].split("\t")
    assert cols0[0] == "1"
    assert cols0[1] == "alpha"
    assert cols0[2] == "chr1"
    assert cols0[3] == "100"            # start_orig
    assert cols0[4] == "200"            # end_orig
    assert cols0[5] == "0"              # start_flanked (clamped)
    assert cols0[6] == "500"            # end_flanked
    assert cols0[7] == "chr1:0-500"
    assert cols0[8].endswith("/png/chr1-0-500.alpha.png"), cols0[8]
    assert cols0[9].endswith("/sample.hg38.html"), cols0[9]
    assert cols0[10] == "1"             # html_table_row matches bed_row_idx


def test_build_pngs_with_igver_propagates_igver_failure(tmp_path, monkeypatch):
    # If igver itself returns non-zero, the driver must SystemExit so the
    # caller (and verify_cohort) sees the build as failed — silent success
    # would let an empty PNG dir slip into a "verified" cohort.
    bed = tmp_path / "sites.hg38.bed"
    _write_bed(bed, [("chr1", 100, 200, "alpha")])
    html_path = tmp_path / "sample.hg38.html"; html_path.write_text("<html/>")
    log = logging.getLogger("test")

    # /usr/bin/false always exits non-zero — perfect stand-in for a failing igver.
    monkeypatch.setenv("IGVER_CMD", "/usr/bin/false")
    with patch.object(b.shutil, "which", return_value="/usr/bin/false"):
        with pytest.raises(SystemExit) as exc:
            b.build_pngs_with_igver(
                sites=bed, tracks=["/data/sample.bam"], genome="hg38", flanking=0,
                out_dir=tmp_path / "out", log=log, html_path=html_path,
            )
    assert exc.value.code != 0


def test_build_pngs_with_igver_errors_on_empty_bed(tmp_path):
    bed = tmp_path / "sites.hg38.bed"
    bed.write_text("# header only\n")
    log = logging.getLogger("test")
    with pytest.raises(SystemExit, match="no data rows"):
        b.build_pngs_with_igver(
            sites=bed, tracks=["/data/sample.bam"], genome="hg38", flanking=0,
            out_dir=tmp_path / "out", log=log, html_path=tmp_path / "x.html",
        )
