"""Unit tests for the PNG-side checks in verify_cohort.py.

Author: Samuel Ahuno
Purpose:
  When build_igvreports.py runs with --also-png, verify_cohort.py picks up
  the manifest TSV and runs three additional checks. These tests synthesize
  a valid manifest + matching PNG files in tmp_path, then mutate one thing
  at a time to confirm each check fires on the right defect.

Run:
  pytest tests/unit/test_verify_cohort_png.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import verify_cohort as vc  # noqa: E402


def _write_sites_bed(path: Path, rows: list[tuple]) -> None:
    with path.open("w") as fh:
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def _write_manifest(path: Path, entries: list[dict]) -> None:
    """Write a manifest TSV matching the schema build_pngs_with_igver emits."""
    header = ("#bed_row_idx\tuid\tchrom\tstart_orig\tend_orig\t"
              "start_flanked\tend_flanked\tregion\tpng_path\thtml_path\thtml_table_row\n")
    with path.open("w") as fh:
        fh.write(header)
        for e in entries:
            fh.write(
                f"{e['bed_row_idx']}\t{e['uid']}\t{e['chrom']}\t"
                f"{e['start_orig']}\t{e['end_orig']}\t"
                f"{e['start_flanked']}\t{e['end_flanked']}\t"
                f"{e['region']}\t{e['png_path']}\t{e['html_path']}\t"
                f"{e['html_table_row']}\n"
            )


def _make_png(path: Path, size_bytes: int = 50_000) -> None:
    """Create a fake PNG file of the requested size (default 50 KB, above
    the 10 KB threshold)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * (size_bytes - 8))


@pytest.fixture
def cohort(tmp_path):
    """Two-region cohort with a valid manifest + matching PNGs."""
    bed = tmp_path / "sites.hg38.bed"
    _write_sites_bed(bed, [
        ("chr1", 100, 200, "alpha"),
        ("chr2", 300, 400, "beta"),
    ])
    html = tmp_path / "sample.hg38.html"
    html.write_text("<html/>")

    png_dir = tmp_path / "png_sample.hg38" / "png"
    png1 = png_dir / "chr1-0-500.alpha.png"
    png2 = png_dir / "chr2-0-700.beta.png"
    _make_png(png1)
    _make_png(png2)

    manifest = tmp_path / "png_sample.hg38" / "manifest.tsv"
    _write_manifest(manifest, [
        {"bed_row_idx": 1, "uid": "alpha", "chrom": "chr1",
         "start_orig": 100, "end_orig": 200,
         "start_flanked": 0, "end_flanked": 500,
         "region": "chr1:0-500", "png_path": str(png1.resolve()),
         "html_path": str(html.resolve()), "html_table_row": 1},
        {"bed_row_idx": 2, "uid": "beta", "chrom": "chr2",
         "start_orig": 300, "end_orig": 400,
         "start_flanked": 0, "end_flanked": 700,
         "region": "chr2:0-700", "png_path": str(png2.resolve()),
         "html_path": str(html.resolve()), "html_table_row": 2},
    ])
    return {"bed": bed, "html": html, "manifest": manifest, "png_dir": png_dir,
            "png1": png1, "png2": png2, "tmp": tmp_path}


# ----- find_png_manifest -----


def test_find_png_manifest_returns_path_when_present(cohort, tmp_path):
    # Manifest lives at <reports_dir>/png_<sample>.<genome>/manifest.tsv
    # so we point reports_dir at tmp_path and check `sample` for genome `hg38`.
    found = vc.find_png_manifest(tmp_path, "sample", "hg38")
    assert found == cohort["manifest"]


def test_find_png_manifest_returns_none_when_absent(tmp_path):
    assert vc.find_png_manifest(tmp_path, "sample", "hg38") is None


# ----- P1: png count matches BED -----


def test_p1_pass(cohort):
    c = vc.check_png_count_matches_bed("sample", cohort["manifest"], cohort["bed"])
    assert c.status == "PASS"
    assert c.observed == "2"
    assert c.expected == "2"


def test_p1_fail_when_manifest_short(cohort):
    # Truncate the BED so it has 3 rows but manifest only has 2.
    _write_sites_bed(cohort["bed"], [
        ("chr1", 100, 200, "alpha"),
        ("chr2", 300, 400, "beta"),
        ("chr3", 500, 600, "gamma"),
    ])
    c = vc.check_png_count_matches_bed("sample", cohort["manifest"], cohort["bed"])
    assert c.status == "FAIL"
    assert c.observed == "2"
    assert c.expected == "3"


# ----- P2: pngs exist and non-empty -----


def test_p2_pass(cohort):
    c = vc.check_pngs_exist_and_nonempty("sample", cohort["manifest"])
    assert c.status == "PASS"


def test_p2_fail_on_missing_png(cohort):
    cohort["png1"].unlink()
    c = vc.check_pngs_exist_and_nonempty("sample", cohort["manifest"])
    assert c.status == "FAIL"
    assert "missing" in c.details


def test_p2_fail_on_tiny_png(cohort):
    # Re-write png1 as a 2 KB file — below the 10 KB threshold.
    cohort["png1"].write_bytes(b"\x00" * 2048)
    c = vc.check_pngs_exist_and_nonempty("sample", cohort["manifest"])
    assert c.status == "FAIL"
    assert "below threshold" in c.details


def test_p2_threshold_can_be_lowered(cohort):
    # The lab's smallest legitimate igver PNG can be ~5 KB on a no-data
    # region. Users should be able to opt down without rewriting the check.
    cohort["png1"].write_bytes(b"\x00" * 6144)
    c = vc.check_pngs_exist_and_nonempty("sample", cohort["manifest"], min_size_kb=5.0)
    assert c.status == "PASS"


# ----- P3: html-row alignment -----


def test_p3_pass(cohort):
    c = vc.check_png_html_row_alignment("sample", cohort["manifest"], cohort["html"])
    assert c.status == "PASS"


def test_p3_fail_when_html_path_diverges(cohort, tmp_path):
    # Pass a different HTML path than the manifest references — should fail.
    other_html = tmp_path / "other.hg38.html"
    other_html.write_text("<html/>")
    c = vc.check_png_html_row_alignment("sample", cohort["manifest"], other_html)
    assert c.status == "FAIL"
    assert "different HTML" in c.details


def test_p3_fail_when_row_indices_not_contiguous(cohort, tmp_path):
    # Rewrite the manifest with non-contiguous html_table_row indices.
    png1, png2 = cohort["png1"], cohort["png2"]
    _write_manifest(cohort["manifest"], [
        {"bed_row_idx": 1, "uid": "alpha", "chrom": "chr1",
         "start_orig": 100, "end_orig": 200, "start_flanked": 0, "end_flanked": 500,
         "region": "chr1:0-500", "png_path": str(png1.resolve()),
         "html_path": str(cohort["html"].resolve()), "html_table_row": 1},
        {"bed_row_idx": 2, "uid": "beta", "chrom": "chr2",
         "start_orig": 300, "end_orig": 400, "start_flanked": 0, "end_flanked": 700,
         "region": "chr2:0-700", "png_path": str(png2.resolve()),
         "html_path": str(cohort["html"].resolve()), "html_table_row": 5},  # gap
    ])
    c = vc.check_png_html_row_alignment("sample", cohort["manifest"], cohort["html"])
    assert c.status == "FAIL"
    assert "contiguous" in c.details
