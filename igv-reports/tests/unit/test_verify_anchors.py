"""Unit tests for verify_anchors.py — parser layer only.

Author: Samuel Ahuno
Purpose:
  Fast pytest suite covering the pure-Python parsing/decision logic in
  verify_anchors.py. No subprocess, no real BAM, no external dependencies.
  Runs in ~1 s on any machine with pytest.

  These tests catch the parser regressions that bit during the original
  iteration: status-taxonomy conflation between SKIP and FAIL, mis-tabbed
  TSV rows being silently mis-parsed, decode_status confusing tolerance
  with notes when columns are out of order.

Run:
  cd <plugin-root>/igv-reports
  pytest tests/unit/ -v
"""

from __future__ import annotations

import base64
import gzip
import json
import sys
from pathlib import Path

import pytest

# Make scripts/ importable without installing the skill as a package.
SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import verify_anchors as va  # noqa: E402


# ---------------------------------------------------------------------------
# load_anchors
# ---------------------------------------------------------------------------

def _write_tsv(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "anchors.tsv"
    p.write_text(body)
    return p


def test_load_anchors_full_row(tmp_path):
    """All 10 columns populated, including notes."""
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\ttumor\tchr2\t25246500\t25246501\t56\t0.05\t\t\tDNMT3A\n"
    ))
    rows = va.load_anchors(p)
    assert len(rows) == 1
    r = rows[0]
    assert (r.sample, r.track_name, r.chrom, r.start, r.end) == ("s1", "tumor", "chr2", 25246500, 25246501)
    assert r.expected == 56
    assert r.tolerance == "0.05"
    assert r.min_count == ""
    assert r.max_count == ""
    assert r.notes == "DNMT3A"


def test_load_anchors_min_max_row(tmp_path):
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\ttumor\tchrX\t100\t200\t50\t\t20\t100\thigh-conf\n"
    ))
    rows = va.load_anchors(p)
    assert rows[0].min_count == "20"
    assert rows[0].max_count == "100"


def test_load_anchors_missing_header_errors(tmp_path):
    """Data row before any header must abort with a clear error."""
    p = _write_tsv(tmp_path, "s1\ttumor\tchr1\t0\t100\t10\t\t\t\t\n")
    with pytest.raises(SystemExit, match="data row before header"):
        va.load_anchors(p)


def test_load_anchors_bad_tolerance_fails_fast(tmp_path):
    """Mis-tabbed row where notes value falls into tolerance must fail at
    load time with a hint, not crash later inside decide_status."""
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\ttumor\tchr2\t100\t200\t10\tNOT_A_NUMBER\t\t\tDNMT3A\n"
    ))
    with pytest.raises(SystemExit) as excinfo:
        va.load_anchors(p)
    msg = str(excinfo.value)
    assert "malformed anchor row" in msg
    assert "awk" in msg  # hint about -F'\t'


def test_load_anchors_bad_min_fails_fast(tmp_path):
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\ttumor\tchr2\t100\t200\t10\t\tNAH\t\t\n"
    ))
    with pytest.raises(SystemExit, match="malformed anchor row"):
        va.load_anchors(p)


def test_load_anchors_missing_file(tmp_path):
    with pytest.raises(SystemExit, match="anchors TSV not found"):
        va.load_anchors(tmp_path / "does_not_exist.tsv")


def test_load_anchors_skips_blank_lines(tmp_path):
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "\n"
        "s1\ttumor\tchr1\t0\t100\t10\t\t\t\t\n"
        "\n"
    ))
    rows = va.load_anchors(p)
    assert len(rows) == 1


# ---------------------------------------------------------------------------
# decide_status
# ---------------------------------------------------------------------------

def _anchor(expected=10, tolerance="", min_count="", max_count=""):
    return va.AnchorRow(
        sample="s", track_name="t", chrom="chr1", start=0, end=100,
        expected=expected, tolerance=tolerance,
        min_count=min_count, max_count=max_count,
    )


def test_decide_status_pass_within_default_tolerance():
    a = _anchor(expected=100)
    status, _ = va.decide_status(a, observed=104, default_tol=0.05)
    assert status == "PASS"


def test_decide_status_fail_outside_default_tolerance():
    a = _anchor(expected=100)
    status, details = va.decide_status(a, observed=110, default_tol=0.05)
    assert status == "FAIL"
    assert "diff_ratio" in details


def test_decide_status_per_row_tolerance_overrides_default():
    """Row tolerance 0.20 should pass observed=115 even though default 0.05 wouldn't."""
    a = _anchor(expected=100, tolerance="0.20")
    status, _ = va.decide_status(a, observed=115, default_tol=0.05)
    assert status == "PASS"


def test_decide_status_min_bound_pass():
    a = _anchor(expected=50, min_count="20")
    status, details = va.decide_status(a, observed=50, default_tol=0.05)
    assert status == "PASS"
    assert "min=20 OK" in details


def test_decide_status_min_bound_fail():
    a = _anchor(expected=50, min_count="100")
    status, details = va.decide_status(a, observed=50, default_tol=0.05)
    assert status == "FAIL"
    assert "min=100 FAIL" in details


def test_decide_status_min_max_combined():
    a = _anchor(min_count="20", max_count="80")
    status, _ = va.decide_status(a, observed=50, default_tol=0.05)
    assert status == "PASS"
    status, _ = va.decide_status(a, observed=10, default_tol=0.05)
    assert status == "FAIL"
    status, _ = va.decide_status(a, observed=100, default_tol=0.05)
    assert status == "FAIL"


def test_decide_status_bounds_override_tolerance():
    """When min/max present, tolerance is ignored."""
    # observed within tolerance of expected, but violates min
    a = _anchor(expected=50, tolerance="0.50", min_count="100")
    status, _ = va.decide_status(a, observed=52, default_tol=0.05)
    assert status == "FAIL"  # min wins over tolerance


def test_decide_status_zero_expected_exact():
    a = _anchor(expected=0)
    status, _ = va.decide_status(a, observed=0, default_tol=0.05)
    assert status == "PASS"
    status, _ = va.decide_status(a, observed=1, default_tol=0.05)
    assert status == "FAIL"


# ---------------------------------------------------------------------------
# decode_track_slice
# ---------------------------------------------------------------------------

def test_decode_track_slice_roundtrip(tmp_path):
    """data: URL → bytes round-trip preserves the payload."""
    payload = b"BAM\x01some bytes here"
    url = "data:application/gzip;base64," + base64.b64encode(payload).decode()
    dest = tmp_path / "out.bin"
    va.decode_track_slice(url, dest)
    assert dest.read_bytes() == payload


def test_decode_track_slice_other_mediatype_accepted(tmp_path):
    """We don't validate the mediatype — payload bytes are what matter."""
    payload = b"\x1f\x8b\x08compressed body"
    url = "data:application/octet-stream;base64," + base64.b64encode(payload).decode()
    dest = tmp_path / "out.bin"
    va.decode_track_slice(url, dest)
    assert dest.read_bytes() == payload


def test_decode_track_slice_not_a_data_url_raises(tmp_path):
    with pytest.raises(ValueError, match="not a data: base64 URL"):
        va.decode_track_slice("http://example.com/blob.bam", tmp_path / "out.bin")


# ---------------------------------------------------------------------------
# find_track
# ---------------------------------------------------------------------------

def test_find_track_hit():
    session = {"tracks": [
        {"name": "ann.bed"},
        {"name": "sample.sorted", "url": "data:..."},
    ]}
    t = va.find_track(session, "sample.sorted")
    assert t is not None and t["url"] == "data:..."


def test_find_track_miss():
    session = {"tracks": [{"name": "other"}]}
    assert va.find_track(session, "missing") is None


def test_find_track_empty():
    assert va.find_track({}, "x") is None
    assert va.find_track({"tracks": []}, "x") is None


# ---------------------------------------------------------------------------
# locate_session_entry — status taxonomy split (was the v1 regression)
# ---------------------------------------------------------------------------

def _make_table_json(rows):
    return {"headers": ["Chrom", "Start", "End", "Name"], "rows": rows}


def _make_session_dict(entries):
    """Build a sessionDictionary mapping str(idx) -> a gzipped+b64 data URL
    that decodes to the given entry dict."""
    out = {}
    for idx, entry in entries.items():
        raw = gzip.compress(json.dumps(entry).encode())
        out[str(idx)] = "data:application/gzip;base64," + base64.b64encode(raw).decode()
    return out


def test_locate_session_entry_ok():
    tj = _make_table_json([["chr2", 25246501, 25246501, "x"]])
    sd = _make_session_dict({0: {"tracks": [{"name": "t"}]}})
    outcome, sess, det = va.locate_session_entry(sd, tj, "chr2", 25246500, 25246501)
    assert outcome == "ok"
    assert sess == {"tracks": [{"name": "t"}]}
    assert det == ""


def test_locate_session_entry_absent_returns_skip_signal():
    """Anchor for a region that's not in the HTML — caller should SKIP."""
    tj = _make_table_json([["chr2", 25246501, 25246501, "x"]])
    sd = _make_session_dict({0: {"tracks": []}})
    outcome, _, det = va.locate_session_entry(sd, tj, "chr2", 99999999, 99999999)
    assert outcome == "absent"
    assert "no tableJson row matched" in det


def test_locate_session_entry_broken_missing_session():
    """Row in tableJson but no corresponding sessionDictionary entry — FAIL."""
    tj = _make_table_json([["chr2", 25246501, 25246501, "x"]])
    sd = {}  # no entries at all
    outcome, _, det = va.locate_session_entry(sd, tj, "chr2", 25246500, 25246501)
    assert outcome == "broken"
    assert "no entry for row index" in det


def test_locate_session_entry_broken_undecodable():
    """Row + session entry present but the session blob can't be gunzipped — FAIL."""
    tj = _make_table_json([["chr2", 25246501, 25246501, "x"]])
    sd = {"0": "data:application/gzip;base64,NOT_VALID_BASE64"}
    outcome, _, det = va.locate_session_entry(sd, tj, "chr2", 25246500, 25246501)
    assert outcome == "broken"
    assert "failed to gunzip/decode" in det


def test_locate_session_entry_broken_bad_headers():
    """tableJson missing the Chrom/Start/End columns we need."""
    tj = {"headers": ["foo", "bar"], "rows": [["x", "y"]]}
    sd = {}
    outcome, _, det = va.locate_session_entry(sd, tj, "chr2", 100, 200)
    assert outcome == "broken"
    assert "missing expected column" in det


# ---------------------------------------------------------------------------
# sample_bam_paths — samplesheet column handling
# ---------------------------------------------------------------------------

def test_sample_bam_paths_tumor_only():
    row = {"sample": "s1", "bam_tumor": "/x/tumor.sorted.bam"}
    out = va.sample_bam_paths(row)
    assert out == [("tumor.sorted", Path("/x/tumor.sorted.bam"))]


def test_sample_bam_paths_tumor_and_normal():
    row = {"sample": "s1", "bam_tumor": "/x/t.bam", "bam_normal": "/x/n.bam"}
    out = va.sample_bam_paths(row)
    names = [n for n, _ in out]
    assert names == ["t", "n"]


def test_sample_bam_paths_extras_filtered_to_bam_cram():
    row = {
        "sample": "s1",
        "bam_tumor": "/x/t.bam",
        "extra_tracks": "/y/extra.bam,/y/annot.bed,/y/other.cram",
    }
    out = va.sample_bam_paths(row)
    names = [n for n, _ in out]
    # bam_tumor + the .bam + the .cram from extras; .bed should be filtered out
    assert names == ["t", "extra", "other"]


def test_sample_bam_paths_blank_row():
    row = {"sample": "s1"}
    assert va.sample_bam_paths(row) == []


# ---------------------------------------------------------------------------
# write_anchors round-trip
# ---------------------------------------------------------------------------

def test_write_load_round_trip(tmp_path):
    anchors_in = [
        va.AnchorRow(sample="s1", track_name="t1", chrom="chr1",
                     start=0, end=100, expected=42, notes="hi"),
        va.AnchorRow(sample="s2", track_name="t2", chrom="chr2",
                     start=200, end=300, expected=7, min_count="3", max_count="20"),
    ]
    out = tmp_path / "anchors.tsv"
    va.write_anchors(anchors_in, out)
    rows = va.load_anchors(out)
    assert len(rows) == 2
    assert rows[0].notes == "hi"
    assert rows[1].min_count == "3"
    assert rows[1].max_count == "20"
