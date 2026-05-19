"""Unit tests for verify_report.py — pure-Python parser helpers.

Author: Samuel Ahuno
Purpose:
  Covers the HTML-extraction helpers and individual checks in
  verify_report.py without needing a real create_report HTML on disk:
  synthesized fixtures in tmp_path exercise every parser branch.

Run:
  cd claude/skills/igv-reports
  pytest tests/unit/ -v
"""

from __future__ import annotations

import base64
import gzip
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import verify_report as vr  # noqa: E402


# ---------------------------------------------------------------------------
# load_sites_bed
# ---------------------------------------------------------------------------

def test_load_sites_bed_basic(tmp_path):
    p = tmp_path / "sites.bed"
    p.write_text(
        "#chrom\tstart\tend\tname\n"
        "chr2\t25246500\t25246501\tDNMT3A\n"
        "chr7\t148884000\t148884001\tEZH2\n"
    )
    rows = vr.load_sites_bed(p)
    assert len(rows) == 2
    assert rows[0] == {"chrom": "chr2", "start": 25246500, "end": 25246501, "name": "DNMT3A"}
    assert rows[1]["name"] == "EZH2"


def test_load_sites_bed_skips_track_line(tmp_path):
    p = tmp_path / "sites.bed"
    p.write_text(
        'track name=foo description="bar"\n'
        "chr1\t100\t200\n"
    )
    rows = vr.load_sites_bed(p)
    assert len(rows) == 1
    assert rows[0]["name"] is None  # 3-col bed; no name


def test_load_sites_bed_blank_lines_ok(tmp_path):
    p = tmp_path / "sites.bed"
    p.write_text(
        "#header\n"
        "\n"
        "chr1\t100\t200\tx\n"
        "\n"
    )
    assert len(vr.load_sites_bed(p)) == 1


def test_load_sites_bed_too_few_cols(tmp_path):
    p = tmp_path / "sites.bed"
    p.write_text("chr1\t100\n")
    with pytest.raises(SystemExit, match="<3 columns"):
        vr.load_sites_bed(p)


def test_load_sites_bed_non_numeric(tmp_path):
    p = tmp_path / "sites.bed"
    p.write_text("chr1\tNOPE\t200\n")
    with pytest.raises(SystemExit, match="non-numeric"):
        vr.load_sites_bed(p)


# ---------------------------------------------------------------------------
# parse_table_json + parse_session_dictionary (regex extraction)
# ---------------------------------------------------------------------------

def _fake_html(table_json: dict, session_dict: dict) -> str:
    """Build a minimal HTML whose JS literals match what create_report emits."""
    return (
        "<html><body><script>\n"
        f"var tableJson = {json.dumps(table_json)};\n"
        f"var sessionDictionary = {json.dumps(session_dict)};\n"
        "</script></body></html>\n"
    )


def test_parse_table_json_extracts_dict():
    tj = {"headers": ["Chrom"], "rows": [["chr1"]]}
    html = _fake_html(tj, {})
    out = vr.parse_table_json(html)
    assert out == tj


def test_parse_session_dictionary_extracts_dict():
    sd = {"0": "data:application/gzip;base64,xxx"}
    html = _fake_html({"headers": [], "rows": []}, sd)
    out = vr.parse_session_dictionary(html)
    assert out == sd


def test_parse_table_json_missing_returns_none():
    assert vr.parse_table_json("<html>nothing here</html>") is None


def test_parse_balanced_blob_handles_braces_in_strings():
    """The brace-balancing scanner must not be tricked by '{' inside string literals."""
    html = "tableJson = {\"name\": \"value with { brace }\", \"n\": 1};"
    out = vr.parse_table_json(html)
    assert out["name"] == "value with { brace }"
    assert out["n"] == 1


def test_parse_balanced_blob_handles_escaped_quotes():
    """Backslash-escaped quotes must not flip the in_str state prematurely."""
    html = 'tableJson = {"name": "has \\" quote", "n": 2};'
    out = vr.parse_table_json(html)
    assert out["name"] == 'has " quote'
    assert out["n"] == 2


# ---------------------------------------------------------------------------
# decode_session_entry
# ---------------------------------------------------------------------------

def _make_data_url(payload: dict) -> str:
    raw = gzip.compress(json.dumps(payload).encode())
    return "data:application/gzip;base64," + base64.b64encode(raw).decode()


def test_decode_session_entry_roundtrip():
    payload = {"tracks": [{"name": "t", "url": "data:..."}]}
    url = _make_data_url(payload)
    assert vr.decode_session_entry(url) == payload


def test_decode_session_entry_bad_prefix_returns_none():
    assert vr.decode_session_entry("http://example.com/file.bam") is None


def test_decode_session_entry_bad_base64_returns_none():
    """Non-fatal — corrupted entries return None so caller can SKIP gracefully."""
    assert vr.decode_session_entry("data:application/gzip;base64,!!!notb64!!!") is None


# ---------------------------------------------------------------------------
# expected_track_labels (covers the Path.stem rule)
# ---------------------------------------------------------------------------

def test_expected_track_labels_from_paths():
    """Positional --tracks mode: igv-reports auto-names by Path.stem (strips
    ONE final suffix). Verified against create_report 1.16.2 in the script."""
    labs = vr.expected_track_labels(
        ["/x/sample.5mC.bedgraph", "/y/gencode.v47.annotation.gff3.gz", "/z/x.bam"],
        track_config=None,
    )
    assert labs == ["sample.5mC", "gencode.v47.annotation.gff3", "x"]


def test_expected_track_labels_from_track_config(tmp_path):
    """--track-config mode: use the `name` field from the JSON, not the path."""
    cfg = tmp_path / "tracks.json"
    cfg.write_text(json.dumps([
        {"name": "tumor", "url": "/x/tumor.bam"},
        {"name": "normal", "url": "/x/normal.bam"},
        {"url": "/x/no-name-track.bam"},  # entries without `name` are silently dropped
    ]))
    labs = vr.expected_track_labels([], track_config=cfg)
    assert labs == ["tumor", "normal"]


def test_expected_track_labels_empty():
    assert vr.expected_track_labels([], None) == []
    assert vr.expected_track_labels(None, None) == []


# ---------------------------------------------------------------------------
# Individual checks — drive them with synthetic inputs
# ---------------------------------------------------------------------------

def test_check_html_exists_pass(tmp_path):
    p = tmp_path / "r.html"; p.write_text("x")
    c = vr.check_html_exists(p)
    assert c.status == "PASS"


def test_check_html_exists_fail(tmp_path):
    c = vr.check_html_exists(tmp_path / "missing.html")
    assert c.status == "FAIL"


def test_check_html_min_size_pass_fail(tmp_path):
    p = tmp_path / "r.html"
    p.write_bytes(b"x" * (2 * 1024 * 1024))  # 2 MB
    assert vr.check_html_min_size(p, 1.0).status == "PASS"
    assert vr.check_html_min_size(p, 3.0).status == "FAIL"


def test_check_region_count_pass():
    bed = [{"chrom": "chr1", "start": 0, "end": 100, "name": "x"}]
    tj = {"headers": ["Chrom"], "rows": [["chr1"]]}
    assert vr.check_region_count(bed, tj).status == "PASS"


def test_check_region_count_fail_when_html_short():
    bed = [{"chrom": "chr1", "start": 0, "end": 100, "name": None}] * 3
    tj = {"headers": ["Chrom"], "rows": [["chr1"]]}
    c = vr.check_region_count(bed, tj)
    assert c.status == "FAIL"
    assert c.observed == "1"
    assert c.expected == "3"


def test_check_region_count_missing_table_json_fails():
    assert vr.check_region_count([{"chrom": "x", "start": 0, "end": 1, "name": None}], None).status == "FAIL"


def test_check_region_coords_match_with_name():
    """HTML stores 1-based start, BED is 0-based half-open."""
    bed = [{"chrom": "chr2", "start": 25246500, "end": 25246501, "name": "DNMT3A"}]
    tj = {
        "headers": ["Chrom", "Start", "End", "Name"],
        "rows": [["chr2", 25246501, 25246501, "DNMT3A"]],
    }
    assert vr.check_region_coords(bed, tj).status == "PASS"


def test_check_region_coords_name_mismatch_fails():
    bed = [{"chrom": "chr2", "start": 100, "end": 200, "name": "EXPECTED"}]
    tj = {
        "headers": ["Chrom", "Start", "End", "Name"],
        "rows": [["chr2", 101, 200, "DIFFERENT"]],
    }
    c = vr.check_region_coords(bed, tj)
    assert c.status == "FAIL"
    assert "name mismatch" in c.details


def test_check_region_coords_off_by_one_aware():
    """BED 0-based start 100 must match HTML 1-based start 101."""
    bed = [{"chrom": "chr1", "start": 100, "end": 200, "name": None}]
    tj = {
        "headers": ["Chrom", "Start", "End"],
        "rows": [["chr1", 101, 200]],
    }
    assert vr.check_region_coords(bed, tj).status == "PASS"


def test_check_region_sessions_pass():
    tj = {"headers": ["Chrom"], "rows": [["chr1"], ["chr2"]]}
    sd = {"0": "data:...", "1": "data:..."}
    assert vr.check_region_sessions(tj, sd).status == "PASS"


def test_check_region_sessions_missing_key():
    tj = {"headers": ["Chrom"], "rows": [["chr1"], ["chr2"]]}
    sd = {"0": "data:..."}  # missing "1"
    c = vr.check_region_sessions(tj, sd)
    assert c.status == "FAIL"


def test_check_tracks_present_pass():
    """Decode the first session entry and confirm all expected labels in tracks[].name."""
    sd = {"0": _make_data_url({"tracks": [{"name": "tumor"}, {"name": "cpg"}, {"name": "rmsk"}]})}
    c = vr.check_tracks_present(sd, ["tumor", "cpg"])
    assert c.status == "PASS"


def test_check_tracks_present_missing_track():
    sd = {"0": _make_data_url({"tracks": [{"name": "tumor"}]})}
    c = vr.check_tracks_present(sd, ["tumor", "missing_track"])
    assert c.status == "FAIL"
    assert "missing_track" in c.details


def test_check_tracks_present_skip_when_no_labels():
    sd = {"0": _make_data_url({"tracks": []})}
    assert vr.check_tracks_present(sd, []).status == "SKIP"


def test_check_tracks_present_empty_session_fails():
    c = vr.check_tracks_present({}, ["x"])
    assert c.status == "FAIL"


def test_check_tracks_present_undecodable_session_fails():
    c = vr.check_tracks_present({"0": "data:application/gzip;base64,!!!"}, ["x"])
    assert c.status == "FAIL"
    assert "gunzip" in c.details or "decode" in c.details
