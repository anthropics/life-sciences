"""Unit tests for verify_anchors.py — parser layer only.

Author: Samuel Ahuno
Purpose:
  Fast pytest suite covering the pure-Python parsing/decision logic in
  verify_anchors.py. No subprocess, no real BAM, no shared-storage
  dependency. Runs in ~1 s on any machine with pytest.

  These tests catch the parser regressions that bit during the original
  iteration: status-taxonomy conflation between SKIP and FAIL, mis-tabbed
  TSV rows being silently mis-parsed, decode_status confusing tolerance
  with notes when columns are out of order.

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


# ---------------------------------------------------------------------------
# bedGraph / wig anchors (methylation-aware path added 2026-05-19)
# ---------------------------------------------------------------------------

def _write_bedgraph(path: Path, rows: list[tuple]) -> Path:
    """Write a 4-col bedGraph (chrom/start/end/value), no header."""
    path.write_text("".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n" for r in rows))
    return path


def test_is_wig_data_line():
    assert va._is_wig_data_line("chr1\t100\t101\t0.5") is True
    assert va._is_wig_data_line("track name=meth") is False
    assert va._is_wig_data_line("browser dense") is False
    assert va._is_wig_data_line("fixedStep chrom=chr1 start=1 step=1") is False
    assert va._is_wig_data_line("variableStep chrom=chr1") is False
    assert va._is_wig_data_line("# comment") is False
    assert va._is_wig_data_line("") is False
    assert va._is_wig_data_line("   ") is False


def test_bedgraph_count_source_plain_text_in_region(tmp_path):
    # 3 of 4 rows overlap [100, 200); the 4th is on a different chrom.
    bg = _write_bedgraph(tmp_path / "sample.hg38.bedgraph", [
        ("chr1", 100, 101, 0.5),
        ("chr1", 150, 151, 0.8),
        ("chr1", 199, 200, 0.3),  # r_end > q_start? 200 > 100 yes; r_start < q_end? 199 < 200 yes
        ("chr2", 100, 101, 0.9),  # different chrom
    ])
    assert va.bedgraph_count_source(bg, "chr1", 100, 200) == 3


def test_bedgraph_count_source_excludes_out_of_region(tmp_path):
    # Rows must overlap [start, end). Boundary cases.
    bg = _write_bedgraph(tmp_path / "sample.hg38.bedgraph", [
        ("chr1", 50, 100, 0.1),    # r_end == q_start -> doesn't overlap (half-open)
        ("chr1", 100, 150, 0.2),   # r_start == q_start -> overlaps
        ("chr1", 195, 200, 0.3),   # r_start < q_end == 200 -> overlaps
        ("chr1", 200, 250, 0.4),   # r_start == q_end -> doesn't overlap (half-open)
        ("chr1", 1000, 1001, 0.5), # way out
    ])
    assert va.bedgraph_count_source(bg, "chr1", 100, 200) == 2


def test_bedgraph_count_source_skips_headers_and_comments(tmp_path):
    bg = tmp_path / "sample.hg38.bedgraph"
    bg.write_text(
        "#header comment\n"
        "track name=test\n"
        "browser dense\n"
        "chr1\t100\t101\t0.5\n"
        "chr1\t150\t151\t0.6\n"
    )
    assert va.bedgraph_count_source(bg, "chr1", 0, 1000) == 2


def test_bedgraph_count_source_handles_gzipped_input(tmp_path):
    # Plain-gzip (not bgzip+tabix). Linear-scan path.
    import gzip
    bg = tmp_path / "sample.hg38.bedgraph.gz"
    with gzip.open(bg, "wt") as fh:
        fh.write("chr3\t100\t101\t0.5\n")
        fh.write("chr3\t150\t151\t0.6\n")
        fh.write("chr3\t999\t1000\t0.7\n")
    assert va.bedgraph_count_source(bg, "chr3", 100, 200) == 2
    assert va.bedgraph_count_source(bg, "chr3", 0, 10000) == 3
    assert va.bedgraph_count_source(bg, "chr4", 0, 10000) == 0


def test_bedgraph_count_source_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="bedGraph track not found"):
        va.bedgraph_count_source(tmp_path / "does_not_exist.bg", "chr1", 0, 100)


def test_bedgraph_count_slice_decodes_gzipped_payload():
    # Mimics how igv_reports/datauri.py encodes a wig/bedGraph slice:
    # gzip(text) base64-encoded. verify_anchors only sees the gzipped
    # bytes after base64 decoding, so we test the bytes-in entry point.
    text = (
        "track name=meth\n"
        "chr1\t100\t101\t0.5\n"
        "chr1\t150\t151\t0.6\n"
        "chr1\t200\t201\t0.7\n"
    )
    assert va.bedgraph_count_slice(gzip.compress(text.encode())) == 3


def test_bedgraph_count_slice_falls_back_to_uncompressed():
    # Some create_report versions write small wig slices uncompressed —
    # the fallback path must accept raw text bytes.
    text = "chr1\t100\t101\t0.5\nchr1\t200\t201\t0.6\n"
    assert va.bedgraph_count_slice(text.encode()) == 2


def test_bedgraph_count_slice_zero_when_empty():
    # No data rows in the slice = silent empty-methylation-slice failure.
    # Caller (verify_one_html) compares to expected via decide_status.
    assert va.bedgraph_count_slice(gzip.compress(b"track name=meth\n")) == 0
    assert va.bedgraph_count_slice(b"") == 0


# ---------------------------------------------------------------------------
# Anchor schema: track_type column with backwards compat
# ---------------------------------------------------------------------------

def test_load_anchors_legacy_no_track_type_defaults_to_bam(tmp_path):
    # Pre-2026-05-19 anchor files lack the track_type column. Loader must
    # accept them and default each row to track_type='bam'.
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\ttumor\tchr2\t100\t200\t42\t\t\t\t\n"
    ))
    rows = va.load_anchors(p)
    assert rows[0].track_type == "bam"


def test_load_anchors_with_track_type_bedgraph(tmp_path):
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\ttrack_type\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\tmeth_track\tbedgraph\tchr2\t100\t200\t8\t\t\t\tDNMT3A_CpGs\n"
    ))
    rows = va.load_anchors(p)
    assert rows[0].track_type == "bedgraph"
    assert rows[0].expected == 8
    assert rows[0].notes == "DNMT3A_CpGs"


def test_load_anchors_rejects_unknown_track_type(tmp_path):
    p = _write_tsv(tmp_path, (
        "#sample\ttrack_name\ttrack_type\tchrom\tstart\tend\texpected\ttolerance\tmin\tmax\tnotes\n"
        "s1\tt1\tcraaam\tchr1\t0\t100\t5\t\t\t\t\n"
    ))
    with pytest.raises(SystemExit, match="unknown track_type 'craaam'"):
        va.load_anchors(p)


def test_write_load_round_trip_preserves_track_type(tmp_path):
    anchors_in = [
        va.AnchorRow(sample="s1", track_name="tumor", track_type="bam",
                     chrom="chr1", start=0, end=100, expected=42),
        va.AnchorRow(sample="s1", track_name="tumor.5mC", track_type="bedgraph",
                     chrom="chr1", start=0, end=100, expected=12),
    ]
    out = tmp_path / "anchors.tsv"
    va.write_anchors(anchors_in, out)
    rows = va.load_anchors(out)
    assert [r.track_type for r in rows] == ["bam", "bedgraph"]


# ---------------------------------------------------------------------------
# sample_bedgraph_paths: samplesheet → (track_name, bedgraph_path) iteration
# ---------------------------------------------------------------------------

def test_sample_bedgraph_paths_picks_bedgraph_from_extras():
    row = {"sample": "s1", "extra_tracks": "/data/x.5mC.bedgraph,/data/x.5hmC.bg"}
    pairs = va.sample_bedgraph_paths(row)
    assert pairs == [("x.5mC", Path("/data/x.5mC.bedgraph")),
                     ("x.5hmC", Path("/data/x.5hmC.bg"))]


def test_sample_bedgraph_paths_strips_gz_suffix_from_track_name():
    # Path.stem of foo.bedgraph.gz is "foo.bedgraph"; igv-reports renders
    # it as just "foo", so we strip one more level.
    row = {"sample": "s1", "extra_tracks": "/data/foo.bedgraph.gz"}
    pairs = va.sample_bedgraph_paths(row)
    assert pairs[0][0] == "foo"


def test_sample_bedgraph_paths_skips_non_bedgraph_extras():
    # bam/vcf in extra_tracks are NOT bedgraphs — sample_bam_paths handles them.
    row = {"sample": "s1", "extra_tracks": "/data/x.5mC.bedgraph,/data/y.bam,/data/z.vcf"}
    pairs = va.sample_bedgraph_paths(row)
    assert pairs == [("x.5mC", Path("/data/x.5mC.bedgraph"))]


def test_sample_bedgraph_paths_empty_when_no_extras():
    assert va.sample_bedgraph_paths({"sample": "s1"}) == []
    assert va.sample_bedgraph_paths({"sample": "s1", "extra_tracks": ""}) == []
