"""Smoke tests for verify_anchors.py — exercises the samtools subprocess
path and the end-to-end slice-decode-and-count flow against the committed
COLO829 BAM fixture.

Author: Samuel Ahuno
Purpose:
  The unit tests cover the parsers with synthetic inputs. These smoke
  tests confirm the verifier's samtools shellouts actually work end-to-end:

    1. `samtools_count` returns the right integer for a real BAM region.
    2. `samtools_index` produces a usable index on a fresh BAM.
    3. The full decode round-trip (read fixture BAM → base64 wrap into a
       fake data: URL → decode_track_slice writes it back out → re-index
       → re-count) preserves the original count exactly.

  Skipped (not failed) if samtools is unavailable — so this works in a
  CI sandbox without the SIF or PATH samtools.

Run:
  cd claude/skills/igv-reports
  pytest tests/smoke/ -v
"""

from __future__ import annotations

import base64
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import verify_anchors as va  # noqa: E402

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "tiny_colo829.hg38.bam"

# Anchor sanity counts — must match fixtures/README.md and build_fixtures.sh.
ANCHOR_CHR2 = ("chr2:25246500-25246501", 5)
ANCHOR_CHR7 = ("chr7:148884000-148884001", 9)


def _samtools_cmd() -> list[str] | None:
    """Resolve samtools the same way verify_anchors does, but return None
    instead of raising when nothing is available. Lets us SKIP gracefully."""
    try:
        return va.resolve_samtools(None)
    except SystemExit:
        # No SIF and no PATH samtools — environment can't run smoke tests.
        return None


@pytest.fixture(scope="module")
def samtools_cmd():
    cmd = _samtools_cmd()
    if cmd is None:
        pytest.skip("no samtools available (set SAMTOOLS_SIF or install samtools)")
    if not FIXTURE.exists():
        pytest.skip(
            f"fixture missing: {FIXTURE} — regenerate with "
            "bash tests/fixtures/build_fixtures.sh"
        )
    return cmd


# ---------------------------------------------------------------------------
# samtools_count + samtools_index against the committed fixture
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("region,expected", [ANCHOR_CHR2, ANCHOR_CHR7])
def test_samtools_count_matches_fixture_anchor(samtools_cmd, region, expected):
    observed = va.samtools_count(samtools_cmd, FIXTURE, region)
    assert observed == expected, (
        f"fixture anchor drift: {region} should be {expected}, got {observed}. "
        "Either the committed BAM was regenerated with different params "
        "(see tests/fixtures/build_fixtures.sh) or the count filter changed."
    )


def test_samtools_index_creates_usable_index(samtools_cmd, tmp_path):
    """Copy the fixture to tmp_path WITHOUT its .bai, then have verify_anchors
    re-index it. After indexing, samtools_count must succeed."""
    bam_copy = tmp_path / "no_bai.bam"
    shutil.copy(FIXTURE, bam_copy)
    # confirm no index exists yet
    assert not (tmp_path / "no_bai.bam.bai").exists()
    va.samtools_index(samtools_cmd, bam_copy)
    assert (tmp_path / "no_bai.bam.bai").exists()
    # count works now
    assert va.samtools_count(samtools_cmd, bam_copy, ANCHOR_CHR2[0]) == ANCHOR_CHR2[1]


def test_samtools_count_missing_bam_raises(samtools_cmd, tmp_path):
    """A missing BAM path should produce a clear RuntimeError, not crash silently.

    Note: samtools tolerates malformed region strings (treats them as unknown
    references and returns 0 with a stderr warning + exit 0). The only
    reliable error trigger is a missing/unreadable BAM file."""
    missing = tmp_path / "does_not_exist.bam"
    with pytest.raises(RuntimeError, match="samtools view -c failed"):
        va.samtools_count(samtools_cmd, missing, "chr1:1-100")


# ---------------------------------------------------------------------------
# Slice decode round-trip (the central correctness claim of verify_anchors)
# ---------------------------------------------------------------------------

def test_full_decode_roundtrip(samtools_cmd, tmp_path):
    """End-to-end: emulate what igv-reports does to embed a BAM slice in an
    HTML data: URL, then have verify_anchors decode it back out and confirm
    the read count is preserved.

    This is the critical correctness claim: if the verifier's slice decode
    silently corrupts the BAM bytes, every anchor verify would silently
    pass when it shouldn't. Catching that here means we trust the
    integration tests further down."""
    # 1. Read the committed BAM raw, wrap it in a data: URL exactly the
    #    way igv_reports/datauri.py does (mediatype application/gzip
    #    because BAM is BGZF gzip — see decision in datauri.get_data_uri).
    raw = FIXTURE.read_bytes()
    data_url = "data:application/gzip;base64," + base64.b64encode(raw).decode()

    # 2. Decode it via the production code path.
    decoded = tmp_path / "decoded.bam"
    va.decode_track_slice(data_url, decoded)

    # 3. Bytes must match exactly.
    assert decoded.read_bytes() == raw

    # 4. samtools should treat the decoded file as a real BAM — index it
    #    and count the same anchors.
    va.samtools_index(samtools_cmd, decoded)
    for region, expected in (ANCHOR_CHR2, ANCHOR_CHR7):
        assert va.samtools_count(samtools_cmd, decoded, region) == expected


def test_resolve_samtools_explicit_sif_missing(tmp_path):
    """Passing a non-existent SIF path explicitly must fail loudly, not fall
    back silently to PATH samtools."""
    fake_sif = tmp_path / "does_not_exist.sif"
    with pytest.raises(SystemExit, match="samtools SIF not found"):
        va.resolve_samtools(fake_sif)


def test_resolve_samtools_env_var(tmp_path, monkeypatch):
    """$SAMTOOLS_SIF env var honored when no --samtools-sif passed."""
    # Use a real-ish path that doesn't exist to verify it's the *path* the
    # env-resolution picks up, not some unrelated SIF.
    fake = tmp_path / "env_sif.sif"
    monkeypatch.setenv("SAMTOOLS_SIF", str(fake))
    with pytest.raises(SystemExit, match=str(fake)):
        va.resolve_samtools(None)
