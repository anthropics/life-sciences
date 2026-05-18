"""Unit tests for generate_tracks_json.py — annotation-default resolver.

Author: Samuel Ahuno
Purpose:
  Exercises the `default:` shortcut path added in the methylation-pathway
  polish round. Without these tests a future Claude session could easily
  break the resolver by adding a 6th key without updating the lookup.

Covers:
  * Happy path: each known default key resolves against a synthetic cfg.
  * hg38 gencode-sibling preference (.gff3.gz over .gtf.gz when present).
  * indexURL: included when .tbi exists; omitted otherwise.
  * Unknown default key -> SystemExit with valid-keys hint.
  * Missing genome in cfg -> SystemExit.
  * Missing path on disk -> SystemExit.
  * build_annotation_tracks() routes `default:` entries through the resolver
    and preserves backwards compat for explicit `url:` entries.
  * `default:` entry without top-level `genome:` -> SystemExit.

Run:
  cd igv-reports-skill && pytest tests/unit/test_generate_tracks_json.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import generate_tracks_json as g  # noqa: E402


def _fake_cfg(genome: str, paths: dict[str, str]) -> dict:
    """Build a minimal databases YAML mirror keyed by genome.

    `paths` maps YAML-keys (CpGIslands, gtf, repMaskerBed, EPDnewCoding,
    EPDnewNonCoding) to filesystem paths."""
    return {"reference_genomes": {"local": {genome: paths}}}


def _touch(path: Path) -> Path:
    """Create an empty file at `path`, parents auto-created."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path


# ----- happy path: each known default key resolves -----


def test_resolve_cgi(tmp_path):
    cgi = _touch(tmp_path / "hg38_CpGIslands.bed")
    cfg = _fake_cfg("hg38", {"CpGIslands": str(cgi)})
    t = g.resolve_annotation_default("cgi", "hg38", cfg)
    assert t["url"] == str(cgi)
    assert t["display_name"] == "CpG islands"
    assert t["format"] == "bed"
    assert t["color"] == "rgb(0,158,115)"
    assert "indexURL" not in t  # no .tbi alongside


def test_resolve_repmasker(tmp_path):
    rmsk = _touch(tmp_path / "rmsk.bed.gz")
    _touch(tmp_path / "rmsk.bed.gz.tbi")
    cfg = _fake_cfg("hg38", {"repMaskerBed": str(rmsk)})
    t = g.resolve_annotation_default("repmasker", "hg38", cfg)
    assert t["url"] == str(rmsk)
    assert t["indexURL"] == str(rmsk) + ".tbi"
    assert t["displayMode"] == "COLLAPSED"


def test_resolve_gencode_hg38_prefers_gff3_sibling(tmp_path):
    # When the YAML's `gtf` points at a .gtf.gz, but a sibling
    # gencode.v47.annotation.gff3.gz + .tbi exists in the same dir,
    # the resolver should switch to the bgzip+tabix .gff3.gz file.
    gtf = _touch(tmp_path / "gencode.v47.annotation.gtf.gz")
    sibling = _touch(tmp_path / "gencode.v47.annotation.gff3.gz")
    _touch(tmp_path / "gencode.v47.annotation.gff3.gz.tbi")
    cfg = _fake_cfg("hg38", {"gtf": str(gtf)})
    t = g.resolve_annotation_default("gencode", "hg38", cfg)
    assert t["url"] == str(sibling), "expected hg38 gencode to prefer .gff3.gz sibling"
    assert t["indexURL"] == str(sibling) + ".tbi"


def test_resolve_gencode_mm10_uses_gtf(tmp_path):
    # The sibling-preference logic only fires for hg38. For mm10 the
    # resolver should use the YAML-named gtf path verbatim.
    gtf = _touch(tmp_path / "gencode.vM25.annotation.gtf.gz")
    cfg = _fake_cfg("mm10", {"gtf": str(gtf)})
    t = g.resolve_annotation_default("gencode", "mm10", cfg)
    assert t["url"] == str(gtf)


def test_resolve_epdnew_coding_and_noncoding(tmp_path):
    coding = _touch(tmp_path / "Hs_EPDnew.hg38.bed.gz")
    noncoding = _touch(tmp_path / "HsNC_EPDnew.hg38.bed.gz")
    cfg = _fake_cfg("hg38", {
        "EPDnewCoding": str(coding),
        "EPDnewNonCoding": str(noncoding),
    })
    tc = g.resolve_annotation_default("epdnew_coding", "hg38", cfg)
    tn = g.resolve_annotation_default("epdnew_noncoding", "hg38", cfg)
    assert tc["url"] == str(coding)
    assert tn["url"] == str(noncoding)
    # Distinct Okabe-Ito colors so coding vs non-coding read separately.
    assert tc["color"] != tn["color"]


# ----- error paths -----


def test_unknown_default_key_lists_valid_keys():
    with pytest.raises(SystemExit) as exc:
        g.resolve_annotation_default("DOES_NOT_EXIST", "hg38", _fake_cfg("hg38", {}))
    # Error should enumerate the valid keys so the user can fix the typo
    # without having to read the source.
    msg = str(exc.value)
    for key in ("cgi", "gencode", "repmasker", "epdnew_coding", "epdnew_noncoding"):
        assert key in msg


def test_missing_genome_in_cfg(tmp_path):
    cgi = _touch(tmp_path / "hg38_CpGIslands.bed")
    cfg = _fake_cfg("hg38", {"CpGIslands": str(cgi)})
    with pytest.raises(SystemExit, match="no entry for genome 'GRCh37'"):
        g.resolve_annotation_default("cgi", "GRCh37", cfg)


def test_missing_yaml_key_for_genome(tmp_path):
    # mm39 famously has no repMaskerBed configured — `default: repmasker`
    # must fail with a clear error rather than silently emitting no track.
    cfg = _fake_cfg("mm39", {"CpGIslands": "/tmp/fake_cgi"})
    with pytest.raises(SystemExit, match="repMaskerBed"):
        g.resolve_annotation_default("repmasker", "mm39", cfg)


def test_missing_path_on_disk(tmp_path):
    cfg = _fake_cfg("hg38", {"CpGIslands": str(tmp_path / "nonexistent.bed")})
    with pytest.raises(SystemExit, match="resolved path missing on disk"):
        g.resolve_annotation_default("cgi", "hg38", cfg)


# ----- build_annotation_tracks() integration -----


def test_build_annotation_tracks_shortcut(tmp_path):
    cgi = _touch(tmp_path / "hg38_CpGIslands.bed")
    cfg = _fake_cfg("hg38", {"CpGIslands": str(cgi)})
    spec = {"genome": "hg38", "annotation": [{"default": "cgi"}]}
    out = g.build_annotation_tracks(spec, tmp_path, cfg)
    assert len(out) == 1
    assert out[0]["name"] == "CpG islands"
    assert out[0]["url"] == str(cgi)
    assert out[0]["type"] == "annotation"
    assert out[0]["color"] == "rgb(0,158,115)"


def test_build_annotation_tracks_shortcut_with_overrides(tmp_path):
    # The user can override the canned display name + color while still
    # using `default:` for path resolution.
    cgi = _touch(tmp_path / "hg38_CpGIslands.bed")
    cfg = _fake_cfg("hg38", {"CpGIslands": str(cgi)})
    spec = {
        "genome": "hg38",
        "annotation": [{
            "default": "cgi",
            "name": "My CpG view",
            "color": "rgb(0,0,0)",
            "displayMode": "COLLAPSED",
        }],
    }
    out = g.build_annotation_tracks(spec, tmp_path, cfg)
    assert out[0]["name"] == "My CpG view"
    assert out[0]["color"] == "rgb(0,0,0)"
    assert out[0]["displayMode"] == "COLLAPSED"
    # url still resolved by the shortcut
    assert out[0]["url"] == str(cgi)


def test_build_annotation_tracks_explicit_path_unchanged(tmp_path):
    # Backwards-compat: an explicit `url:` entry must not need a cfg and
    # must produce the same shape as before this round of changes.
    explicit = _touch(tmp_path / "my_custom.bed")
    spec = {
        "annotation": [{
            "name": "My custom track",
            "url": str(explicit),
            "format": "bed",
            "color": "rgb(1,2,3)",
        }],
    }
    out = g.build_annotation_tracks(spec, tmp_path, {})
    assert out[0] == {
        "name": "My custom track",
        "url": str(explicit),
        "format": "bed",
        "type": "annotation",
        "displayMode": "EXPANDED",
        "color": "rgb(1,2,3)",
    }


def test_build_annotation_tracks_mixed(tmp_path):
    # Explicit + shortcut entries can coexist; order is preserved.
    cgi = _touch(tmp_path / "hg38_CpGIslands.bed")
    explicit = _touch(tmp_path / "custom.bed")
    cfg = _fake_cfg("hg38", {"CpGIslands": str(cgi)})
    spec = {
        "genome": "hg38",
        "annotation": [
            {"name": "Custom first", "url": str(explicit)},
            {"default": "cgi"},
        ],
    }
    out = g.build_annotation_tracks(spec, tmp_path, cfg)
    assert [t["name"] for t in out] == ["Custom first", "CpG islands"]


def test_shortcut_without_top_level_genome(tmp_path):
    spec = {"annotation": [{"default": "cgi"}]}  # missing `genome:`
    with pytest.raises(SystemExit, match="top-level `genome:`"):
        g.build_annotation_tracks(spec, tmp_path, {})
