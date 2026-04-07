#!/usr/bin/env python3
"""
analyze_fasta.py - General-purpose FASTA file analyzer.
Automatically detects nucleotide vs protein sequences and generates a comprehensive report.

Usage:
    python3 analyze_fasta.py input.fasta
    python3 analyze_fasta.py input.fasta --json
    python3 analyze_fasta.py input.fasta --html
"""

import sys
import json
import html as html_mod
import argparse
from pathlib import Path
from collections import Counter
from datetime import datetime

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ──────────────────────────────────────────────
# Type detection
# ──────────────────────────────────────────────

def detect_sequence_type(seq_str):
    """Detect whether a sequence is nucleotide or protein."""
    nuc_chars = set("ATCGUNatcgun")
    sample = seq_str[:500].replace("-", "").replace(".", "")
    nuc_ratio = sum(1 for c in sample if c in nuc_chars) / max(len(sample), 1)
    return "nucleotide" if nuc_ratio > 0.85 else "protein"


# ──────────────────────────────────────────────
# Nucleotide analysis
# ──────────────────────────────────────────────

def analyze_nucleotide(record):
    seq = record.seq
    seq_str = str(seq).upper()
    length = len(seq)
    comp = Counter(seq_str)

    result = {
        "id": record.id,
        "description": record.description,
        "length_bp": length,
        "gc_content": round(gc_fraction(seq) * 100, 2),
        "composition": {
            "A": comp.get("A", 0),
            "T": comp.get("T", 0),
            "G": comp.get("G", 0),
            "C": comp.get("C", 0),
            "N": comp.get("N", 0),
        },
        "composition_pct": {
            "A": round(comp.get("A", 0) / max(length, 1) * 100, 2),
            "T": round(comp.get("T", 0) / max(length, 1) * 100, 2),
            "G": round(comp.get("G", 0) / max(length, 1) * 100, 2),
            "C": round(comp.get("C", 0) / max(length, 1) * 100, 2),
        },
        "at_content": round((comp.get("A", 0) + comp.get("T", 0)) / max(length, 1) * 100, 2),
    }

    # Dinucleotides
    dinuc = Counter(seq_str[i:i+2] for i in range(len(seq_str) - 1))
    top_dinuc = dinuc.most_common(5)
    result["top_dinucleotides"] = {d: c for d, c in top_dinuc}

    # Find simple ORFs (ATG ... stop)
    stops = {"TAA", "TAG", "TGA"}
    orfs = []
    for frame in range(3):
        i = frame
        while i < length - 2:
            codon = seq_str[i:i+3]
            if codon == "ATG":
                start = i
                j = i + 3
                while j < length - 2:
                    c = seq_str[j:j+3]
                    if c in stops:
                        orf_len = j - start + 3
                        if orf_len >= 300:
                            orfs.append({
                                "frame": f"+{frame+1}",
                                "start": start + 1,
                                "end": j + 3,
                                "length_bp": orf_len,
                                "length_aa": orf_len // 3,
                            })
                        break
                    j += 3
                i = j + 3
            else:
                i += 3

    result["orfs_found"] = len(orfs)
    result["orfs"] = orfs[:10]

    try:
        mw = molecular_weight(seq, "DNA")
        result["molecular_weight_da"] = round(mw, 1)
    except Exception:
        pass

    return result


# ──────────────────────────────────────────────
# Protein analysis
# ──────────────────────────────────────────────

def analyze_protein(record):
    seq_str = str(record.seq).upper().replace("X", "").replace("*", "")
    length = len(seq_str)
    comp = Counter(seq_str)

    result = {
        "id": record.id,
        "description": record.description,
        "length_aa": length,
        "composition": {aa: comp.get(aa, 0) for aa in sorted(comp.keys())},
        "composition_pct": {
            aa: round(count / max(length, 1) * 100, 2)
            for aa, count in sorted(comp.items())
        },
    }

    try:
        pa = ProteinAnalysis(seq_str)
        result["molecular_weight_da"] = round(pa.molecular_weight(), 1)
        result["isoelectric_point"] = round(pa.isoelectric_point(), 2)
        result["aromaticity"] = round(pa.aromaticity(), 4)
        result["instability_index"] = round(pa.instability_index(), 2)
        result["stability"] = "stable" if pa.instability_index() < 40 else "unstable"
        result["gravy"] = round(pa.gravy(), 4)
        result["hydrophobicity"] = "hydrophilic" if pa.gravy() < 0 else "hydrophobic"

        helix, turn, sheet = pa.secondary_structure_fraction()
        result["secondary_structure_pct"] = {
            "helix": round(helix * 100, 1),
            "turn": round(turn * 100, 1),
            "sheet": round(sheet * 100, 1),
        }

        aa_pct = pa.amino_acids_percent
        charged = aa_pct.get("R", 0) + aa_pct.get("K", 0) + aa_pct.get("D", 0) + aa_pct.get("E", 0)
        result["charged_residues_pct"] = round(charged, 1)

        aromatic = aa_pct.get("F", 0) + aa_pct.get("W", 0) + aa_pct.get("Y", 0)
        result["aromatic_residues_pct"] = round(aromatic, 1)

    except Exception as e:
        result["analysis_error"] = str(e)

    return result


# ──────────────────────────────────────────────
# Main analysis
# ──────────────────────────────────────────────

def analyze_fasta(filepath):
    path = Path(filepath)
    if not path.exists():
        return {"error": f"File not found: {filepath}"}

    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        return {"error": "No sequences found in FASTA file"}

    seq_type = detect_sequence_type(str(records[0].seq))

    report = {
        "file": path.name,
        "file_path": str(path.resolve()),
        "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "total_sequences": len(records),
        "sequence_type": seq_type,
        "sequences": [],
        "summary": {},
    }

    lengths = []
    for record in records:
        if seq_type == "nucleotide":
            analysis = analyze_nucleotide(record)
        else:
            analysis = analyze_protein(record)
        report["sequences"].append(analysis)
        lengths.append(len(record.seq))

    report["summary"] = {
        "total_sequences": len(records),
        "total_residues": sum(lengths),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "avg_length": round(sum(lengths) / len(lengths), 1),
        "n50": calculate_n50(lengths),
    }

    if seq_type == "nucleotide":
        gc_values = [s.get("gc_content", 0) for s in report["sequences"]]
        report["summary"]["avg_gc_content"] = round(sum(gc_values) / len(gc_values), 2)
        report["summary"]["total_orfs"] = sum(s.get("orfs_found", 0) for s in report["sequences"])

    return report


def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    running = 0
    for length in sorted_lengths:
        running += length
        if running >= total / 2:
            return length
    return 0


# ──────────────────────────────────────────────
# Text format
# ──────────────────────────────────────────────

def format_text_report(report):
    lines = []
    lines.append("=" * 60)
    lines.append(f"  FASTA ANALYSIS: {report['file']}")
    lines.append("=" * 60)
    lines.append("")

    s = report["summary"]
    lines.append(f"Sequence type:       {report['sequence_type']}")
    lines.append(f"Total sequences:     {s['total_sequences']}")
    lines.append(f"Total residues:      {s['total_residues']:,}")
    lines.append(f"Min length:          {s['min_length']:,}")
    lines.append(f"Max length:          {s['max_length']:,}")
    lines.append(f"Avg length:          {s['avg_length']:,}")
    lines.append(f"N50:                 {s['n50']:,}")

    if report["sequence_type"] == "nucleotide":
        lines.append(f"Avg GC content:      {s.get('avg_gc_content', 'N/A')}%")
        lines.append(f"ORFs found:          {s.get('total_orfs', 0)}")

    lines.append("")
    lines.append("-" * 60)

    for i, seq in enumerate(report["sequences"][:20]):
        lines.append("")
        lines.append(f"  Sequence {i+1}: {seq.get('id', 'unknown')}")
        lines.append(f"  {seq.get('description', '')}")
        lines.append("")

        if report["sequence_type"] == "nucleotide":
            lines.append(f"    Length:           {seq['length_bp']:,} bp")
            lines.append(f"    GC content:       {seq['gc_content']}%")
            lines.append(f"    AT content:       {seq['at_content']}%")
            comp = seq["composition_pct"]
            lines.append(f"    Composition:      A={comp['A']}%  T={comp['T']}%  G={comp['G']}%  C={comp['C']}%")
            if seq.get("molecular_weight_da"):
                lines.append(f"    Molecular weight: {seq['molecular_weight_da']:,.1f} Da")
            if seq["orfs_found"] > 0:
                lines.append(f"    ORFs (>=100aa):   {seq['orfs_found']}")
                for orf in seq["orfs"][:5]:
                    lines.append(f"      {orf['frame']} pos {orf['start']}-{orf['end']} ({orf['length_aa']} aa)")
        else:
            lines.append(f"    Length:            {seq['length_aa']} aa")
            if seq.get("molecular_weight_da"):
                lines.append(f"    Molecular weight:  {seq['molecular_weight_da']:,.1f} Da")
            if seq.get("isoelectric_point"):
                lines.append(f"    Isoelectric point: {seq['isoelectric_point']}")
            if seq.get("instability_index"):
                lines.append(f"    Instability index: {seq['instability_index']} ({seq.get('stability', '')})")
            if seq.get("gravy"):
                lines.append(f"    GRAVY:             {seq['gravy']} ({seq.get('hydrophobicity', '')})")
            if seq.get("aromaticity"):
                lines.append(f"    Aromaticity:       {seq['aromaticity']}")
            if seq.get("secondary_structure_pct"):
                ss = seq["secondary_structure_pct"]
                lines.append(f"    Secondary struct:  helix={ss['helix']}%  sheet={ss['sheet']}%  turn={ss['turn']}%")
            if seq.get("charged_residues_pct"):
                lines.append(f"    Charged residues:  {seq['charged_residues_pct']}%")

        lines.append("    " + "-" * 40)

    if len(report["sequences"]) > 20:
        lines.append(f"\n  ... and {len(report['sequences']) - 20} more sequences.")

    lines.append("")
    lines.append("=" * 60)
    return "\n".join(lines)


# ──────────────────────────────────────────────
# HTML format
# ──────────────────────────────────────────────

def generate_html_report(report, json_str):
    """Generate a complete HTML report with analysis results, JSON viewer, and methodology section."""

    s = report["summary"]
    seq_type = report["sequence_type"]
    is_nuc = seq_type == "nucleotide"

    summary_rows = f"""
        <tr><td>Sequence type</td><td><span class="badge {'nuc' if is_nuc else 'prot'}">{seq_type}</span></td></tr>
        <tr><td>Total sequences</td><td>{s['total_sequences']}</td></tr>
        <tr><td>Total residues</td><td>{s['total_residues']:,} {'bp' if is_nuc else 'aa'}</td></tr>
        <tr><td>Min length</td><td>{s['min_length']:,}</td></tr>
        <tr><td>Max length</td><td>{s['max_length']:,}</td></tr>
        <tr><td>Avg length</td><td>{s['avg_length']:,}</td></tr>
        <tr><td>N50</td><td>{s['n50']:,}</td></tr>"""

    if is_nuc:
        summary_rows += f"""
        <tr><td>Avg GC content</td><td>{s.get('avg_gc_content', 'N/A')}%</td></tr>
        <tr><td>ORFs found (&ge;100 aa)</td><td>{s.get('total_orfs', 0)}</td></tr>"""

    seq_cards = ""
    for i, seq in enumerate(report["sequences"][:20]):
        sid = html_mod.escape(seq.get("id", "unknown"))
        sdesc = html_mod.escape(seq.get("description", ""))

        if is_nuc:
            comp = seq["composition_pct"]
            gc = seq["gc_content"]
            at = seq["at_content"]
            bar_a, bar_t, bar_g, bar_c = comp["A"], comp["T"], comp["G"], comp["C"]

            orfs_html = ""
            if seq["orfs_found"] > 0:
                orf_rows = ""
                for orf in seq["orfs"][:8]:
                    orf_rows += f'<tr><td>{orf["frame"]}</td><td>{orf["start"]:,} - {orf["end"]:,}</td><td>{orf["length_bp"]:,} bp</td><td>{orf["length_aa"]:,} aa</td></tr>'
                orfs_html = f"""
                <h4>ORFs detected ({seq['orfs_found']})</h4>
                <table class="inner-table">
                    <thead><tr><th>Frame</th><th>Position</th><th>Length</th><th>Protein</th></tr></thead>
                    <tbody>{orf_rows}</tbody>
                </table>"""

            dinuc_bars = ""
            for dn, count in (seq.get("top_dinucleotides") or {}).items():
                max_dn = max(seq["top_dinucleotides"].values()) if seq["top_dinucleotides"] else 1
                w = count / max_dn * 100
                dinuc_bars += f'<div class="mini-bar-row"><span class="mini-label">{dn}</span><div class="mini-track"><div class="mini-fill nuc-fill" style="width:{w}%"></div></div><span class="mini-val">{count:,}</span></div>'

            seq_cards += f"""
            <div class="seq-card">
                <div class="seq-header">
                    <span class="seq-num">#{i+1}</span>
                    <div><h3>{sid}</h3><p class="seq-desc">{sdesc}</p></div>
                </div>
                <div class="seq-grid">
                    <div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Length</span><span class="metric-value">{seq['length_bp']:,} bp</span></div>
                            <div class="metric"><span class="metric-label">GC content</span><span class="metric-value">{gc}%</span></div>
                            <div class="metric"><span class="metric-label">AT content</span><span class="metric-value">{at}%</span></div>
                        </div>
                        <h4>Base composition</h4>
                        <div class="comp-bar">
                            <div class="comp-seg seg-a" style="width:{bar_a}%" title="A: {bar_a}%">A</div>
                            <div class="comp-seg seg-t" style="width:{bar_t}%" title="T: {bar_t}%">T</div>
                            <div class="comp-seg seg-g" style="width:{bar_g}%" title="G: {bar_g}%">G</div>
                            <div class="comp-seg seg-c" style="width:{bar_c}%" title="C: {bar_c}%">C</div>
                        </div>
                        <div class="comp-legend">
                            <span><i class="dot dot-a"></i>A {bar_a}%</span>
                            <span><i class="dot dot-t"></i>T {bar_t}%</span>
                            <span><i class="dot dot-g"></i>G {bar_g}%</span>
                            <span><i class="dot dot-c"></i>C {bar_c}%</span>
                        </div>
                    </div>
                    <div><h4>Top dinucleotides</h4>{dinuc_bars}</div>
                </div>
                {orfs_html}
            </div>"""
        else:
            mw = seq.get("molecular_weight_da", "N/A")
            mw_str = f"{mw:,.1f} Da" if isinstance(mw, (int, float)) else mw
            pi = seq.get("isoelectric_point", "N/A")
            stab = seq.get("instability_index", "N/A")
            stab_label = seq.get("stability", "")
            gravy = seq.get("gravy", "N/A")
            hydro = seq.get("hydrophobicity", "")
            arom = seq.get("aromaticity", "N/A")
            ss = seq.get("secondary_structure_pct") or {}
            charged = seq.get("charged_residues_pct", "N/A")
            aromatic_pct = seq.get("aromatic_residues_pct", "N/A")

            ss_html = ""
            if ss:
                ss_html = f"""
                <h4>Predicted secondary structure</h4>
                <div class="ss-bars">
                    <div class="ss-row"><span class="ss-label">Helix</span><div class="mini-track"><div class="mini-fill helix-fill" style="width:{ss.get('helix',0)}%"></div></div><span class="mini-val">{ss.get('helix',0)}%</span></div>
                    <div class="ss-row"><span class="ss-label">Sheet</span><div class="mini-track"><div class="mini-fill sheet-fill" style="width:{ss.get('sheet',0)}%"></div></div><span class="mini-val">{ss.get('sheet',0)}%</span></div>
                    <div class="ss-row"><span class="ss-label">Turn</span><div class="mini-track"><div class="mini-fill turn-fill" style="width:{ss.get('turn',0)}%"></div></div><span class="mini-val">{ss.get('turn',0)}%</span></div>
                </div>"""

            seq_cards += f"""
            <div class="seq-card">
                <div class="seq-header">
                    <span class="seq-num">#{i+1}</span>
                    <div><h3>{sid}</h3><p class="seq-desc">{sdesc}</p></div>
                </div>
                <div class="seq-grid">
                    <div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Length</span><span class="metric-value">{seq['length_aa']:,} aa</span></div>
                            <div class="metric"><span class="metric-label">Molecular weight</span><span class="metric-value">{mw_str}</span></div>
                            <div class="metric"><span class="metric-label">Isoelectric point</span><span class="metric-value">{pi}</span></div>
                        </div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Stability</span><span class="metric-value">{stab} <small>({stab_label})</small></span></div>
                            <div class="metric"><span class="metric-label">GRAVY</span><span class="metric-value">{gravy} <small>({hydro})</small></span></div>
                            <div class="metric"><span class="metric-label">Aromaticity</span><span class="metric-value">{arom}</span></div>
                        </div>
                        <div class="metric-group">
                            <div class="metric"><span class="metric-label">Charged residues</span><span class="metric-value">{charged}%</span></div>
                            <div class="metric"><span class="metric-label">Aromatic residues</span><span class="metric-value">{aromatic_pct}%</span></div>
                        </div>
                    </div>
                    <div>{ss_html}</div>
                </div>
            </div>"""

    overflow_note = ""
    if len(report["sequences"]) > 20:
        overflow_note = f'<p class="overflow-note">Showing 20 of {len(report["sequences"])} sequences.</p>'

    json_escaped = html_mod.escape(json_str)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FASTA Analysis - {html_mod.escape(report['file'])}</title>
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
        :root {{
            --bg:#fff; --surface:#f8f9fa; --surface-2:#f1f3f5; --border:#dee2e6; --border-lt:#e9ecef;
            --text:#212529; --text-2:#495057; --text-dim:#6c757d; --text-cap:#adb5bd;
            --accent:#4263eb; --accent-light:#dbe4ff; --accent-dark:#364fc7;
            --green:#2b8a3e; --red:#c92a2a; --amber:#e67700;
            --green-bg:#ebfbee; --red-bg:#fff5f5; --amber-bg:#fff9db;
        }}
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:'Inter',system-ui,sans-serif; background:var(--bg); color:var(--text); line-height:1.7; }}
        .top-bar {{ background:var(--text); padding:0.5rem 2rem; display:flex; align-items:center; justify-content:space-between; }}
        .top-bar-logo {{ font-weight:700; font-size:0.95rem; color:#fff; letter-spacing:0.02em; }}
        .top-bar-info {{ font-size:0.72rem; color:rgba(255,255,255,0.4); }}
        .accent-strip {{ height:3px; background:var(--accent); }}
        .container {{ max-width:920px; margin:0 auto; padding:2rem 2rem 4rem; }}

        .hero {{ text-align:center; padding:2rem 0 1.8rem; border-bottom:1px solid var(--border); margin-bottom:2.5rem; }}
        .hero h1 {{ font-size:1.4rem; font-weight:700; color:var(--text); margin-bottom:0.5rem; }}
        .hero-meta {{ font-size:0.78rem; color:var(--text-cap); margin-top:0.4rem; }}
        .hero-file {{ font-family:'JetBrains Mono',monospace; font-size:0.85rem; color:var(--accent); background:var(--accent-light); padding:0.15rem 0.6rem; border-radius:3px; }}

        .sn {{ display:inline-flex; align-items:center; gap:0.5rem; font-size:0.68rem; font-weight:600; color:var(--accent); text-transform:uppercase; letter-spacing:0.1em; margin-bottom:0.5rem; }}
        .sn::before {{ content:''; width:20px; height:2px; background:var(--accent); border-radius:1px; }}
        h2 {{ font-size:1.2rem; font-weight:700; color:var(--text); margin-bottom:0.6rem; }}
        h3 {{ font-size:0.95rem; font-weight:600; color:var(--text); margin-bottom:0.3rem; }}
        h4 {{ font-size:0.78rem; font-weight:600; color:var(--text-dim); text-transform:uppercase; letter-spacing:0.05em; margin:1rem 0 0.4rem; }}
        p {{ color:var(--text-2); font-size:0.9rem; margin-bottom:0.6rem; }}
        section {{ margin-bottom:2.5rem; }}
        small {{ color:var(--text-dim); }}

        .badge {{ display:inline-block; padding:0.1rem 0.45rem; border-radius:3px; font-family:'JetBrains Mono',monospace; font-size:0.7rem; font-weight:600; }}
        .badge.nuc {{ background:var(--accent-light); color:var(--accent); }}
        .badge.prot {{ background:var(--amber-bg); color:var(--amber); }}

        .summary-table {{ width:100%; border-collapse:collapse; margin-top:0.8rem; }}
        .summary-table td {{ padding:0.5rem 0.7rem; font-size:0.85rem; border-bottom:1px solid var(--border-lt); }}
        .summary-table td:first-child {{ font-weight:600; color:var(--text); width:45%; }}
        .summary-table td:last-child {{ color:var(--text-2); font-family:'JetBrains Mono',monospace; font-size:0.8rem; }}
        .summary-table tr:hover td {{ background:var(--surface); }}

        .seq-card {{ background:var(--surface); border:1px solid var(--border-lt); border-radius:6px; padding:1.4rem 1.6rem; margin-bottom:1rem; }}
        .seq-header {{ display:flex; align-items:baseline; gap:0.7rem; margin-bottom:0.8rem; padding-bottom:0.7rem; border-bottom:1px solid var(--border-lt); }}
        .seq-num {{ font-size:1.3rem; font-weight:700; color:var(--border); }}
        .seq-desc {{ font-size:0.75rem; color:var(--text-cap); margin:0; word-break:break-all; }}
        .seq-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:1.3rem; }}
        @media(max-width:640px) {{ .seq-grid {{ grid-template-columns:1fr; }} }}

        .metric-group {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(130px,1fr)); gap:0.4rem; margin-bottom:0.6rem; }}
        .metric {{ background:var(--bg); border:1px solid var(--border-lt); border-radius:4px; padding:0.45rem 0.6rem; }}
        .metric-label {{ display:block; font-size:0.65rem; font-weight:600; text-transform:uppercase; letter-spacing:0.04em; color:var(--text-cap); }}
        .metric-value {{ font-family:'JetBrains Mono',monospace; font-size:0.85rem; font-weight:600; color:var(--text); }}

        .comp-bar {{ display:flex; height:24px; border-radius:4px; overflow:hidden; margin:0.4rem 0; }}
        .comp-seg {{ display:flex; align-items:center; justify-content:center; font-size:0.68rem; font-weight:600; color:#fff; min-width:18px; }}
        .seg-a {{ background:#4263eb; }} .seg-t {{ background:#e03131; }} .seg-g {{ background:#2f9e44; }} .seg-c {{ background:#e8590c; }}
        .comp-legend {{ display:flex; gap:0.8rem; font-size:0.72rem; color:var(--text-dim); flex-wrap:wrap; }}
        .dot {{ display:inline-block; width:7px; height:7px; border-radius:2px; margin-right:3px; vertical-align:middle; }}
        .dot-a {{ background:#4263eb; }} .dot-t {{ background:#e03131; }} .dot-g {{ background:#2f9e44; }} .dot-c {{ background:#e8590c; }}

        .mini-bar-row {{ display:flex; align-items:center; gap:0.4rem; margin-bottom:0.25rem; }}
        .mini-label {{ font-family:'JetBrains Mono',monospace; font-size:0.7rem; width:26px; text-align:right; color:var(--text-dim); }}
        .mini-track {{ flex:1; height:14px; background:var(--surface-2); border-radius:3px; overflow:hidden; }}
        .mini-fill {{ height:100%; border-radius:3px; }}
        .nuc-fill {{ background:var(--accent); opacity:0.5; }}
        .mini-val {{ font-family:'JetBrains Mono',monospace; font-size:0.65rem; color:var(--text-cap); width:45px; }}

        .inner-table {{ width:100%; border-collapse:collapse; margin-top:0.4rem; font-size:0.8rem; }}
        .inner-table th {{ text-align:left; padding:0.35rem 0.5rem; font-size:0.65rem; font-weight:600; letter-spacing:0.05em; text-transform:uppercase; color:var(--text-dim); background:var(--surface-2); border-bottom:2px solid var(--accent); }}
        .inner-table td {{ padding:0.3rem 0.5rem; border-bottom:1px solid var(--border-lt); color:var(--text-2); font-family:'JetBrains Mono',monospace; font-size:0.75rem; }}

        .ss-bars {{ margin:0.25rem 0; }}
        .ss-row {{ display:flex; align-items:center; gap:0.4rem; margin-bottom:0.25rem; }}
        .ss-label {{ font-size:0.72rem; width:50px; text-align:right; color:var(--text-dim); }}
        .helix-fill {{ background:var(--red); opacity:0.5; }}
        .sheet-fill {{ background:var(--accent); opacity:0.5; }}
        .turn-fill {{ background:var(--green); opacity:0.5; }}

        .overflow-note {{ font-size:0.8rem; color:var(--text-cap); font-style:italic; text-align:center; margin-top:0.8rem; }}

        .json-toggle {{ background:var(--text); color:#fff; border:none; padding:0.4rem 0.9rem; border-radius:4px; font-family:'Inter',sans-serif; font-size:0.8rem; font-weight:600; cursor:pointer; margin-bottom:0.7rem; }}
        .json-toggle:hover {{ background:#343a40; }}
        .json-block {{ background:#212529; border-radius:4px; padding:1rem 1.3rem; overflow-x:auto; display:none; }}
        .json-block pre {{ font-family:'JetBrains Mono',monospace; font-size:0.72rem; line-height:1.6; color:#e9ecef; white-space:pre-wrap; word-break:break-word; }}

        .how-section {{ background:var(--surface); border:1px solid var(--border-lt); border-radius:6px; padding:1.8rem 2rem; position:relative; }}
        .how-section::before {{ content:''; position:absolute; top:0; left:0; width:100%; height:3px; background:var(--accent); border-radius:6px 6px 0 0; }}
        .how-section h2 {{ margin-bottom:0.8rem; }}
        .how-grid {{ display:grid; grid-template-columns:1fr 1fr; gap:1.2rem; margin-top:0.8rem; }}
        @media(max-width:640px) {{ .how-grid {{ grid-template-columns:1fr; }} }}
        .how-card {{ background:var(--bg); border:1px solid var(--border-lt); border-radius:4px; padding:1rem; }}
        .how-card h3 {{ font-size:0.85rem; margin-bottom:0.35rem; }}
        .how-card p {{ font-size:0.82rem; margin-bottom:0.25rem; }}
        .how-card code {{ font-family:'JetBrains Mono',monospace; font-size:0.76em; background:var(--surface-2); padding:0.08em 0.3em; border-radius:2px; }}
        .how-card ul {{ list-style:none; padding:0; margin:0.25rem 0 0; }}
        .how-card ul li {{ font-size:0.8rem; padding:0.15rem 0 0.15rem 0.8rem; position:relative; color:var(--text-2); }}
        .how-card ul li::before {{ content:''; position:absolute; left:0; top:0.55rem; width:4px; height:4px; border-radius:50%; background:var(--accent); }}

        .pipe-flow {{ display:flex; align-items:center; justify-content:center; gap:0.35rem; flex-wrap:wrap; margin:0.8rem 0; }}
        .pipe-node {{ padding:0.4rem 0.7rem; border-radius:4px; font-size:0.75rem; font-weight:600; text-align:center; background:var(--surface-2); color:var(--text-2); border:1px solid var(--border-lt); }}
        .pipe-node.n-accent {{ background:var(--accent-light); color:var(--accent); border-color:rgba(66,99,235,0.2); }}
        .pipe-arrow {{ color:var(--text-cap); font-size:0.9rem; }}

        .divider {{ width:40px; height:2px; background:var(--border); margin:2rem auto; }}
        .footer {{ margin-top:2rem; padding:1rem 0; border-top:1px solid var(--border-lt); text-align:center; }}
        .footer p {{ font-size:0.75rem; color:var(--text-cap); margin:0; }}
    </style>
</head>
<body>

<div class="top-bar">
    <div class="top-bar-logo">FASTA Analyzer</div>
    <div class="top-bar-info">analyze_fasta.py &middot; Biopython</div>
</div>
<div class="accent-strip"></div>

<div class="container">

    <header class="hero">
        <h1>FASTA Analysis Report</h1>
        <div><span class="hero-file">{html_mod.escape(report['file'])}</span></div>
        <div class="hero-meta">Generated on {report.get('analysis_date', '')} &middot; analyze_fasta.py</div>
    </header>

    <section>
        <div class="sn">1. Summary</div>
        <h2>File statistics</h2>
        <table class="summary-table"><tbody>{summary_rows}</tbody></table>
    </section>

    <section>
        <div class="sn">2. Per-sequence analysis</div>
        <h2>Sequence details</h2>
        {seq_cards}
        {overflow_note}
    </section>

    <div class="divider"></div>

    <section>
        <div class="sn">3. Raw data</div>
        <h2>JSON output</h2>
        <p>Complete analysis output in JSON format, as consumed by Claude for biological interpretation.</p>
        <button class="json-toggle" onclick="var b=document.getElementById('json-block');b.style.display=b.style.display==='block'?'none':'block';this.textContent=b.style.display==='block'?'Hide JSON':'Show JSON'">Show JSON</button>
        <div class="json-block" id="json-block"><pre>{json_escaped}</pre></div>
    </section>

    <div class="divider"></div>

    <section>
        <div class="sn">4. How it works</div>
        <div class="how-section">
            <h2>Skill: analyze-fasta</h2>
            <p>This tool combines a <strong>Python script</strong> that performs bioinformatics analysis with <strong>Claude (Anthropic AI)</strong> that interprets results in biological context. The script computes; Claude interprets.</p>

            <div class="pipe-flow">
                <div class="pipe-node">FASTA file</div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-accent">analyze_fasta.py<br><small>Biopython</small></div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node">JSON + HTML</div>
                <div class="pipe-arrow">&rarr;</div>
                <div class="pipe-node n-accent">Claude<br><small>Interpretation</small></div>
            </div>

            <div class="how-grid">
                <div class="how-card">
                    <h3>Libraries</h3>
                    <ul>
                        <li><strong>Biopython</strong> (<code>Bio.SeqIO</code>) &mdash; FASTA parsing</li>
                        <li><strong>Bio.SeqUtils</strong> &mdash; GC content, molecular weight</li>
                        <li><strong>Bio.SeqUtils.ProtParam</strong> &mdash; protein physicochemical analysis (pI, GRAVY, stability, secondary structure)</li>
                        <li><strong>NumPy</strong> &mdash; numerical computations</li>
                        <li><strong>Python stdlib</strong> &mdash; Counter, json, argparse, pathlib</li>
                    </ul>
                </div>
                <div class="how-card">
                    <h3>What it analyzes</h3>
                    <ul>
                        <li><strong>Auto-detection</strong> of sequence type (nucleotide vs protein)</li>
                        <li><strong>Nucleotides:</strong> GC%, base composition, dinucleotides, ORFs in 3 frames, molecular weight, N50</li>
                        <li><strong>Proteins:</strong> MW, pI, instability index, GRAVY, aromaticity, secondary structure, charged residues</li>
                        <li><strong>Multi-sequence:</strong> aggregate statistics, N50 for assembly quality</li>
                    </ul>
                </div>
                <div class="how-card">
                    <h3>ORF detection</h3>
                    <p>Searches 3 forward reading frames (+1, +2, +3) for ORFs starting with <code>ATG</code> and ending at stop codons (<code>TAA</code>, <code>TAG</code>, <code>TGA</code>), minimum 100 amino acids (300 bp).</p>
                    <p><strong>Limitation:</strong> does not search the reverse complement or model introns. For complete gene prediction, use Augustus, GenScan, or Prodigal.</p>
                </div>
                <div class="how-card">
                    <h3>Claude's role</h3>
                    <p>The <code>analyze-fasta</code> skill instructs Claude to:</p>
                    <ul>
                        <li>Run the script with <code>--json</code></li>
                        <li>Interpret GC%, ORFs, pI, GRAVY in biological context</li>
                        <li>Infer organism, function, or sequence type</li>
                        <li>Compare against reference values</li>
                        <li>Suggest follow-up analyses (BLAST, phylogeny, etc.)</li>
                    </ul>
                </div>
            </div>

            <h4 style="margin-top:1.2rem">Usage</h4>
            <table class="summary-table" style="margin-top:0.4rem">
                <tr><td><code>python3 analyze_fasta.py input.fasta</code></td><td>Text report to terminal</td></tr>
                <tr><td><code>python3 analyze_fasta.py input.fasta --json</code></td><td>JSON output for pipelines</td></tr>
                <tr><td><code>python3 analyze_fasta.py input.fasta --html</code></td><td>Generate this HTML report</td></tr>
            </table>
        </div>
    </section>

    <footer class="footer">
        <p>Generated by analyze_fasta.py &middot; Biopython {get_biopython_version()} &middot; Python {sys.version.split()[0]}</p>
    </footer>
</div>

</body>
</html>"""


def get_biopython_version():
    try:
        import Bio
        return Bio.__version__
    except Exception:
        return "?"


# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="General-purpose FASTA file analyzer")
    parser.add_argument("fasta", help="FASTA file to analyze")
    parser.add_argument("--json", action="store_true", help="Output in JSON format")
    parser.add_argument("--html", nargs="?", const="auto", default=None, metavar="OUTPUT",
                        help="Generate HTML report. No argument: auto-generates in output/output_<name>.html")
    args = parser.parse_args()

    report = analyze_fasta(args.fasta)

    if "error" in report:
        print(f"ERROR: {report['error']}", file=sys.stderr)
        sys.exit(1)

    json_str = json.dumps(report, indent=2, ensure_ascii=False)

    if args.html is not None:
        if args.html == "auto":
            fasta_stem = Path(args.fasta).stem
            project_root = Path(args.fasta).resolve().parent
            output_dir = project_root / "output"
            output_dir.mkdir(exist_ok=True)
            output_path = output_dir / f"output_{fasta_stem}.html"
        else:
            output_path = Path(args.html)
        html_content = generate_html_report(report, json_str)
        output_path.write_text(html_content, encoding="utf-8")
        print(f"HTML report generated: {output_path}")
    elif args.json:
        print(json_str)
    else:
        print(format_text_report(report))


if __name__ == "__main__":
    main()
