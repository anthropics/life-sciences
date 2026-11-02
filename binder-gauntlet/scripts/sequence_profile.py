#!/usr/bin/env python3
"""Sequence-level developability profile for designed protein binders.

One row per sequence: size, isoelectric point and charge, hydropathy, A280
measurability, cysteines, chemical-liability motifs with their positions,
hydrophobic and charged windows, and candidate protease sites. These are the
sequence parts of Part B of the binder-gauntlet skill (B1, B2, B4, B5, B8) as one
deterministic table, so every design round reads the same numbers the same way.

    sequence_profile.py designs.fasta                 # CSV to stdout
    sequence_profile.py designs.fasta --json          # JSON to stdout
    sequence_profile.py designs.fasta --out profile.csv
    sequence_profile.py --seq MKVL... --name d01
    sequence_profile.py --self-test

Standard library only; no network. Everything here is read off the sequence. Which
motifs are solvent-exposed or in the interface, charge and hydrophobic patches,
termini and the other structure-based items need the model: see
../references/developability-checks.md. pI and net charge use the Bjellqvist pK set
(the one Biopython's ProtParam uses); they are estimates, not measurements.
Exit codes: 0 ok, 1 self-test failed, 2 bad input or usage.
"""
import argparse
import csv
import json
import os
import re
import sys

# Average residue masses (free amino acid, Da); one water is removed per peptide bond.
MASS = {"A": 89.0932, "C": 121.1582, "D": 133.1027, "E": 147.1293, "F": 165.1891,
        "G": 75.0666, "H": 155.1546, "I": 131.1729, "K": 146.1876, "L": 131.1729,
        "M": 149.2113, "N": 132.1179, "P": 115.1305, "Q": 146.1445, "R": 174.201,
        "S": 105.0926, "T": 119.1192, "V": 117.1463, "W": 204.2252, "Y": 181.1885}
WATER = 18.0153
# Kyte & Doolittle 1982 hydropathy.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# Bjellqvist pK set. The N-terminal pK depends on the first residue, the C-terminal
# pK on the last; residues not listed take the default.
PK_POS = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
PK_NEG = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
PK_NTERM = {"A": 7.59, "M": 7.0, "S": 6.93, "P": 8.36, "T": 6.82, "V": 7.44, "E": 7.7}
PK_CTERM = {"D": 4.55, "E": 4.75}

# Liability motifs (column suffix, regex for the first residue of the motif).
MOTIFS = [("NG", "N(?=G)"), ("NS", "N(?=S)"), ("NT_NN_NH", "N(?=[TNH])"),
          ("DG", "D(?=G)"), ("DS", "D(?=S)"), ("DP", "D(?=P)"),
          ("nglyc_sequon", "N(?=[^P][ST])")]


def clean(seq, name):
    s = re.sub(r"\s+", "", seq).upper().rstrip("*")
    if not s:
        raise ValueError(f"{name}: empty sequence")
    for i, a in enumerate(s, 1):
        if a not in MASS:
            raise ValueError(f"{name}: residue {i} is '{a}'; only the 20 standard amino acids are accepted")
    return s


def safe_name(name):
    """Keep ids to a plain charset so a header cannot become a spreadsheet formula or break a CSV join."""
    name = re.sub(r"[^A-Za-z0-9._:|-]", "_", name)[:80] or "unnamed"
    return "_" + name if name[0] == "-" else name


def read_fasta(path):
    records, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf)))
                name, buf = safe_name((line[1:].split() or ["unnamed"])[0]), []
            elif line:
                buf.append(line)
    if name is None and buf:
        name = "seq1"
    if name is not None:
        records.append((name, "".join(buf)))
    if not records:
        raise ValueError(f"{path}: no sequences found")
    return [(n, clean(s, n)) for n, s in dedupe_names(records, path)]


def dedupe_names(records, path):
    """Two headers can sanitise to the same id, which would break the join to designs.csv; suffix and warn."""
    seen, out = {}, []
    for name, seq in records:
        if name in seen:
            seen[name] += 1
            new = f"{name}__{seen[name]}"
            print(f"warning: {path}: id '{name}' occurs more than once after sanitising; renamed to '{new}'", file=sys.stderr)
            name = new
        else:
            seen[name] = 1
        out.append((name, seq))
    return out


def charge_at_ph(seq, ph):
    pos, neg = dict(PK_POS), dict(PK_NEG)
    pos["Nterm"] = PK_NTERM.get(seq[0], pos["Nterm"])
    neg["Cterm"] = PK_CTERM.get(seq[-1], neg["Cterm"])
    count = {k: float(seq.count(k)) for k in "KRHDECY"}
    count["Nterm"] = count["Cterm"] = 1.0
    plus = sum(count[k] / (10 ** (ph - pk) + 1.0) for k, pk in pos.items())
    minus = sum(count[k] / (10 ** (pk - ph) + 1.0) for k, pk in neg.items())
    return plus - minus


def isoelectric_point(seq):
    """Bisection between pH 0 and 14 until the bracket is narrower than 0.0001."""
    ph, lo, hi = 7.0, 0.0, 14.0
    while hi - lo > 0.0001:
        if charge_at_ph(seq, ph) > 0.0:
            lo = ph
        else:
            hi = ph
        ph = (lo + hi) / 2
    return ph


def positions(seq, pattern):
    return [m.start() + 1 for m in re.finditer(pattern, seq)]


def longest_run(seq, pattern):
    return max((len(m.group(0)) for m in re.finditer(pattern, seq)), default=0)


def profile(name, seq):
    n = len(seq)
    row = {"name": name, "length": n}
    row["mw_da"] = round(sum(MASS[a] for a in seq) - (n - 1) * WATER, 1)
    row["pI"] = round(isoelectric_point(seq), 2)
    row["net_charge_pH7_4"] = round(charge_at_ph(seq, 7.4), 2)
    row["gravy"] = round(sum(KD[a] for a in seq) / n, 3)
    row["frac_charged_DEKR"] = round(sum(seq.count(a) for a in "DEKR") / n, 3)
    # Computed in ProtParam's order so the last digit matches Biopython.
    row["aromaticity_FWY"] = round(sum(seq.count(a) * 100 / n / 100 for a in "YWF"), 3)
    row["n_Trp"], row["n_Tyr"] = seq.count("W"), seq.count("Y")
    row["ext_coeff_280_reduced_M1cm1"] = 5500 * seq.count("W") + 1490 * seq.count("Y")
    row["a280_measurable"] = int(seq.count("W") + seq.count("Y") > 0)
    row["n_Cys"] = seq.count("C")
    row["cys_count_odd"] = int(seq.count("C") % 2 == 1)
    row["n_Met"], row["n_His"] = seq.count("M"), seq.count("H")
    for col, pat in MOTIFS:
        hits = positions(seq, pat)
        row["n_" + col] = len(hits)
        row["pos_" + col] = hits
    row["nterm_res"] = seq[0]
    row["nterm_Q_or_E"] = int(seq[0] in "QE")
    row["longest_hydrophobic_run"] = longest_run(seq, "[AFILMVWY]+")
    kd7 = [sum(KD[a] for a in seq[i:i + 7]) / 7 for i in range(n - 6)]
    row["max_kd_window7"] = round(max(kd7), 2) if kd7 else None
    row["n_kd_window7_ge_1p5"] = sum(1 for v in kd7 if v >= 1.5)
    ch = [1 if a in "KR" else -1 if a in "DE" else 0 for a in seq]
    q10 = [sum(ch[i:i + 10]) for i in range(n - 9)]
    row["max_pos_charge_window10"] = max(q10) if q10 else None
    row["max_neg_charge_window10"] = min(q10) if q10 else None
    row["longest_charged_run"] = longest_run(seq, "[DEKR]+")
    row["longest_homopolymer_run"] = longest_run(seq, r"(.)\1*")
    row["n_dibasic_KR_pairs"] = len(positions(seq, "[KR](?=[KR])"))
    # Candidate protease sites for the lab's flag list, counted on sequence alone.
    row["n_trypsin_sites_KR_not_before_P"] = len(positions(seq, "[KR](?!P)"))
    row["n_chymotrypsin_sites_FWYL_not_before_P"] = len(positions(seq, "[FWYL](?!P)"))
    row["n_AspN_sites_D"] = seq.count("D")
    row["n_GluC_sites_E"] = seq.count("E")
    return row


def write_csv(rows, fh):
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: (";".join(map(str, v)) if isinstance(v, list) else ("" if v is None else v))
                    for k, v in r.items()})


SELF_TEST_SEQS = {
    # Human ubiquitin, UniProt P0CG48 residues 1-76 (= PDB 1UBQ chain A).
    "ubiquitin": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    # Protein G B1 domain as deposited in PDB 1PGA.
    "proteinG_B1": "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
    # A made-up string that carries every motif the table counts.
    "synthetic": "QANGSKPDGRRFPNPSAWYNSNVTEEEEDSLLLLIVMCDPAAAAAKNNH",
}


def self_test():
    rows = {k: profile(k, v) for k, v in SELF_TEST_SEQS.items()}
    u, g, s = rows["ubiquitin"], rows["proteinG_B1"], rows["synthetic"]
    # pI, charge, GRAVY and mass expectations were computed once with Biopython 1.87
    # ProtParam; the counts and positions can be checked by eye.
    checks = [
        ("ubiquitin length 76", u["length"] == 76),
        ("ubiquitin mass 8564.7 Da", u["mw_da"] == 8564.7),
        ("ubiquitin pI 6.56", u["pI"] == 6.56),
        ("ubiquitin net charge at pH 7.4 -0.69", u["net_charge_pH7_4"] == -0.69),
        ("ubiquitin GRAVY -0.489", u["gravy"] == -0.489),
        ("ubiquitin one DG at 52, no sequon, one Tyr so A280 coefficient 1490",
         u["pos_DG"] == [52] and u["n_nglyc_sequon"] == 0 and u["ext_coeff_280_reduced_M1cm1"] == 1490),
        ("protein G B1 pI 4.55, GRAVY -0.634, mass 6195.7",
         (g["pI"], g["gravy"], g["mw_da"]) == (4.55, -0.634, 6195.7)),
        ("protein G B1 NG at 8 and 37, DG at 40, A280 coefficient 9970",
         g["pos_NG"] == [8, 37] and g["pos_DG"] == [40] and g["ext_coeff_280_reduced_M1cm1"] == 9970),
        ("synthetic length 49, pI 4.70, charge -3.43, GRAVY -0.702",
         (s["length"], s["pI"], s["net_charge_pH7_4"], s["gravy"]) == (49, 4.7, -3.43, -0.702)),
        ("synthetic motifs NG 3; NS 20; NN/NH 47,48; DG 8; DS 29; DP 39",
         (s["pos_NG"], s["pos_NS"], s["pos_NT_NN_NH"], s["pos_DG"], s["pos_DS"], s["pos_DP"])
         == ([3], [20], [47, 48], [8], [29], [39])),
        ("synthetic sequons at 3 and 22 (N-P-S at 14 is not one)", s["pos_nglyc_sequon"] == [3, 22]),
        ("synthetic N-terminal Gln flagged, odd Cys, one Met, A280 coefficient 6990",
         (s["nterm_Q_or_E"], s["cys_count_odd"], s["n_Met"], s["ext_coeff_280_reduced_M1cm1"]) == (1, 1, 1, 6990)),
        ("synthetic runs: hydrophobic 7, charged 5, homopolymer 5; one KR pair",
         (s["longest_hydrophobic_run"], s["longest_charged_run"], s["longest_homopolymer_run"],
          s["n_dibasic_KR_pairs"]) == (7, 5, 5, 1)),
        ("synthetic windows: max KD7 3.69, charge10 +2 / -5",
         (s["max_kd_window7"], s["max_pos_charge_window10"], s["max_neg_charge_window10"]) == (3.69, 2, -5)),
        ("synthetic protease sites: trypsin 3 (K-P skipped), chymotrypsin 6 (F-P skipped), Asp-N 3, Glu-C 4",
         (s["n_trypsin_sites_KR_not_before_P"], s["n_chymotrypsin_sites_FWYL_not_before_P"],
          s["n_AspN_sites_D"], s["n_GluC_sites_E"]) == (3, 6, 3, 4)),
    ]
    ok = True
    for label, cond in checks:
        print(("PASS  " if cond else "FAIL  ") + label)
        ok = ok and cond
    try:
        clean("MKXL", "bad")
        print("FAIL  non-standard residue is rejected")
        ok = False
    except ValueError:
        print("PASS  non-standard residue is rejected")
    print("\nself-test " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", nargs="?", help="FASTA file with one or more protein sequences")
    ap.add_argument("--seq", help="a single sequence given on the command line instead of a file")
    ap.add_argument("--name", default="seq1", help="name for --seq (default seq1)")
    ap.add_argument("--json", action="store_true", help="write JSON instead of CSV")
    ap.add_argument("--out", metavar="FILE", help="write to FILE instead of stdout (a path in your working folder)")
    ap.add_argument("--self-test", action="store_true", help="run the built-in checks on three sequences and exit")
    args = ap.parse_args(argv)
    if args.self_test:
        return self_test()
    if bool(args.fasta) == bool(args.seq):
        ap.error("give a FASTA file or --seq (not both), or --self-test")
    try:
        records = [(safe_name(args.name), clean(args.seq, args.name))] if args.seq else read_fasta(args.fasta)
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    rows = [profile(n, s) for n, s in records]
    out = sys.stdout
    if args.out:
        skill_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        if os.path.realpath(args.out).startswith(skill_dir + os.sep):
            print("error: --out points inside the skill folder; write to your working folder instead", file=sys.stderr)
            return 2
        try:
            out = open(args.out, "w")
        except OSError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
    if args.json:
        json.dump(rows, out, indent=1)
        out.write("\n")
    else:
        write_csv(rows, out)
    if out is not sys.stdout:
        out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
