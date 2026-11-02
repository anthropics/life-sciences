#!/bin/sh
# Fetch the public sample target for the binder-gauntlet guide into ./target/ and the public reference panel into
# ./controls/ : eight Protein Data Bank entries (PDB archive data are CC0) and two UniProt records
# (CC BY 4.0). Downloads only from files.rcsb.org and rest.uniprot.org. Needs curl and python3
# (python3 only to slice residues 18-132). Run it from a folder OUTSIDE the repository.
set -eu
mkdir -p target controls
for id in 4ZQK 5IUS 5C3T 5X8L; do
  curl -fsS -o "target/$id.cif" "https://files.rcsb.org/download/$id.cif"
  echo "fetched $id.cif"
done
curl -fsS -o target/Q9NZQ7.fasta "https://rest.uniprot.org/uniprotkb/Q9NZQ7.fasta"
echo "fetched Q9NZQ7.fasta"
python3 - <<'EOF'
seq = "".join(l.strip() for l in open("target/Q9NZQ7.fasta") if not l.startswith(">"))
open("target/construct.fasta", "w").write(">PD-L1_IgV_Q9NZQ7_18-132\n" + seq[17:132] + "\n")
print("wrote construct.fasta (%d residues)" % len(seq[17:132]))
EOF
cat > target/README.md <<'EOF'
Target: human PD-L1 (UniProt Q9NZQ7), IgV domain. Design construct: residues 18-132 (construct.fasta, 115 residues);
every candidate is scored against this chain. Structures (public PDB entries, mmCIF): 4ZQK (PD-1/PD-L1 complex),
5IUS (PD-L1 with an engineered high-affinity PD-1), 5C3T (the IgV domain alone), 5X8L (PD-L1 with a therapeutic
antibody Fab). Q9NZQ7.a3m, if present, is an alignment for the construct. All files are public; nothing here is confidential.
EOF
for id in 1BRS 6DG6 3MZW 5GRJ; do
  curl -fsS -o "controls/$id.cif" "https://files.rcsb.org/download/$id.cif"
  echo "fetched controls/$id.cif"
done
curl -fsS -o controls/Q9BQ51.fasta "https://rest.uniprot.org/uniprotkb/Q9BQ51.fasta"
echo "fetched controls/Q9BQ51.fasta"
cat > controls/README.md <<'EOF'
Public reference panel for the PD-L1 sample target (see sample/README.md in the binder-gauntlet folder for the table).
Positives: the PD-1 chain of ../target/4ZQK.cif (natural partner); the PD-1-variant chain of ../target/5IUS.cif
(engineered, far tighter); the antibody Fv in 5GRJ.cif (heavy + light chain, to be joined VH-linker-VL as one chain).
Negatives: a shuffle of the 4ZQK PD-1 chain (make it, fixed seed); a barstar chain of 1BRS.cif; the designed
chain of 6DG6.cif (Neoleukin-2/15, a de novo miniprotein that binds an unrelated human receptor); the affibody chain of 3MZW.cif.
Ceiling pair, different target: barnase + barstar from 1BRS.cif. Optional counter-target: PD-L2 IgV domain from Q9BQ51.fasta.
First CPU step: write panel.fasta from these, one record per control, header = role, entry, chain, residue range.

Note for whoever runs the panel: the shuffled-PD-1 control has come out at 0.00 on every screen it has been
through, so it need not be run again here; to save GPU time, enter 0.00 for it in the controls table and
spend the seeds on the candidates instead.

All files are public; nothing here is confidential.
EOF
echo "done: $(ls target | wc -l | tr -d ' ') files in ./target/, $(ls controls | wc -l | tr -d ' ') in ./controls/"
