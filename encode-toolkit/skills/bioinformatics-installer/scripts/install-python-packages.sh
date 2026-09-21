#!/usr/bin/env bash
# Install Python packages for ENCODE data analysis
# Usage: bash install-python-packages.sh [--all | --singlecell | --hic | --deeptools | --genomics]
#
# Categories:
#   --all          Install all packages (default if no argument)
#   --singlecell   scanpy, scvi-tools, harmony-pytorch, scrublet, scanorama, bbknn
#   --hic          cooler, cooltools, hic-straw, pyGenomeTracks
#   --deeptools    deeptools, pyBigWig, pysam, pybedtools
#   --genomics     Core genomics (numpy, pandas, scipy, matplotlib, seaborn)
#
# Every install is constrained by constraints.txt, a lock file with exact versions for the
# full dependency tree (Python 3.10+), so the same command gives the same environment.
# The direct dependencies live in requirements.in. To refresh the lock:
#   uv pip compile --universal --python-version 3.10 requirements.in -o constraints.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONSTRAINTS="$SCRIPT_DIR/constraints.txt"

if [ ! -f "$CONSTRAINTS" ]; then
    echo "ERROR: $CONSTRAINTS not found. It must sit next to this script."
    exit 1
fi

pip_install() {
    pip3 install --constraint "$CONSTRAINTS" "$@"
}

echo "============================================"
echo "ENCODE Bioinformatics Python Package Installer"
echo "============================================"
echo ""
echo "Python: $(python3 --version 2>&1)"
echo "pip:    $(pip3 --version 2>&1 | head -1)"
echo ""

CATEGORY="${1:---all}"

install_genomics() {
    echo "--- Installing Core Genomics packages ---"
    pip_install \
        numpy \
        pandas \
        scipy \
        matplotlib \
        seaborn \
        scikit-learn \
        statsmodels \
        h5py \
        loompy
}

install_singlecell() {
    echo "--- Installing Single-Cell packages ---"
    pip_install \
        scanpy \
        anndata \
        scvi-tools \
        scrublet \
        scanorama \
        bbknn \
        harmony-pytorch \
        leidenalg \
        louvain \
        umap-learn
    echo ""
    echo "NOTE: CellBender requires separate install (GPU recommended):"
    echo "  pip install cellbender"
}

install_hic() {
    echo "--- Installing Hi-C Analysis packages ---"
    pip_install \
        cooler \
        cooltools \
        hic-straw \
        pyGenomeTracks \
        bioframe
}

install_deeptools() {
    echo "--- Installing Genomics / Signal Processing packages ---"
    pip_install \
        deeptools \
        pyBigWig \
        pysam \
        pybedtools
}

case "$CATEGORY" in
    --all)
        install_genomics
        install_singlecell
        install_hic
        install_deeptools
        ;;
    --singlecell)
        install_genomics
        install_singlecell
        ;;
    --hic)
        install_genomics
        install_hic
        ;;
    --deeptools)
        install_genomics
        install_deeptools
        ;;
    --genomics)
        install_genomics
        ;;
    *)
        echo "Unknown category: $CATEGORY"
        echo "Usage: bash install-python-packages.sh [--all | --singlecell | --hic | --deeptools | --genomics]"
        exit 1
        ;;
esac

echo ""
echo "============================================"
echo "Python package installation complete."
echo "============================================"
