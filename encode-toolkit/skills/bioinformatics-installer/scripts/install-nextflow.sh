#!/usr/bin/env bash
# Install the pinned Nextflow release and check for a container runtime.
# Docker and Singularity are only checked: when one is missing, the script prints how to get it.
# Usage: bash install-nextflow.sh [--docker | --singularity | --both]
#
# Options:
#   --docker       Nextflow, then check for Docker (default, for local/cloud)
#   --singularity  Nextflow, then check for Singularity/Apptainer (for HPC clusters)
#   --both         Nextflow, then check for both

set -euo pipefail

# Pinned Nextflow release. The ENCODE Toolkit pipelines are validated against this version
# with `nextflow lint` and `nextflow run -preview`. To move to another release, change the
# version and replace the checksum with the sha256 published for the
# `nextflow-<version>-dist` asset on https://github.com/nextflow-io/nextflow/releases
NEXTFLOW_VERSION="26.04.6"
NEXTFLOW_SHA256="182a63c74074e2dc7956ffa3c8cd59de952ed2c44394e21faf5e1736b945444c"
NEXTFLOW_URL="https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow-${NEXTFLOW_VERSION}-dist"

sha256_of() {
    if command -v sha256sum &> /dev/null; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

echo "============================================"
echo "ENCODE Nextflow Pipeline Infrastructure Setup"
echo "============================================"
echo ""

MODE="${1:---docker}"

# --- Install Nextflow ---
install_nextflow() {
    echo "--- Installing Nextflow ---"
    # Only an existing install of exactly the pinned release is accepted. Otherwise the pinned,
    # checksum-verified release is installed; a launcher it would overwrite is kept as a backup.
    local existing_version=""
    if command -v nextflow &> /dev/null; then
        # A launcher that cannot start (no Java, broken install) counts as "no usable version"
        existing_version="$(nextflow -version 2>&1 | awk '$1 == "version" {print $2; exit}' || true)"
    fi
    if [ "$existing_version" = "$NEXTFLOW_VERSION" ]; then
        echo "Nextflow ${NEXTFLOW_VERSION} already installed: $(command -v nextflow)"
    else
        if [ -n "$existing_version" ]; then
            echo "Found Nextflow ${existing_version} at $(command -v nextflow); the pipelines are validated against ${NEXTFLOW_VERSION}."
        fi
        # Check Java
        if ! command -v java &> /dev/null; then
            echo "ERROR: Java 17+ is required for Nextflow."
            echo "Install Java first:"
            echo "  macOS:  brew install openjdk@17"
            echo "  Ubuntu: sudo apt-get install -y openjdk-17-jdk"
            echo "  conda:  conda install -c conda-forge 'openjdk>=17'"
            exit 1
        fi

        JAVA_VER=$(java -version 2>&1 | head -1 | awk -F '"' '{print $2}' | awk -F '.' '{print $1}')
        if [ "$JAVA_VER" -lt 17 ] 2>/dev/null; then
            echo "WARNING: Java $JAVA_VER detected. Nextflow ${NEXTFLOW_VERSION} requires Java 17+."
        fi

        # Download the pinned, self-contained release and verify it before it is ever executed
        local tmp_file
        tmp_file="$(mktemp)"
        echo "Downloading Nextflow ${NEXTFLOW_VERSION}..."
        curl -fsSL "$NEXTFLOW_URL" -o "$tmp_file"

        local actual_sha
        actual_sha="$(sha256_of "$tmp_file")"
        if [ "$actual_sha" != "$NEXTFLOW_SHA256" ]; then
            echo "ERROR: checksum mismatch for Nextflow ${NEXTFLOW_VERSION}; refusing to install."
            echo "  expected: $NEXTFLOW_SHA256"
            echo "  actual:   $actual_sha"
            rm -f "$tmp_file"
            exit 1
        fi
        echo "Checksum verified."

        # Install into a directory on PATH when possible, otherwise into ~/.local/bin
        local install_dir
        if [ -w /usr/local/bin ]; then
            install_dir="/usr/local/bin"
        else
            install_dir="$HOME/.local/bin"
            mkdir -p "$install_dir"
        fi
        local nextflow_bin="$install_dir/nextflow"
        if [ -e "$nextflow_bin" ]; then
            mv "$nextflow_bin" "$nextflow_bin.previous"
            echo "Kept the launcher that was there as $nextflow_bin.previous"
        fi
        mv "$tmp_file" "$nextflow_bin"
        chmod 755 "$nextflow_bin"
        echo "Nextflow installed to $nextflow_bin"

        # Put the pinned release ahead of any other Nextflow on PATH
        if [ "$(command -v nextflow || true)" != "$nextflow_bin" ]; then
            echo "Put it first on your PATH: export PATH=\"$install_dir:\$PATH\""
        fi

        # Call the installed file directly: it may not be on PATH yet in this shell
        "$nextflow_bin" -version
    fi
    echo ""
}

# --- Install Docker ---
install_docker() {
    echo "--- Checking Docker ---"
    if command -v docker &> /dev/null; then
        echo "Docker already installed: $(docker --version)"
    else
        OS="$(uname -s)"
        case "$OS" in
            Darwin)
                echo "macOS detected. Install Docker Desktop:"
                echo "  brew install --cask docker"
                echo "  OR download from https://www.docker.com/products/docker-desktop/"
                ;;
            Linux)
                echo "Linux detected. Install Docker Engine:"
                echo "  curl -fsSL https://get.docker.com | sh"
                echo "  sudo usermod -aG docker \$USER"
                echo "  newgrp docker"
                ;;
            *)
                echo "Unsupported OS: $OS"
                echo "See https://docs.docker.com/get-docker/"
                ;;
        esac
    fi
    echo ""
}

# --- Install Singularity ---
install_singularity() {
    echo "--- Checking Singularity ---"
    if command -v singularity &> /dev/null; then
        echo "Singularity already installed: $(singularity --version)"
    elif command -v apptainer &> /dev/null; then
        echo "Apptainer (Singularity successor) already installed: $(apptainer --version)"
    else
        echo "Singularity/Apptainer not found."
        echo "For HPC clusters, check: module avail singularity"
        echo ""
        echo "Install options:"
        echo "  conda: conda install -c conda-forge singularity"
        echo "  Linux: See https://apptainer.org/docs/admin/main/installation.html"
        echo ""
        echo "NOTE: Singularity requires Linux. On macOS, use Docker or a Linux VM."
    fi
    echo ""
}

# --- Execute ---
install_nextflow

case "$MODE" in
    --docker)
        install_docker
        ;;
    --singularity)
        install_singularity
        ;;
    --both)
        install_docker
        install_singularity
        ;;
    *)
        echo "Unknown mode: $MODE"
        echo "Usage: bash install-nextflow.sh [--docker | --singularity | --both]"
        exit 1
        ;;
esac

echo "============================================"
echo "Pipeline infrastructure setup complete."
echo ""
echo "Test with:"
echo "  nextflow run hello"
echo ""
echo "Run ENCODE pipelines with:"
echo "  nextflow run main.nf -profile local    # Docker"
echo "  nextflow run main.nf -profile slurm    # Singularity + SLURM"
echo "============================================"
