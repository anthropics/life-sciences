# Several tools below are distributed only as linux/x86_64 binaries, so the image is
# pinned to that platform (it runs under emulation on arm64 hosts).
FROM --platform=linux/amd64 ubuntu:22.04

LABEL maintainer="Dr. Alex M. Mawla, PhD <https://github.com/ammawla/encode-toolkit>"
LABEL description="ENCODE Toolkit DNase-seq pipeline: BWA + Hotspot2 + HINT-ATAC"

ENV DEBIAN_FRONTEND=noninteractive

# perl is listed explicitly: FastQC and Trim Galore are Perl programs
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl gzip pigz unzip \
    python3 python3-pip python3-dev \
    openjdk-17-jre-headless \
    libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev libcurl4-openssl-dev \
    build-essential libboost-dev perl \
    bc \
    bedops \
    && rm -rf /var/lib/apt/lists/*

# samtools 1.19
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar xjf samtools-1.19.tar.bz2 && cd samtools-1.19 && \
    ./configure --prefix=/usr/local && make -j4 && make install && \
    cd .. && rm -rf samtools-1.19*

# BWA 0.7.18 (0.7.17 does not link with GCC 10 or newer)
RUN wget -q https://github.com/lh3/bwa/archive/refs/tags/v0.7.18.tar.gz && \
    tar xzf v0.7.18.tar.gz && cd bwa-0.7.18 && \
    make -j4 && cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa-0.7.18 v0.7.18.tar.gz

# bedtools 2.31.0
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools-2.31.0.tar.gz && \
    tar xzf bedtools-2.31.0.tar.gz && cd bedtools2 && \
    make -j4 && cp bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedtools2 bedtools-2.31.0.tar.gz

# UCSC Kent tools (bedGraphToBigWig)
RUN wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
    -O /usr/local/bin/bedGraphToBigWig && chmod +x /usr/local/bin/bedGraphToBigWig

# modwt 1.0 (wavelet smoothing; required by Hotspot2's density-peaks step)
RUN wget -q https://github.com/StamLab/modwt/archive/refs/tags/1.0.tar.gz && \
    tar xzf 1.0.tar.gz && make -C modwt-1.0/src && \
    cp modwt-1.0/bin/modwt /usr/local/bin/ && \
    rm -rf modwt-1.0 1.0.tar.gz

# Hotspot2 2.1.2 (latest tagged release). `make` builds the binaries into bin/; the
# wrapper scripts live in scripts/ and must sit next to each other on the PATH.
RUN wget -q https://github.com/Altius/hotspot2/archive/refs/tags/v2.1.2.tar.gz && \
    tar xzf v2.1.2.tar.gz && make -C hotspot2-2.1.2 && \
    cp hotspot2-2.1.2/bin/* hotspot2-2.1.2/scripts/* /usr/local/bin/ && \
    rm -rf hotspot2-2.1.2 v2.1.2.tar.gz

# Picard 3.1.1
RUN wget -q https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar \
    -O /opt/picard.jar && \
    echo '#!/bin/bash\njava -jar /opt/picard.jar "$@"' > /usr/local/bin/picard && \
    chmod +x /usr/local/bin/picard

# Trim Galore + cutadapt
RUN pip3 install --no-cache-dir cutadapt==4.6 && \
    wget -q https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz && \
    tar xzf 0.6.10.tar.gz && ln -s /TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    rm 0.6.10.tar.gz

# HINT (RGT). 0.13.2 depends on PyVCF, which no longer installs with current setuptools;
# 1.0.2 dropped that dependency and keeps the rgt-hint options the workflow uses.
RUN pip3 install --no-cache-dir RGT==1.0.2

# FastQC + MultiQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip -d /opt/ && chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/ && rm fastqc_v0.12.1.zip

RUN pip3 install --no-cache-dir multiqc==1.21

WORKDIR /data
