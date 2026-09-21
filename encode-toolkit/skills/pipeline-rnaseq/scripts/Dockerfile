# ENCODE RNA-seq Pipeline Docker Image
# STAR 2.7.11b, RSEM 1.3.3, kallisto 0.50.1, RSeQC 5.0.3, samtools 1.19

# Several tools below are distributed only as linux/x86_64 binaries, so the image is
# pinned to that platform (it runs under emulation on arm64 hosts).
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
# perl is listed explicitly: FastQC, Trim Galore and the RSEM wrappers are Perl programs
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl git ca-certificates perl unzip bc \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-openssl-dev \
    python3 python3-pip python3-dev openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# STAR 2.7.11b (statically linked release binary)
RUN wget -q https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz && \
    tar xzf 2.7.11b.tar.gz && cp STAR-2.7.11b/bin/Linux_x86_64_static/STAR /usr/local/bin/ && \
    rm -rf STAR-2.7.11b 2.7.11b.tar.gz

# RSEM 1.3.3
RUN wget -q https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.3.tar.gz && \
    tar xzf v1.3.3.tar.gz && cd RSEM-1.3.3 && make -j4 && make install && \
    cd .. && rm -rf RSEM-1.3.3 v1.3.3.tar.gz

# kallisto 0.50.1 (its index format differs from 0.48 and earlier: build the index with 0.50.1)
RUN wget -q https://github.com/pachterlab/kallisto/releases/download/v0.50.1/kallisto_linux-v0.50.1.tar.gz && \
    tar xzf kallisto_linux-v0.50.1.tar.gz && cp kallisto/kallisto /usr/local/bin/ && \
    rm -rf kallisto kallisto_linux-v0.50.1.tar.gz

# samtools 1.19
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar xjf samtools-1.19.tar.bz2 && cd samtools-1.19 && \
    ./configure --prefix=/usr/local && make -j4 && make install && \
    cd .. && rm -rf samtools-1.19*

# UCSC Kent tools
RUN wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
    -O /usr/local/bin/bedGraphToBigWig && chmod +x /usr/local/bin/bedGraphToBigWig

# Python tools: RSeQC, MultiQC, cutadapt (used by Trim Galore)
RUN pip3 install --no-cache-dir \
    RSeQC==5.0.3 multiqc==1.21 cutadapt==4.6

# Trim Galore 0.6.10 (a Perl script; not published on PyPI)
RUN wget -q https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz && \
    tar xzf 0.6.10.tar.gz && ln -s /TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    rm 0.6.10.tar.gz

# FastQC 0.12.1
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip -d /opt/ && chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/ && rm fastqc_v0.12.1.zip

WORKDIR /data
CMD ["/bin/bash"]
