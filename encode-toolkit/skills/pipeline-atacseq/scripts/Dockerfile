# ENCODE ATAC-seq Pipeline Docker Image
# Bowtie2 2.5.4, samtools 1.19, Picard 3.1.1, bedtools 2.31.0, MACS2 2.2.9.1, IDR 2.0.4.2, deepTools 3.5.5

# Several tools below are distributed only as linux/x86_64 binaries, so the image is
# pinned to that platform (it runs under emulation on arm64 hosts).
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
# perl is listed explicitly: FastQC and Trim Galore are Perl programs
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl git ca-certificates perl unzip bc \
    zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-openssl-dev \
    python3 python3-pip python3-dev openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Python tools. numpy is pinned on the same line so that pip reports a conflict rather than
# upgrading it past 1.23, the last release IDR 2.0.4.2 works with.
RUN pip3 install --no-cache-dir \
    numpy==1.23.5 macs2==2.2.9.1 deeptools==3.5.5 multiqc==1.21 cutadapt==4.6

# IDR 2.0.4.2 is not on PyPI (the PyPI package named idr is unrelated). Install the release
# ENCODE uses from its GitHub tag. It imports numpy in setup.py and uses numpy.int at
# runtime, so numpy is installed first and held below 1.24.
RUN pip3 install --no-cache-dir --no-build-isolation \
    https://github.com/kundajelab/idr/archive/refs/tags/2.0.4.2.tar.gz

# Bowtie2 2.5.4
RUN wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip && \
    unzip -q bowtie2-2.5.4-linux-x86_64.zip -d /opt/ && \
    ln -s /opt/bowtie2-2.5.4-linux-x86_64/bowtie2* /usr/local/bin/ && \
    rm bowtie2-2.5.4-linux-x86_64.zip

# samtools 1.19
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar xjf samtools-1.19.tar.bz2 && cd samtools-1.19 && \
    ./configure --prefix=/usr/local && make -j4 && make install && \
    cd .. && rm -rf samtools-1.19*

# bedtools 2.31.0
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools-2.31.0.tar.gz && \
    tar xzf bedtools-2.31.0.tar.gz && cd bedtools2 && \
    make -j4 && cp bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedtools2 bedtools-2.31.0.tar.gz

# Picard 3.1.1 (needs Java 17)
RUN wget -q https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar \
    -O /opt/picard.jar && \
    echo '#!/bin/bash\njava -jar /opt/picard.jar "$@"' > /usr/local/bin/picard && \
    chmod +x /usr/local/bin/picard

# Trim Galore 0.6.10 (a Perl script; not published on PyPI)
RUN wget -q https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz && \
    tar xzf 0.6.10.tar.gz && ln -s /TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    rm 0.6.10.tar.gz

# FastQC 0.12.1
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip -d /opt/ && chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/ && rm fastqc_v0.12.1.zip

# UCSC Kent tools
RUN wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
    -O /usr/local/bin/bedGraphToBigWig && chmod +x /usr/local/bin/bedGraphToBigWig

WORKDIR /data
CMD ["/bin/bash"]
