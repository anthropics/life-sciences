FROM ubuntu:22.04

LABEL maintainer="Dr. Alex M. Mawla, PhD <https://github.com/ammawla/encode-toolkit>"
LABEL description="ENCODE Toolkit Hi-C pipeline: BWA + pairtools + Juicer + cooler"

ENV DEBIAN_FRONTEND=noninteractive

# perl is listed explicitly: FastQC and Trim Galore are Perl programs
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl gzip pigz unzip \
    build-essential perl \
    python3 python3-pip python3-dev \
    openjdk-17-jre-headless \
    libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev libcurl4-openssl-dev \
    libhdf5-dev \
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

# pairtools 1.1.2 (first release line with prebuilt wheels; 1.0.3 no longer builds from source)
RUN pip3 install --no-cache-dir pairtools==1.1.2

# Juicer tools 2.20.00
RUN wget -q https://github.com/aidenlab/Juicebox/releases/download/v2.20.00/juicer_tools.2.20.00.jar \
    -O /opt/juicer_tools.jar

# cooler 0.9.3
RUN pip3 install --no-cache-dir cooler==0.9.3

# FastQC 0.12.1
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip -d /opt/ && chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/ && rm fastqc_v0.12.1.zip

# MultiQC
RUN pip3 install --no-cache-dir multiqc==1.21

# Trim Galore
RUN pip3 install --no-cache-dir cutadapt==4.6 && \
    wget -q https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz && \
    tar xzf 0.6.10.tar.gz && ln -s /TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    rm 0.6.10.tar.gz

WORKDIR /data
