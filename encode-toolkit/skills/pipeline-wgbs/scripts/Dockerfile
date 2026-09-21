# Several tools below are distributed only as linux/x86_64 binaries, so the image is
# pinned to that platform (it runs under emulation on arm64 hosts).
FROM --platform=linux/amd64 ubuntu:22.04

LABEL maintainer="Dr. Alex M. Mawla, PhD <https://github.com/ammawla/encode-toolkit>"
LABEL description="ENCODE Toolkit WGBS pipeline: Bismark + MethylDackel"

ENV DEBIAN_FRONTEND=noninteractive

# perl is listed explicitly: FastQC, Trim Galore and Bismark are Perl programs
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl unzip gzip pigz \
    build-essential perl \
    python3 python3-pip \
    openjdk-17-jre-headless \
    libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev libcurl4-openssl-dev \
    libbigwig-dev \
    && rm -rf /var/lib/apt/lists/*

# samtools 1.19
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar xjf samtools-1.19.tar.bz2 && cd samtools-1.19 && \
    ./configure --prefix=/usr/local && make -j4 && make install && \
    cd .. && rm -rf samtools-1.19*

# htslib (tabix, bgzip) 1.19
RUN wget -q https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar xjf htslib-1.19.tar.bz2 && cd htslib-1.19 && \
    ./configure --prefix=/usr/local && make -j4 && make install && ldconfig && \
    cd .. && rm -rf htslib-1.19*

# Bowtie2 2.5.4
RUN wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip && \
    unzip -q bowtie2-2.5.4-linux-x86_64.zip -d /opt/ && \
    ln -s /opt/bowtie2-2.5.4-linux-x86_64/bowtie2* /usr/local/bin/ && \
    rm bowtie2-2.5.4-linux-x86_64.zip

# Bismark 0.24.2
RUN wget -q https://github.com/FelixKrueger/Bismark/archive/v0.24.2.tar.gz && \
    tar xzf v0.24.2.tar.gz && mv Bismark-0.24.2 /opt/bismark && \
    ln -s /opt/bismark/bismark /usr/local/bin/ && \
    ln -s /opt/bismark/deduplicate_bismark /usr/local/bin/ && \
    ln -s /opt/bismark/bismark_genome_preparation /usr/local/bin/ && \
    rm v0.24.2.tar.gz

# Trim Galore 0.6.10
RUN pip3 install --no-cache-dir cutadapt==4.6 && \
    wget -q https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz && \
    tar xzf 0.6.10.tar.gz && ln -s /TrimGalore-0.6.10/trim_galore /usr/local/bin/ && \
    rm 0.6.10.tar.gz

# MethylDackel 0.6.1
RUN wget -q https://github.com/dpryan79/MethylDackel/archive/0.6.1.tar.gz && \
    tar xzf 0.6.1.tar.gz && cd MethylDackel-0.6.1 && \
    make CFLAGS="-O3 -Wall" LIBBIGWIG="-lBigWig" && cp MethylDackel /usr/local/bin/ && \
    cd .. && rm -rf MethylDackel-0.6.1*

# FastQC 0.12.1
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip -d /opt/ && chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/ && rm fastqc_v0.12.1.zip

# MultiQC
RUN pip3 install --no-cache-dir multiqc==1.21

WORKDIR /data
