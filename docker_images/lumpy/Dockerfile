# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

# Install necessary tools
RUN apt-get install -y \
    curl

# Install conda & lumpy-sv
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN conda config --add channels bioconda
RUN conda install -c bioconda lumpy-sv

# Download lumpy scripts
RUN apt-get install -y \
    wget

RUN wget -q https://github.com/arq5x/lumpy-sv/releases/download/0.3.0/lumpy-sv.tar.gz
RUN mkdir -p /miniconda/bin/scripts

RUN tar -xzf lumpy-sv.tar.gz && rm lumpy-sv.tar.gz && cp -r lumpy-sv/scripts/* /miniconda/bin/scripts && rm -rf lumpy-sv/

# Install vcftools
RUN apt-get update -y -q
RUN apt-get install -y \
    autoconf \
    build-essential \
    g++ \
    git \
    pkg-config \
    zlib1g-dev

RUN git clone https://github.com/vcftools/vcftools.git
WORKDIR vcftools
RUN ./autogen.sh
RUN ./configure --disable-dependency-tracking && make && make install

COPY b37.bed /opt/
COPY hg19.bed /opt/
COPY hg38.bed /opt/

COPY convert_header.py /opt/bin/
COPY lumpy2merge.py /opt/bin/
COPY vcf2bedpe.py /opt/bin/

ENV PATH=$PATH:/miniconda/bin