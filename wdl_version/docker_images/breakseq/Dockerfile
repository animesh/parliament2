# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

# Install necessary tools
RUN apt-get install -y \
    curl

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN conda config --add channels bioconda
RUN conda install -c bioconda breakseq2 samtools

ENV PATH=$PATH:/miniconda/bin

# include breakpoint files required for BreakSeq
COPY breakseq2_bplib_20150129.hs37d5 /breakseq2_bplib_20150129.hs37d5
COPY breakseq2_bplib_20150129.hg19 /breakseq2_bplib_20150129.hg19
COPY bplib.hg38.gff /