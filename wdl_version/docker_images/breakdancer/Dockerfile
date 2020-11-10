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

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN conda config --add channels bioconda
RUN conda install -c bioconda samtools

COPY bam2cfg /usr/bin
COPY breakdancer-max /usr/bin
RUN chmod +x /usr/bin/bam2cfg
RUN chmod +x /usr/bin/breakdancer-max

COPY ctx_to_vcf.py /opt/bin/
COPY merge_files.py /opt/bin/

ENV PATH=$PATH:/miniconda/bin