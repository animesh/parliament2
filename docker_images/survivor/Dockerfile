# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

# Install necessary tools
RUN apt-get update --fix-missing -q
RUN apt-get install -y \
    build-essential \
    git \
    zlib1g-dev

# Install SURVIVOR
RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git
WORKDIR /SURVIVOR/Debug/
RUN git checkout 1d1d33b016dbf818b1678a27dee3d3de7f0fda0b
RUN make
WORKDIR /

# Install vcftools
RUN apt-get update -y -q
RUN apt-get install -y \
    autoconf \
    g++ \
    pkg-config

RUN git clone https://github.com/vcftools/vcftools.git
WORKDIR vcftools
RUN ./autogen.sh
RUN ./configure --disable-dependency-tracking && make && make install
WORKDIR /

RUN apt-get install -y \
    python3

COPY combine_combined.py /opt/bin/
COPY all.phred.txt /opt/bin/
COPY correct_max_position.py /opt/bin/

ENV PATH=$PATH:/SURVIVOR/Debug/