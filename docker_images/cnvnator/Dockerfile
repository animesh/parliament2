# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

# Install necessary tools to download root
RUN apt-get install -y \
    git

RUN rm -rf root/

# Install root
# Step 1: get source code
RUN git clone https://github.com/root-project/root.git
WORKDIR root/
RUN git checkout -b v6-20-04 v6-20-04

# Step 1A: download everything necessary to compile root
RUN apt-get install -y \
    # binutils \
    build-essential \
    cmake \
    # dpkg-dev \
    g++ \
    gcc \
    git \
    # libjpeg \
    # liblzma-dev \
    # libpng \
    libx11-dev \
    libxext-dev \
    libxft-dev \
    libxpm-dev \
    # libzstd-dev \
    python3-dev \
    wget

# Step 2: build commands
WORKDIR /
RUN mkdir build
WORKDIR build
RUN cmake /root
RUN make -j2
WORKDIR /

# Step 3: source root
RUN ["/bin/bash", "-c", "source /build/bin/thisroot.sh"]
ENV ROOTSYS=/build
ENV PATH=$PATH:$ROOTSYS/bin
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/build/lib

# Install CNVnator
# Step 0: download everything needed to download CNVnator
RUN apt-get install -y \
    unzip

# Step 1: download CNVnator
RUN wget -q https://github.com/abyzovlab/CNVnator/releases/download/v0.4.1/CNVnator_v0.4.1.zip
RUN unzip CNVnator_v0.4.1.zip

# Step 1A: download everything needed to install CNVnator
RUN wget -q https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar -xjf samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2
RUN apt-get install -y \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev
WORKDIR samtools-1.10
RUN ./configure --prefix=/
RUN make
RUN make install

# Step 2: install CNVnator
WORKDIR /CNVnator_v0.4.1/src/
RUN ln -s /samtools-1.10/ samtools
RUN apt-get install -y \
    libcurl4-openssl-dev \
    libreadline-dev 
RUN make