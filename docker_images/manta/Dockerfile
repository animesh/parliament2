# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q
RUN apt-get update -q

# Install necessary tools
RUN apt-get install -y -q \
    bzip2 \
    gcc \
    g++ \
    make \
    python \
    wget \
    zlib1g-dev

# Install Manta
ENV MANTA_VERSION=1.6.0
RUN wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
RUN tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2
RUN mkdir build
WORKDIR build
RUN ../manta-${MANTA_VERSION}.release_src/configure --jobs=4
RUN make -j4 install