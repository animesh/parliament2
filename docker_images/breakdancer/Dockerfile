# Set the base image to Ubuntu 14.04
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get install -y apt-utils -q
RUN apt-get upgrade -y -q

# Install necessary tools
RUN apt-get install -y \
    build-essential \
    cmake \
    curl \
    git \
    zlib1g-dev

RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN conda config --add channels bioconda
RUN conda install -c bioconda samtools

# Install breakdancer
RUN git clone --recursive https://github.com/genome/breakdancer.git
WORKDIR breakdancer
RUN mkdir build
WORKDIR build
RUN cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local
RUN make 
RUN make install

WORKDIR /
