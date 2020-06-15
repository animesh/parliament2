# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q
RUN apt-get update -q

# Install necessary tools
RUN apt-get install -y \
    git \
    inkscape \
    libbz2-dev \
    liblzma-dev \
    python3-pip \
    r-base

# install pip libraries required for svviz2
RUN pip3 install pybind11
RUN pip3 install seqlib
RUN pip3 install numpy
RUN pip3 install rpy2
RUN pip3 install pysam cython

# Install svviz2
RUN pip3 install -U git+git://github.com/nspies/svviz2.git
RUN pip3 install --upgrade --force-reinstall genomeview

RUN svviz2