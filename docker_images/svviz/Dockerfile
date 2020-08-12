# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

RUN apt-get update -q --fix-missing
# Install necessary tools
RUN apt-get install -y \
    curl \
    inkscape \
    librsvg2-bin \
    rsync

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install -c bioconda svviz
RUN conda install -c conda-forge biopython

ENV PATH=$PATH:/miniconda/bin