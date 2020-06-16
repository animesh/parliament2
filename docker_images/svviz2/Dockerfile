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
    r-base \
    wget

# install pip libraries required for svviz2
RUN pip3 install pybind11
RUN pip3 install seqlib
RUN pip3 install numpy
RUN pip3 install rpy2

# install htslib
RUN wget -q https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
RUN tar -xjf htslib-1.10.2.tar.bz2 && rm htslib-1.10.2.tar.bz2
# RUN apt-get install -y \
#     libbz2-dev \
#     liblzma-dev \
#     libncurses5-dev
WORKDIR htslib-1.10.2
RUN ./configure --prefix=/
RUN make
RUN make install
WORKDIR /

# RUN pip3 install pysam
# RUN pip3 install cython

# RUN apt-get install -y \
#     curl

# RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
# ENV PATH=/miniconda/bin:${PATH}
# RUN conda update -y conda

# RUN conda config --add channels bioconda
# RUN conda install -c bioconda cython

# Install svviz2
# RUN pip3 install -U git+git://github.com/nspies/genomeview.git
RUN pip3 install -U git+git://github.com/nspies/svviz2.git

RUN svviz2