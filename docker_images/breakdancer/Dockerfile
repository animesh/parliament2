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
    cmake \
    gcc
    git

# Install breakdancer
RUN git clone --recursive https://github.com/genome/breakdancer.git
WORKDIR breakdancer
RUN mkdir build
RUN cmake --help
WORKDIR build
RUN cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local
RUN make 
RUN make install

WORKDIR /
