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
    autoconf \
    automake \
    build-essential \
    cmake \
    default-jdk \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    wget \
    zlib1g-dev

# Install Jasmine
RUN git clone https://github.com/mkirsche/Jasmine.git --branch v1.0.2
RUN Jasmine/build.sh