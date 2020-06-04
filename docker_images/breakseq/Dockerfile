# Set the base image to Ubuntu 14.04
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get install -y apt-utils -q
RUN apt-get upgrade -y -q

# Install necessary tools
RUN apt-get install -y \
    python-pip

# Install BreakSeq2
RUN pip install https://github.com/bioinform/breakseq2/archive/2.2.tar.gz