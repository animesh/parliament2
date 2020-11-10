# Set the base image to Ubuntu 20.04
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update -q
RUN apt-get upgrade -y -q

RUN apt-get update -q --fix-missing
# Install necessary tools
RUN apt-get install -y \
    python3 \
    python3-pip

# Install truvari
RUN pip3 install Truvari