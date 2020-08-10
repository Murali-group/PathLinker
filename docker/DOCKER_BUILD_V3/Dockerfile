FROM ubuntu:20.04

# version 1: create and activate env on run 

# INSTALL TOOLS
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

# INSTALL GIT
RUN apt update \
    && apt install -y git

USER root


# install miniconda and configure bashrc for pathlinker env 
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.6.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    conda init bash

# install PathLinker
WORKDIR /home/PathLinker
COPY *.py /home/PathLinker/

# link example data 
COPY /example /home/PathLinker/example/

COPY docker/tests/minimal_env.yml tests/
RUN conda env create -f tests/minimal_env.yml && \
    echo "conda activate path-linker" >> ~/.bashrc 

USER root

WORKDIR /home/PathLinker