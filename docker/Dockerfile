FROM ubuntu:20.04

# create conda environment in advance and activate on run

# install tools
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

# install miniconda and configure bashrc
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    conda init bash

WORKDIR /home/PathLinker

# create path-linker conda environment
COPY docker/minimal_env.yml ./docker/
RUN conda env create -f docker/minimal_env.yml && \
    echo "conda activate path-linker" >> ~/.bashrc 

# install PathLinker
COPY *.py ./
COPY *.md ./
