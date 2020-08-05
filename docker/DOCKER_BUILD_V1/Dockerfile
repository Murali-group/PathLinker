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

# install PathLinker
WORKDIR /home/PathLinkerHome
RUN git clone https://github.com/ajshedivy/Pathlinker-project.git

# install miniconda and configure bashrc for pathlinker env 
WORKDIR /home/PathLinkerHome/Pathlinker-project/example
RUN wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda env create -f /home/PathLinkerHome/Pathlinker-project/example/minimal_env.yml" >> ~/.bashrc && \
    echo "conda activate path-linker" >> ~/.bashrc 


USER root

WORKDIR /home/PathLinkerHome/Pathlinker-project