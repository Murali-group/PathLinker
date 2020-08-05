FROM ubuntu:20.04

# version 2: clone pathlinker code base and use saved bash script to run pathlinker 

# INSTALL  TOOLS
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

# INSTALL GIT
RUN apt update \
    && apt install -y git

USER root
# install miniconda, Python 

RUN mkdir /home/PathLinkerHome
WORKDIR /home/PathLinkerHome
RUN wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean --all && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc


# RUN bash /home/seluser/Miniconda2-latest-Linux-x86_64.sh -b -p \
#     && eval "$(/home/PathLinkerHome/anaconda3/bin/conda shell.bash hook)" 

# install PathLinker
WORKDIR /home/PathLinkerHome
RUN git clone https://github.com/ajshedivy/Pathlinker-project.git

USER root

WORKDIR /home/PathLinkerHome/Pathlinker-project