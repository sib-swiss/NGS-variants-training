FROM linuxserver/code-server:4.16.1

# Install base utilities
RUN apt-get update && \
    apt-get install -y wget bc tree && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget -O ~/miniconda.sh --quiet \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc

# intall required software with conda
COPY environment.yml ./

RUN conda env create \
    -n ngs-tools \
    -f environment.yml

RUN conda init 

COPY settings.json /config/data/User/
COPY keybindings.json /config/data/User/
RUN mkdir -p /config/project

# copy a different snpeff config file in order to download to /data (in stead of conda directory)
# comment the line below if running this container locally.
COPY snpEff.config /opt/conda/envs/ngs-tools/share/snpeff-5.0-1/snpEff.config
