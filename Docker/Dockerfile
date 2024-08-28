FROM linuxserver/code-server:4.92.2

# Install base utilities
RUN apt-get update && \
    apt-get install -y wget bc tree && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install mamba
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3
RUN /opt/miniforge3/condabin/mamba init

# Put conda in path so we can use conda activate
ENV PATH=/opt/miniforge3/condabin/:$PATH


# intall required software with conda
COPY environment.yml ./

RUN mamba env create \
    -n ngs-tools \
    -f environment.yml

RUN mamba init 

COPY settings.json /config/data/User/
COPY keybindings.json /config/data/User/
RUN mkdir -p /config/project

# copy a different snpeff config file in order to download to /data (in stead of conda directory)
# comment the line below if running this container locally.
COPY snpEff.config /opt/conda/envs/ngs-tools/share/snpeff-5.0-1/snpEff.config
