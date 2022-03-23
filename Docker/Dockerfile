FROM jupyter/minimal-notebook:latest

# change user to root to modify .bashrc for e.g. highlighting in the console
USER root
RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc

# change back to notebook user
USER $NB_UID

# remove auto-generated work directory
RUN rm -r /home/jovyan/work

# intall required software with conda
RUN /opt/conda/bin/conda install -y --quiet -c bioconda \
python=3.8 \
samtools=1.12 \
bwa=0.7.17 \
snpeff=5.0 \
gatk4=4.2.0.0

# copy a different snpeff config file in order to download to /data (in stead of conda directory)
# comment the line below if running this container locally.
COPY snpEff.config /opt/conda/share/snpeff-5.0-0
