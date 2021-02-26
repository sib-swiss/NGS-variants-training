FROM jupyter/minimal-notebook

USER root
RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc
USER $NB_UID

RUN rm -r /home/jovyan/work

RUN /opt/conda/bin/conda install -y --quiet -c bioconda \
python=3.8 \
samtools \
bwa \
snpeff \
gatk4

COPY snpEff.config /opt/conda/share/snpeff-5.0-0
