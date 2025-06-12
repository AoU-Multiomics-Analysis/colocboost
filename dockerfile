FROM mambaorg/micromamba:1.5.3
## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
USER root


RUN micromamba -y -n base install  \
    conda-forge::r-base \ 
    conda-forge::r-tidyverse \
    dnachun::r-colocboost \
    conda-forge::r-data.table \
    conda-forge::r-bedr \
    bioconda::bedtools \
    conda-forge::r-argparse

#RUN apt-get install -y wget

RUN curl  https://raw.githubusercontent.com/evin-padhi/colocboost_WDL/refs/heads/main/utils/colocboost_utils.R -o ./colocboost_utils.R
RUN curl  https://raw.githubusercontent.com/evin-padhi/colocboost_WDL/refs/heads/main/utils/run_colocboost.R -o ./run_colocboost.R
