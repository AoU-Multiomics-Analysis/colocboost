FROM mambaorg/micromamba:1.5.3
## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

micromamba -y -n base install  \
    conda-forge::r-tidyverse \
    dnachun::r-colocboost \
    conda-forge::r-data.table \
    bioconda::r-bedr \
    bioconda::bedtools


RUN wget https://github.com/evin-padhi/colocboost_WDL/blob/main/utils/colocboost_utils.R
RUN wget https://github.com/evin-padhi/colocboost_WDL/blob/main/utils/run_colocboost.R
