FROM mambaorg/micromamba:1.5.3
## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
USER root


RUN micromamba -y -n base install  \
    conda-forge::r-base \ 
    conda-forge::r-tidyverse \
    dnachun::r-colocboost \
    dnachun::r-pecotmr \
    conda-forge::r-data.table \
    conda-forge::r-bedr \
    bioconda::bedtools \
    conda-forge::r-argparse \
    conda-forge::r-janitor \
    bioconda::htslib \ 
    conda-forge::r-zoo \
    dnachun::r-twosamplemr 


#RUN apt-get install -y wget

COPY utils/* . 

RUN mkdir -p /src

WORKDIR /src


