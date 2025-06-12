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
