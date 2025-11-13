FROM ghcr.io/prefix-dev/pixi:latest

RUN pixi global install -c dnachun -c conda-forge -c bioconda r-base=4.4 bedtools htslib wget findutils


RUN pixi global install --environment r-base \
    r-argparse \
    r-argparser \
    r-bedr \
    r-colocboost \
    r-data.table \
    r-enrichr \
    r-gprofiler2 \
    r-janitor \
    r-optparse \
    r-pecotmr \
    r-tidyverse \
    r-susier

RUN find ${HOME}/.pixi/envs/r-base/bin -name '*bioconductor-*-post-link.sh' | \
   xargs -I % bash -c "PREFIX=${HOME}/.pixi/envs/r-base PATH=${HOME}/.pixi/envs/r-base/bin:${PATH} % || true"

COPY utils/* . 
