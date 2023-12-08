FROM rocker/rstudio:4.1.2

LABEL maintainer="ihc.europa@gmail.com"
LABEL description="Base R image with RStudio and python for performing clustering from genotype files "

# Safer bash scripts with 'set -euxo pipefail'
SHELL ["/bin/bash", "-euxo", "pipefail", "-c"]

RUN apt-get update -qq -y \
    && apt-get install --no-install-recommends -qq -y \
        software-properties-common \
        build-essential r-base \ 
        cmake \
        tree \ 
        && add-apt-repository ppa:deadsnakes/ppa \ 
        && apt-get install -y python3.11 \
        python3-pip 

WORKDIR /app

# copy in app folder
COPY  clustering_from_genotypes clustering_from_genotypes

# run packages install
RUN pip install snakemake
RUN Rscript clustering_from_genotypes/src/R/install_packages.R

# change default python3
RUN sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && sudo update-alternatives  --set python3 /usr/bin/python3.11
