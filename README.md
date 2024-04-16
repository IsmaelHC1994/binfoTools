A repository to back up and describe workflows and scripts from several bioinformatic projects. 
In process of updating.

## Clustering genomic profiles

Generates a dendrogram and a barplot for sample clustering according to sample genotypes.

### Usage
0. Clone/Download repository
1. Build and run docker image from repository root:

       sudo docker build -t <user>/clustering-from-genotypes clustering_app/.
   
       sudo docker run -it -d --name dendrograms <user>/clustering-from-genotypes

       sudo docker exec -it dendrograms /bin/bash

2. In docker terminal, launch workflow either using MAKE or Snakemake (docker is configured to start in the root of the app):

       cd clustering_from_genotypes
       make -f src/make/workflow_dendrogram_genotypes.mk all  MAIN_NAME=foo_fop
       # When finished, you can run bash docker_stop.sh to stop and remove docker image.

* Alternatively, launch either MAKE or Snakemake from the app root (GNU Make > 4.0, python > 3.8, snakemake, R > 4.0 with provided packages in src/R/install_packages.R are required):

       make -f [path_to_app]/src/make/workflow_dendrogram_genotypes.mk all MAIN_NAME=foo_fop

       snakemake -c1 -s [path_to_app] src/snakemake/Snakefile all --config MAIN_NAME=foo_fop

## VCF-GTF annotation
Annotates VCF files (files containing mutations found in the genomic data of a specific organism) to include the genes affected by these mutations. Uses parallelization to perform analysis on several VCF files at the same time (with make -f NUM_PROCESS=<>).

1. Launching workflow:

       make -f src/make/workflow.mk all

## Make help

For any make file in this repository, typing  `make -f src/make/<makefile>.mk` in the terminal will print the usage help in the terminal. All rules of the workflow will be described with examples. Alternatively, the usage rule can be used to display the help:  `make -f src/make/workflow_dendrogram_genotypes.mk usage`

- To display what commands make is going to run, the variable -n should be used. For example, for the rule `genotypes`, typing `make -f src/make/workflow_dendrogram_genotypes.mk genotypes -n` will display the following:
  
         mkdir -p results/genotypes python3 src/python/get_genotypes_samples.py -f data/vcf/example.vcf -id data/variants_id/example.txt -o results/genotypes/example.tsv

