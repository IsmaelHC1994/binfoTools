#
# Clusters samples from a genotype file and generates a dendrogram and a barplot showing the accumulation of variants in the population
# 

############## file TARGETS #################################

# main file name
MAIN ?= example

# file with variants id
GENOTYPE ?= results/genotypes/${MAIN}

# number of clusters for NBCLUST (search for optimal number of clusters)
N_CLUSTERS ?= 3

# Output dir for plots
# PLOT_DIR ?= results/genotype_plots/$(basename ${GENOTYPE})/
PLOT_DIR ?= results/genotype_plots/${MAIN}/

############## make  ##############################################

# check the make version
ifeq ($(origin. RECIPEPREFIX), undefined)
  $(error "### Error! Please use GNU Make 4.0 or later ###")
endif

# Makefile customization
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information
usage::
> @echo "#"
> @echo "# dendrogram_from_genotypes.mk: Clusters samples from a genotype file and generates a dendrogram and a barplot showing the accumulation of variants in the population."
> @echo "#"
> @echo "# make -f dendrogram_from_genotypes.mk get GENOTYPE=${GENOTYPE} N_CLUSTERS=${N_CLUSTERS} PLOT_DIR=${PLOT_DIR} "
> @echo "#"

############## rules  ##############################################

# generate plots:
${PLOT_DIR}:
> mkdir -p ${PLOT_DIR}
> Rscript src/R/dendrogram_from_genotypes.R -f ${GENOTYPE} -n ${N_CLUSTERS} -od ${PLOT_DIR}

# list files
get:: ${PLOT_DIR}
> ls -lh $(dir ${PLOT_DIR})*

# open plots
open_plot:: 
> ls ${PLOT_DIR}*.tiff | xargs -n 1 xdg-open 

# remove files generated
get!::
> rm -r ${PLOT_DIR}

# targets that are not files.
.PHONY: usage install get get! open_plot