#
# Workflow to generate a dendrogram and a barplot from clustering samples according to variant genotype.
# 

############## file TARGETS #################################

# filename that will be used for the entire workflow
MAIN_NAME ?= example

# Wide table 
WIDE_TABLE ?= data/wide_tables/${MAIN_NAME}.tsv

# VCF file
VCF ?= data/vcf/${MAIN_NAME}.vcf

# ID file
ID ?= data/variants_id/${MAIN_NAME}.txt

# genotype files
GENOTYPES ?= results/genotypes/${MAIN_NAME}.tsv

# minimum number of clusters for NBCLUST (search for optimal number of clusters)
N_CLUSTERS ?= 3

# Output directory for plots
PLOT_DIR ?= results/genotype_plots/${MAIN_NAME}/

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
> @echo "# workflow_dendrogram_genotypes.mk: Workflow to extract genomic profiles from a VCF file using variants ids (chr_) clusters samples from a genotype file and generates a dendrogram and a barplot showing the accumulation of variants in the population."
> @echo "#"
> @echo "# make all: runs the entire worflow. Adding the -n command will display all steps of the pipeline."
> @echo "#"
> @echo "##### Specific rules usage #####"
> @echo "#"
> @echo "### variants_id ###"
> @echo "# Generates a file with variant ids in the format chr_pos_refAllele_altAllele from a wide table containing the corresponding tab separated fields"
> @echo "# make -f workflow_dendrogram_genotypes.mk variants_id WIDE_TABLE=${WIDE_TABLE} ID=${ID} "
> @echo "#"
> @echo "### genotypes ###"
> @echo "# Retrieves the genomic profile (heterozygous for the variant/homozygous for the variant/ref allele) from a VCF file using a file with variant ids"
> @echo "# make -f workflow_dendrogram_genotypes.mk genotypes VCF=${VCF} ID=${ID} GENOTYPES=${GENOTYPES}"
> @echo "#"
> @echo "### plots ###"
> @echo "# Launches a Rscript to perform hierarchical clustering using a selected number of optimal clusters to look for"
> @echo "# make -f workflow_dendrogram_genotypes.mk plots PLOT_DIR=${PLOT_DIR} N_CLUSTERS=${N_CLUSTERS} GENOTYPES=${GENOTYPES}"
> @echo "#"
> @echo "#"


############## rules  ##############################################


#### Generating variants id from a wide table task:

variants_id:

# make a directory to store variant_ids files:
> mkdir -p data/variants_id

# generate the ids
> python3 src/python/variant_ids_from_wide_format.py -f ${WIDE_TABLE} -o ${ID} --chrFormat 


#### extracting genotypes task:

genotypes:

# make a directory and then generate genotypes file
> mkdir -p results/genotypes
> python3 src/python/get_genotypes_samples.py -f ${VCF} -id ${ID} -o ${GENOTYPES}


#### generate plots task:

plots:

# make a directory and store plots
> mkdir -p ${PLOT_DIR}
> Rscript src/R/dendrogram_from_genotypes.R -f ${GENOTYPES} -n ${N_CLUSTERS} -od ${PLOT_DIR}


#### open plots
# wont work in docker as it requires a display  

open_plots:
> ls ${PLOT_DIR}*.tiff | xargs -n 1 xdg-open 


#### Run all commands in one shot:

all: variants_id genotypes plots

#### remove files generated

clean:
> rm -rf ${ID} ${GENOTYPES} ${PLOT_DIR}


# targets that are not files.
.PHONY: usage install clean variants_id genotypes plots open_plots
