#
# Workflow to annotate a dendrogram and a barplot from clustering samples according to variant genotype.
# 

############## file TARGETS #################################

# filename project that will be used for the entire workflow
MAIN_NAME ?= example

# VCF dir
VCF_dir ?= data/${MAIN_NAME}

# GTF file
GTF_FILE ?= data/${MAIN_NAME}/${MAIN_NAME}.gtf.gz

# GTF dictionary
GTF_DICT ?= results/gtf_dictionaries/${MAIN_NAME}.pk

# output directory
OUTPUT_ANN ?= results/annotated_vcf/${MAIN_NAME}

# Number of processes
NUM_PROCESS ?= 1

# Pair of bases to search for genes up/downstream from mutation position
NUM_PB ?= 200000

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
> @echo "# workflow.mk: Workflow to annotate VCF files in a directory from the information stored in a GTF file."
> @echo "#"
> @echo "# make dictionary_from_gtf: Extracts gene information from a GTF file and saves it as a pickle dictionary."
> @echo "#"
> @echo "# make vcf_annotation: Annotates VCF so INFO field shows genes located near the mutation."
> @echo "#"
> @echo "# make all: runs the entire worflow. Adding the -n command will display all steps of the pipeline."

############## rules  ##############################################

#### Simplify gtf file:

simplify_gtf:

> zcat ${GTF_FILE} | awk -v OFS='\t' '{if($$3=="gene"){print $$1, $$4, $$5, $$10}}' > ${GTF_FILE}.tsv

#### Generate a pickle dictionary from a GTF file:

dictionary_from_gtf:

# make a directory to store the dictionary
> mkdir -p results/gtf_dictionaries

# generate the dictionary
> python3 src/python/pickle_gtf.py -f ${GTF_FILE}.tsv -o ${GTF_DICT}

#### Annotate VCF files with GTF information

vcf_annotation:

# make dictionary to store outputs
> mkdir -p ${OUTPUT_ANN}

> python3 src/python/vcf_gtf_annotation.py -VCF ${VCF_dir} -gtf ${GTF_DICT} -o ${OUTPUT_ANN} --num_process ${NUM_PROCESS} --pb ${NUM_PB}

#### Run all commands in one shot:

all: dictionary_from_gtf vcf_annotation

#### remove files generated

clean:
> rm -rf ${GTF_DICT} ${OUTPUT_ANN}


# targets that are not files.
.PHONY: usage install clean dictionary_from_gtf vcf_annotation 