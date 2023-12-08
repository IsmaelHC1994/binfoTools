#
# Extract samples genotypes from a VCF using a list of variants id
# 

############## file TARGETS #################################

# VCF file
VCF ?= data/vcf/example.vcf 

# ID file with variants id (chr_pos_ref_alt)
ID ?= data/variants_id/id_example.txt

# genotype files
GENOTYPES ?= results/genotypes/example_genotypes.tsv

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
> @echo "# genotypes_from_variants.mk: Uses a VCF with a variant id file to extract sample genotype."
> @echo "#"
> @echo "# make -f genotypes_from_variants.mk get WIDE_TABLE=${VCF} ID=${ID} GENOTYPES=${GENOTYPES} "
> @echo "#"

############## rules  ##############################################

# obtain genotypes
${GENOTYPES}:
> mkdir -p results/genotypes
> python3 src/python/get_genotypes_samples.py -f ${VCF} -id ${ID} -o ${GENOTYPES}

# list files:
get:: ${GENOTYPES}
> ls -lh $(dir ${GENOTYPES})*

# remove files generated
get!::
> rm -f ${GENOTYPES}

# targets that are not files.
.PHONY: usage install get get!