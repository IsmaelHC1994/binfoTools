#
# Filter VCF to generate a smaller VCF with extracted variants according to a list of variants ID (chr_pos_ref_alt) 
# 

############## file TARGETS #################################

# The directory that stores the output VCF
OUT_DIR ?= results/filtered_vcf

# id name
NAME ?= example

# VCF file
VCF ?= data/vcf/${NAME}.vcf

# file with variants id
ID ?= data/id_example.txt

# Output name
OUT ?= ${OUT_DIR}/${NAME}_filtered.vcf


############## make  ##############################################


# Check the make version.
ifeq ($(origin .RECIPEPREFIX), undefined)
  $(error "### Error! Please use GNU Make 4.0 or later ###")
endif

# Makefile customizations.
.RECIPEPREFIX = >
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Print usage information
usage::
> @echo "#"
> @echo "# extract_variants_by_id.mk: Filter VCF to generate a smaller VCF with extracted variants according to a list of variants ID (chr_pos_ref_alt) "
> @echo "#"
> @echo "# make get VCF=${VCF} ID=${ID}"
> @echo "#"

############## rules  ##############################################


# Filter VCF
${OUT}:
> mkdir -p ${OUT_DIR}
> python3 src/python/filter_vcf_by_id.py -f ${VCF} -id ${ID} -o ${OUT}

# List the data.
get:: ${OUT}
> @ls -lh $(dir ${OUT})${NAME}*

# Remove files generated
get!::
> rm -f ${OUT} -r ${OUT_DIR}

# Targets that are not files.
.PHONY: usage install get get!