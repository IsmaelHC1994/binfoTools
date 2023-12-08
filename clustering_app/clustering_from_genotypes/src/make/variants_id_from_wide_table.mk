#
# Generates a file with a variant id (chr_pos_ref_alt) per line, obtained from a wide tablet supplied.
# 

############## file TARGETS #################################

# Wide table 
WIDE_TABLE ?= data/example_wide_table.tsv

# ID file to output
ID ?= data/id_example_.txt

############## make  ##############################################


# Check the make version.
ifeq ($(origin .RECIPEPREFIX), undefined)
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
> @echo "# variants_id_from_wide_table.mk: Generates a file with a variant id (chr_pos_ref_alt) per line, obtained from a wide tablet supplied."
> @echo "#"
> @echo "# make -f variants_id_from_wide_table.mk get WIDE_TABLE=${WIDE_TABLE} ID=${ID} "
> @echo "#"

############## rules  ##############################################

# generate IDS
${ID}:
> python3 src/python/variant_ids_from_wide_format.py -f ${WIDE_TABLE} -o ${ID} --chrFormat 

# list the data.
get:: ${ID}
> @ls -lh $(dir ${ID})*

# remove files generated
get!::
> rm -f ${ID}

# targets that are not files.
.PHONY: usage install get get!