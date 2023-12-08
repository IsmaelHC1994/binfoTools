import argparse
from pathlib import Path
from os import path

################################# argparser ####################################

parser = argparse.ArgumentParser()

parser.add_argument("--TSV", "-f", type=str, required=True, help="table with variants in wide format (tab separated with the first 4 cols: chromosome, position, reference allele, alternative allele)")
parser.add_argument("--output", "-o", required=True, type=str, help="output file")
parser.add_argument("--chrFormat", "-c", type=bool, required=False, action=argparse.BooleanOptionalAction, help=".adds 'chr' before chromosome number. by default, leaves chr field as only the number")
args = parser.parse_args()


################################ DEF ###########################################

def generate_ids(dirfile,diroutput, add_chr='True'):
    """ generate a variant ids file (chr_pos_ref_alt) from a tsv file that have the first 4 columns in the next order: chr, pos, reference allele, alternative allele

    Args:
        dirfile (str): route to tsv file that contains the variants in wide format
        add_chr (bool): if True, "chr" is added to chromosome number
        diroutput (str): route to save output file containing only variants id
    """    
    out_lines = []
    
    # diroutput = path.dirname(dirfile) + Path(dirfile).stem + "_variant_ids.txt"

    if add_chr:
        chr_format = 'chr'
    else:
        chr_format = '' 
    
    with open(dirfile, 'r') as file_in:
        
        for line in file_in.readlines()[1:]:
            chr, pos, ref, alt, *_ = line.strip().split("\t")
            out_lines.append(("_").join([chr_format + chr, pos, ref, alt]) + "\n")
            
    with open(diroutput, 'w') as file_out:
        
        file_out.writelines(out_lines)

################################ MAIN ###################################

generate_ids(args.TSV, args.output, args.chrFormat)
