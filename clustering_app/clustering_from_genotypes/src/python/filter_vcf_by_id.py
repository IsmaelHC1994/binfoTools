from class_vcf import Vcf
from class_variant import Variant
import argparse

################################# argparser ####################################

parser = argparse.ArgumentParser()

parser.add_argument("--VCF", "-f", type=str, required=True)
parser.add_argument("--ID", "-id", type=str, required=True)
parser.add_argument("--OUT", "-o", type=str, required=True)
args = parser.parse_args()


################################ DEF ###########################################

def find_me_variants(dirfile, diroutput, id_file):
    vfile = Vcf(dirfile)
    (vfile.get_variants_using_id(id_file)).save_vcf(diroutput)

################################ MAIN #########################################

find_me_variants(args.VCF, args.OUT, args.ID)
