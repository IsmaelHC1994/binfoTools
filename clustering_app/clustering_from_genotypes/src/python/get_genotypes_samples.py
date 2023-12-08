from class_vcf import Vcf
import argparse

################################# argparser ####################################

parser = argparse.ArgumentParser()

parser.add_argument("--VCF", "-f", type=str, required=True)
parser.add_argument("--ID", "-id", type=str, required=True)
parser.add_argument("--OUT", "-o", type=str, required=True)
args = parser.parse_args()


################################ DEF ###########################################

def extract_genotypes_from_variant_ids(dirfile, diroutput, id_file):
    vfile = Vcf(dirfile)
    (vfile.get_variants_using_id(id_file)).check_all_samples_presence_absence(diroutput, True)

################################ MAIN #########################################

extract_genotypes_from_variant_ids(args.VCF, args.OUT, args.ID)