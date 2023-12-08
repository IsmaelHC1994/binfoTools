# load pickle module
import pickle, argparse, gzip
from collections import defaultdict

def parse_gtf_file(gtf_file: str) -> dict:
    """ Parses a gtf file to obtain a dictionary in the form of (start, end): gene id.

    Args:
        gtf_file (str): full route to gtf file

    Returns:
        dict: dictionary in the form of (start, end): gene id
    """
    
    gtf_dict = defaultdict(dict)
    
    with open(gtf_file, 'rt') as gtf_f:
        
            for line in gtf_f :
                chr, start, end, ensembl_id = line.strip().split('\t')
                
                gtf_dict[chr][(int(start), int(end))] = ensembl_id.replace('"', '').replace(';', '')
            
    return gtf_dict


if __name__ == '__main__':
    
    pars = argparse.ArgumentParser(description= 'Generate GTF pickle dictionary')
    pars.add_argument("-f", 
                    help="gtf file",
                    required=True)
    pars.add_argument("-o", 
                      help='output dictionary and filename', required=True)
    args = pars.parse_args()
    
    gtf_dict = parse_gtf_file(args.f)
    out = open(args.o, 'wb')
    pickle.dump(gtf_dict, out)
    out.close()