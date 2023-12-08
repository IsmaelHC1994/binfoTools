import argparse
from pathlib import Path
from class_variant import Variant
from class_vcf import Vcf
from os import path

################################# argparser ####################################

# parser = argparse.ArgumentParser()

# parser.add_argument("--VCF", "-f", type=str, required=True, help="vcf file")
# parser.add_argument("--sample-groups", "-sg", required=True, help="file with samples as rownames and the cluster they belong to in the first column")
# args = parser.parse_args()


################################ DEF ###########################################
# dir_file:str , dir_groups: str
def recount_variants_per_group():
    """_summary_

    Args:
        dir_file (str): _description_
        dir_groups (str): _description_
    """       
    dir_output = path.dirname('/home/ihenarejos/workspace/phd_summary/pipe_dendrogram/wes_of_sig_genes_val_exp.vcf') + "/sample_group_analysis.tsv"

    # variant info/gene file
    genes = {}
    with open('/home/ihenarejos/workspace/phd_summary/pipe_dendrogram/variants_sig_genes_wide_table.tsv', 'r') as info_f:
        lines = info_f.readlines()[1:]
        for line in lines:
            chr, pos, ref, alt, type, gene, *_ = line.strip().split("\t")
            variant_id = "_".join(['chr'+chr, pos, ref, alt])
            if gene not in genes:
                genes[gene] = []
            genes[gene].append((variant_id, type))
            

    # groups file
    groups = {str(i):set() for i in range(1, 5)}
    with open('/home/ihenarejos/workspace/phd_summary/pipe_dendrogram/clusters_val_exp.tsv', 'r') as clusters:
        lines = clusters.readlines()[1:]
        for line in lines:
            sample, group = line.strip().split("\t")
            groups[group].add(sample)
    
    # vcf file
    vcf_file = Vcf('/home/ihenarejos/workspace/phd_summary/pipe_dendrogram/wes_of_sig_genes_val_exp.vcf')
    
    
    recount = []
    for variant in vcf_file.list_variants:
        recount.append(Variant(variant, sample_list=vcf_file.sample_list).count_presence_per_group(groups))
    
    gene_groups = {gene:[0, 0, 0, 0] for gene in genes}
    for i in recount:
        variant, *group_res = i.split("\t")
        group_total = sum([int(i.strip()) for i in group_res])
        if group_total > 40:
            continue
        for gene in genes:
            gene_values = genes[gene]
            for v, _ in gene_values:
                if variant == v:
                    gene_groups[gene] = [a + int(b.strip()) for a, b in zip(gene_groups[gene], group_res)]
    
    print(gene_groups)
    # with open(dir_output, 'w') as out:    
    #     for variant in vcf_file.list_variants:
    #         recount = Variant(variant, sample_list=vcf_file.sample_list).count_presence_per_group(groups)
    #         out.write(recount)
                
################################ MAIN ###################################

recount_variants_per_group()