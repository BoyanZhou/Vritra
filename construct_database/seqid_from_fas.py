import os
import sys
# import logging.handlers
# from ete3 import NCBITaxa
# import re
# from Bio import SeqIO


"""
给定fas，提取出其中的seq id
# metagenomic requires, then diamond v2.0.15.153 will be loaded automatically
module add python/cpu/3.7.2
# !!! but for this script, we must use python 3.6.5, because only under this python, ete3 is available 
module add python/cpu/3.6.5
"""


def protein_id_from_fas(fas_path, fas_type):
    """
    :param fas_path:
    :param fas_type: fas_type must be 'InterPro' or 'UniProt'.
    :return: protein_id_list like, ["A0A371S335", "A0A089PYJ7"]
    """
    if fas_type not in ['InterPro', 'UniProt']:
        raise ValueError("Invalid option. fas_type must be 'InterPro' or 'UniProt'.")

    # 1. get the headers of all sequences, startwith(">")
    # # UniRef90 format: ">UniRef90_A0A371S335 Diaminopropionate ammonia-lyase n=2 Tax=Bacillaceae TaxID=186817 RepID=A0A371S335_9BACI"
    # # InterPro format: ">A0A089PYJ7|unreviewed|D-phenylhydantoinase|taxID:158822"
    protein_id_list = []    # like, ["A0A371S335", "A0A089PYJ7"]

    with open(fas_path, "r") as fas_f:
        for line in fas_f:
            if line.startswith(">"):
                if fas_type == "InterPro":
                    protein_id = line.split("|")[0][1:].strip()
                else:
                    protein_id = line.split(" ")[0]     # >UniRef90_A0A371S335
                    if protein_id.startswith(">UniRef"):
                        protein_id = protein_id.split("_")[-1]
                    else:
                        continue
                protein_id_list.append(protein_id)
    return protein_id_list


if __name__ == "__main__":
    mode = sys.argv[1]

    if mode == "id_from_fas":
        fas_path, fas_type = ["/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/uric_acid_gene_set/uric_acid_gene_set_to_uniref90_P66899_all_seq_after_iteration.fas", "UniProt"]
        Id_list_temp = protein_id_from_fas(fas_path, fas_type)
        # full_mapping_file_path = "/gpfs/data/lilab/home/zhoub03/software/UniRef/idmapping_selected.tab"
        # output_subset_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/uric_acid_gene_set/uric_acid_gene_set_to_uniref90_P66899_all_seq_after_iteration_id_mapping.tab"
        # subset_mapping_file(full_mapping_file_path, output_subset_path, Id_list_temp, "UniRef90")


