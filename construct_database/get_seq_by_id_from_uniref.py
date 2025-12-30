"""
此脚本用来从UniRef文件中提取所需要蛋白质序列，序列编号的来源是blast的结果文件的第二列
module add diamond/0.9.18
module add python/cpu/3.6.5

get_seq_from_fas_given_id
get_seq_from_fas_avoid_id

get_prot_id_from_blast_res


"""

import os
import re
import sys


# Auxiliary Function_1: Output fas to path and return seq ID-length dict
def get_seq_from_fas_given_id(protein_fas_path, output_fas_path, prot_id_list):
    """
    Given the UniRef protein ID, extract protein fas from the large set of fas
    :param protein_fas_path = "/gpfs/data/lilab/home/zhoub03/software/UniRef/uniref90.fasta"
    :param output_fas_path: 
    :param prot_id_list=["UniRef100_Q197F8", "UniRef100_Q197F7", "UniRef100_Q6GZX2"]
    :return: 
    """
    # we tried the following method, but seems too time consuming
    # record_dict = SeqIO.index(protein_fas_path, "fasta")    #
    # record_dict.close()
    prot_id_set = set(prot_id_list)
    output_f = open(output_fas_path, "w")

    with open(protein_fas_path, "r") as protein_f:
        # >UniRef100_Q197F7 Uncharacterized protein 003L n=1 Tax=Invertebrate iridescent virus 3 TaxID=345201 RepID=003L_IIV3
        whether_output = False
        while True:
            line = protein_f.readline()
            if line == '':  # End of file (EOF)
                break
            if line.startswith(">"):
                # Header of sequence
                prot_id = line.split(" ")[0][1:]    # "UniRef100_Q197F7"
                if prot_id in prot_id_set:
                    whether_output = True
                    output_f.write(line)
                else:
                    whether_output = False
            else:
                # Protein sequence: MARPLLGKTSSVRRRLESLS
                if whether_output:
                    output_f.write(line)
    output_f.close()


# Auxiliary Function_1b: Output fas to path and return seq ID-length dict
def get_seq_from_fas_avoid_id(protein_fas_path, output_fas_path, avoid_id_set):
    """
    Given the UniRef protein ID set, extract protein fas from the large set of fas excluding these ids
    :param protein_fas_path = "/gpfs/data/lilab/home/zhoub03/software/UniRef/uniref90.fasta"
    :param output_fas_path:
    :param prot_id_list=["UniRef100_Q197F8", "UniRef100_Q197F7", "UniRef100_Q6GZX2"]
    :return:
    """
    # we tried the following method, but seems too time consuming
    # record_dict = SeqIO.index(protein_fas_path, "fasta")    #
    # record_dict.close()
    output_f = open(output_fas_path, "w")
    remaining_id_list = []

    with open(protein_fas_path, "r") as protein_f:
        # >UniRef100_Q197F7 Uncharacterized protein 003L n=1 Tax=Invertebrate iridescent virus 3 TaxID=345201 RepID=003L_IIV3
        whether_output = False
        while True:
            line = protein_f.readline()
            if line == '':  # End of file (EOF)
                break
            if line.startswith(">"):
                # Header of sequence
                prot_id = line.split(" ")[0][1:]    # "UniRef100_Q197F7"
                if prot_id in avoid_id_set:
                    whether_output = False
                else:
                    whether_output = True
                    remaining_id_list.append(prot_id)
                    output_f.write(line)
            else:
                # Protein sequence: MARPLLGKTSSVRRRLESLS
                if whether_output:
                    output_f.write(line)
    output_f.close()
    return remaining_id_list


def deduplicate_fasta(input_fasta_path, filtered_fasta_path):
    """
    Processes a FASTA file to keep only unique contig IDs.
    Args:
        input_fasta_path (str): Path to the input FASTA file.
        filtered_fasta_path (str): Path to the output filtered FASTA file.
    Returns:
    """
    contig_id_counts = {}
    input_f = open(input_fasta_path, "r")
    output_f = open(filtered_fasta_path, "w")
    whether_dup = False
    for line in input_f:
        if line.startswith(">"):
            uniref_id = line.split(" ")[0][1:]  # uniref id without ">"
            if uniref_id not in contig_id_counts:
                contig_id_counts.update({uniref_id: 1})
                whether_dup = False
            else:
                contig_id_counts[uniref_id] += 1
                whether_dup = True
        if not whether_dup:
            output_f.write(line)
    input_f.close()
    output_f.close()
    # report whether unique and duplicated id
    unique_id = []
    dup_id_dict = {}
    for uniref_id, count in contig_id_counts.items():
        if count == 1:
            unique_id.append(uniref_id)
        else:
            dup_id_dict.update({uniref_id: count})
    print(f"The unique uniref IDs are {unique_id}")
    print(f"The duplicated uniref IDs are {dup_id_dict}")
    return unique_id, dup_id_dict
# Example usage:
# deduplicate_fasta("frc_iteration_final_set.fas", "frc_iteration_final_set_dedup.fas")
# deduplicate_fasta("oxc_iteration_final_set.fas", "oxc_iteration_final_set_dedup.fas")
# deduplicate_fasta("oxdd_iteration_final_set.fas", "oxdd_iteration_final_set_dedup.fas")




# Auxiliary Function_2 no filter (Single Seq Blast):
# Get the matched protein UniRef90 ID from the blast results
def get_prot_id_from_blast_res(blast_res_path):
    prot_id_list = []
    with open(blast_res_path, "r") as blast_res_f:
        for line in blast_res_f:
            # "A0A014MFG0|unreviewed|Amidohydrolase-related    UniRef90_UPI0015C61A5F  70.068  441     132     0       1       441     1       441     0.0     646"
            cols = line.split("\t")
            prot_id_list.append(cols[1])
    return prot_id_list


# Auxiliary Function_3 (Single Seq Blast, for initial query protein ID):
# Get the matched protein UniRef90 ID from the blast results after certain filter
def get_prot_id_from_blast_res_filtered(blast_res_path, query_len_dict, matched_len_dict, identity_threshold=50, coverage_threshold=0.8):
    prot_id_list = []
    with open(blast_res_path, "r") as blast_res_f:
        for line in blast_res_f:
            # "A0A014MFG0|unreviewed|Amidohydrolase-related    UniRef90_UPI0015C61A5F  70.068  441     132     0       1       441     1       441     0.0     646"
            cols = line.split("\t")
            identity = float(cols[2])
            query_coverage = float(cols[3]) / query_len_dict[cols[0]]
            matched_coverage = float(cols[3]) / matched_len_dict[cols[1]]
            if identity >= identity_threshold and query_coverage >= coverage_threshold and matched_coverage >= coverage_threshold:
                prot_id_list.append(cols[1])
    return prot_id_list


# Auxiliary Function_4 (Multi Seq Blast):
# Get the sequence set of both query and matched sequences from blast result file
def get_blasted_and_matched_prot_id(blast_res_path, prot_len_dict, identity_threshold=50, coverage_threshold=0.8):
    prot_id_query = set()
    prot_id_macthed = set()
    with open(blast_res_path, "r") as blast_res_f:
        for line in blast_res_f:
            # "A0A014MFG0|unreviewed|Amidohydrolase-related    UniRef90_UPI0015C61A5F  70.068  441     132     0       1       441     1       441     0.0     646"
            cols = line.split("\t")
            identity = float(cols[2])
            query_coverage = float(cols[3])/prot_len_dict[cols[0]]
            matched_coverage = float(cols[3])/prot_len_dict[cols[1]]
            if identity >= identity_threshold and query_coverage >= coverage_threshold and matched_coverage >= coverage_threshold:
                prot_id_query.add(cols[0])
                prot_id_macthed.add(cols[1])
    return prot_id_query, prot_id_macthed


# Auxiliary Function_5 (Single or Multi Seq):
# Get the sequence ID-length dict
def get_seq_length_from_fas(fas_path, seq_type):
    """
    :param fas_path
    :param seq_type: must be "nucl" or "prot"
    :return: dict of {protein ID: length}
    """
    if seq_type not in ["nucl", "prot"]:
        raise TypeError(f"Error! The sequence type in {fas_path} must be 'nucl' or 'prot'!")
    contig_lengths = {}
    with open(fas_path, "r") as fas_f:
        current_id = None
        current_length = 0
        for line in fas_f:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if current_id is not None:
                    contig_lengths[current_id] = current_length
                current_id = line.split("\t")[0].split(" ")[0][1:]   # Get the sequence ID without ">"
                current_length = 0      # Reset the sequence length
            else:
                current_length += len(line)     # Add the length of the sequence line
        if current_id is not None:              # Add the last contig
            contig_lengths[current_id] = current_length
    # check whether it is nucleotide sequence
    if seq_type == "nucl":
        contig_lengths = {i: int(j/3) for i, j in contig_lengths.items()}
    return contig_lengths


# Auxiliary Function_6 (Single or Multi Seq):
# Get the sequence UniRef ID
def get_seq_id_from_fas(fas_path):
    """
    :param fas_path
    :return: dict of {protein ID: length}
    """
    uniref_id_list = []
    with open(fas_path, "r") as fas_f:
        for line in fas_f:
            if line.startswith(">"):  # Header line
                uniref_id = line.strip().split(" ")[0][1:]  # remove ">"
                uniref_id_list.append(uniref_id)
    return uniref_id_list



"""
def get_blasted_prot_id_from_files(blasted_file_path_list):
    total_blasted_prot = set()
    for blasted_file_path in blasted_file_path_list:
        blasted_prot_set = set()
        with open(blasted_file_path, "r") as blasted_f:
            for line in blasted_f:
                blasted_prot_set.add(line.split("\t")[0])
        # merge the prot from single file to the total set
        total_blasted_prot = total_blasted_prot | blasted_prot_set
    return total_blasted_prot

    seq_record_dict = SeqIO.to_dict(SeqIO.parse(protein_fas_path, "fasta"))     # {seq_id: seq_record}
    # SeqRecord(seq=Seq('MGDIMRPVPFKQLLCWIAEEYRSQWTIFGIPESQFFIKENGKSIQIFDESCATP...RPV', SingleLetterAlphabet()), id='X1ALD3|unreviewed|Selenate', name='X1ALD3|unreviewed|Selenate',
    # description='X1ALD3|unreviewed|Selenate reductase subunit YgfK|taxID:412755', dbxrefs=[])
    output_f = open(output_fas_path, "w")
    pattern = r'taxID:(\d+)'
    for seq_id, seq_record in seq_record_dict.items():
        match = re.search(pattern, seq_record.description)
        # Check if a match is found
        if match:
            # the seq matches the given taxID
            if match.group(1) == taxID:
                SeqIO.write(seq_record, output_f, "fasta")
        else:
            print(f"{seq_record.description} has no taxon ID")
"""


if __name__ == "__main__":
    pass