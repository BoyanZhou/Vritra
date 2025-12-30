"""
原先想使用 pairwise alignment，但是速度太慢，处理数千的 pairwise 都可能需要几十个小时；
所以sequence较多的group使用blast，较少的则直接使用Biopython

module add diamond/0.9.18
module add python/cpu/3.6.5

"""


import pandas as pd
from Bio import SeqIO
import os
from Bio import pairwise2
import time


#################################################################
# combined function of Pairwise Identity By Biopython and Blast #
#################################################################
def pairwise_identity_by_two_method(all_seq_dict, output_fas_path, target_uniref_id_list, n_seq_threshold=20, if_symmetric=True):
    if len(target_uniref_id_list) <= n_seq_threshold:
        # small size using Biopython
        seq_dict = {i: all_seq_dict[i] for i in target_uniref_id_list}
        identity_res_df = seqs_pairwise_identity(seq_dict, symmetric=if_symmetric)
    else:
        # larger size using blast
        # # Subtract all centroid sequences within species
        subset_fas_given_uniref100(all_seq_dict, output_fas_path, target_uniref_id_list)
        # # Blast to it self to get pairwise identity file
        species_blast_res_path = blast_to_self(output_fas_path, num_alignments=len(target_uniref_id_list))
        identity_res_df = pd.read_table(species_blast_res_path, header=None)
        identity_res_df = identity_res_df.iloc[:, 0:3]
        identity_res_df.columns = ['Seq1', 'Seq2', 'Identity']
        # Attention, A-B and B-A may be different in blast result, remove A-A and B-B
        identity_res_df = identity_res_df.loc[identity_res_df['Seq1'] != identity_res_df['Seq2']]
    return identity_res_df


##################################
# Pairwise Identity By Biopython #
##################################
def seqs_pairwise_identity(seq_dict, symmetric=True):
    """
    Given the sequences within group, calculate pairwise identity
    seq_dict = seq_record_within_uniref90
    :param seq_dict: 'UniRef100_A0A0F0H8L9': SeqRecord(seq=Seq('MAKALEGVKVLDMTHVQSGPSSTQLLAWLGADVIKLETPGRGDITRGQLRDLPD...GVI'),
    id='UniRef100_A0A0F0H8L9', name='UniRef100_A0A0F0H8L9', description='UniRef100_A0A0F0H8L9 Formyl-CoA:oxalate
    CoA-transferase n=26 Tax=Lentzea aerocolonigenes (Lechevalieria aerocolonigenes) (Saccharothrixaerocolonigenes)
    TaxID=68170 RepID=A0A0F0H8L9', dbxrefs=[])
    :param symmetric: record both A-B and B-A
    :return:
    """
    time_start = time.time()
    pairwise_identity_list = []     # store the result in the format, ["UniRef100_A", ]
    uniref100_list = list(seq_dict)
    for i in range(len(uniref100_list) - 1):
        uniref100_i = uniref100_list[i]
        for j in range(i+1, len(uniref100_list)):
            uniref100_j = uniref100_list[j]
            iden_ij = identity_cal(seq_dict[uniref100_i].seq, seq_dict[uniref100_j].seq)
            pairwise_identity_list.append([uniref100_i, uniref100_j, iden_ij])
            if symmetric:
                pairwise_identity_list.append([uniref100_j, uniref100_i, iden_ij])
    time_cost = time.time() - time_start
    print(f"It took {time_cost} seconds to calculate the pairwise identity between {len(uniref100_list)} sequences.")
    pairwise_identity_df = pd.DataFrame(pairwise_identity_list, columns=['Seq1', 'Seq2', 'Identity'])
    return pairwise_identity_df
# identity_res_df = seqs_pairwise_identity(seq_record_within_uniref90)


# identity calculation by Biopython
def identity_cal(seq1, seq2):
    # Define two protein sequences
    # seq1 = Seq("MKTLLILAVVALLTQVQCE")
    # seq2 = Seq("MKALLILAVVALSTQVQC-")
    # Perform a global alignment
    alignments = pairwise2.align.globalxx(seq1, seq2)
    # Select the first alignment (highest score)
    aligned_seq1, aligned_seq2, score, start, end = alignments[0]
    # Calculate pairwise identity
    identical_positions = sum(i == j for i, j in zip(aligned_seq1, aligned_seq2))
    sequence_identity = (identical_positions * 2 / (len(aligned_seq1) + len(aligned_seq2))) * 100
    return sequence_identity


##############################
# Pairwise Identity By Blast #
##############################
# Step 1: Extract sequence given UniRef100 ID, using Biopython, !!! for small size fasta !!! (not used in new scripts)
def subset_fas_given_uniref100(all_seq_dict, output_fas_path, target_uniref_id_list):
    """
    :param all_seq_dict: the fas has been read into a dict by Biopython
    :param output_fas_path:
    :param target_uniref_id_list:
    :return:
    """
    # input_fas_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/FRC/frc_blast_iteration_final_set/frc_target_UniRef100.fas"
    # all_seq_dict = SeqIO.to_dict(SeqIO.parse(input_fas_path, "fasta"))  # {'UniRef100_G7ZBD1': seq_record}
    result_file = open(output_fas_path, "w")
    found_uniref_id_list = []
    not_found_uniref_id_list = []
    for uniref100_id in target_uniref_id_list:
        if uniref100_id in all_seq_dict:
            found_uniref_id_list.append(uniref100_id)
            SeqIO.write(all_seq_dict[uniref100_id], result_file, "fasta")
        else:
            not_found_uniref_id_list.append(uniref100_id)
    result_file.close()
    print(f"In identity_within_group, {len(found_uniref_id_list)} UniRef IDs are found. "
          f"{len(not_found_uniref_id_list)} UniRef IDs are not found. They are {not_found_uniref_id_list}")


"""
input_fas_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/FRC/frc_blast_iteration_final_set/frc_target_UniRef100.fas"
output_fas_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/FRC/frc_blast_iteration_final_set/UniRef100_under_UniRef90_O87838.fas"
subset_fas_given_uniref100(input_fas_path, output_fas_path, list(largest_Uniref90["UniRef100_ID"]))
"""


def blast_to_self(fas_path, num_alignments=100):
    # move adding blast to environment to parent scripts "recluster_species_UniRef"
    # Blast_bin_path = "/gpfs/data/lilab/home/zhoub03/software/blast/ncbi-blast-2.16.0+/bin"
    # if Blast_bin_path not in os.environ["PATH"]:
    #     os.environ["PATH"] = os.environ["PATH"] + ":" + Blast_bin_path
    # fas_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/FRC/frc_blast_iteration_final_set/UniRef100_under_UniRef90_O87838.fas"
    os.system(f"makeblastdb -in {fas_path} -dbtype prot -parse_seqids")
    fas_prefix = ".".join(os.path.split(fas_path)[-1].split(".")[:-1])
    fas_dir = os.path.split(fas_path)[0]
    blast_res_path = os.path.join(fas_dir, f'{fas_prefix}_blast_self.txt')
    blast_command = f"blastp -query {fas_path} -db {fas_path} -outfmt 6 -num_alignments {num_alignments} > {blast_res_path}"
    os.system(blast_command)
    return blast_res_path

