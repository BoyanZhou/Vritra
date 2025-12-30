"""
module add python/cpu/3.6.5

Script Name: get_centroid.py
Author: Boyan Zhou
Date: 2025-02-02
Description: [Get the centroid seq for given seqs; used in "recluster_species_UniRef.py"]

Function Overview:
    1. centroid_seq():
        - Purpose: Get the centroid (UniRef100_ID), need to calculate the pairwise identities
        - Dependencies: identity_within_group

    2. centroid_seq_given_identity():
        - Purpose: Get the centroid (UniRef100_ID) given the pairwise identities
        - Dependencies: None
"""

import os
import pandas as pd
from Bio import SeqIO
import recluster_species_UniRef.identity_within_group as iwg
import sys


# key function1
def centroid_seq(seq_record_dict, cluster_name, intermediate_folder):
    """
    :param seq_record_dict: {uniref100_id: seq_record_by_Bio}
    :param cluster_name: only for temp filename, can be a UniRef90
    :param intermediate_folder: intermediate result by Blast
    :return: the single UniRef100 ID of the centroid sequence
    """
    # Pre-A: whether only one sequence, no need to calculate pairwise identity
    if len(seq_record_dict) == 0:
        raise ValueError("The input of centroid_seq is empty")
    if len(seq_record_dict) == 1:
        return list(seq_record_dict.keys())[0]

    # A: Get the pairwise identity pandas data frame
    if len(seq_record_dict) < 11:
        # 1. using Biopython, data frame with three columns, ['Seq1', 'Seq2', 'Identity']
        identity_res_df = iwg.seqs_pairwise_identity(seq_record_dict, symmetric=True)
    else:
        # 2. use Blast
        cluster_fas_path = os.path.join(intermediate_folder, f"{cluster_name}.fas")
        # # Subtract all sequences of the cluster to fasta file
        with open(cluster_fas_path, "w") as cluster_f:
            for uniref100_id, seq_record in seq_record_dict.items():
                SeqIO.write(seq_record, cluster_f, "fasta")
        # # Blast to it self to get pairwise identity file
        cluster_blast_res_path = iwg.blast_to_self(cluster_fas_path, num_alignments=len(seq_record_dict))
        identity_res_df = pd.read_table(cluster_blast_res_path, header=None)
        identity_res_df = identity_res_df.iloc[:, 0:3]
        identity_res_df.columns = ['Seq1', 'Seq2', 'Identity']
    # B: Get the centroid seq, highest identity to other sequences
    uniref100_mean_iden_dict = {}       # the average identity to other seq within the cluster
    for uniref100_id in seq_record_dict.keys():
        uniref100_mean_iden_dict.update({uniref100_id: identity_res_df.loc[identity_res_df['Seq1'] == uniref100_id]['Identity'].mean(skipna=True)})
    highest_iden = max(uniref100_mean_iden_dict.values())
    uniref100_id_highest_iden_list = [key for key, value in uniref100_mean_iden_dict.items() if value == highest_iden]
    # C: if there multiple seq with the same highest identity to other seqs, choose the longest
    if len(uniref100_id_highest_iden_list) == 1:
        centroid_uniref100 = uniref100_id_highest_iden_list[0]
    elif len(uniref100_id_highest_iden_list) > 1:
        uniref100_id_length_dict = {i: len(seq_record_dict[i].seq) for i in uniref100_id_highest_iden_list}
        centroid_uniref100 = max(uniref100_id_length_dict, key=uniref100_id_length_dict.get)
    else:
        print(f"In get_centroid, uniref100_mean_iden_dict is {uniref100_mean_iden_dict}, identity_res_df is {identity_res_df},"
              f"seq_record_dict is {seq_record_dict}")
        raise ValueError("Invalid value provided")
    return centroid_uniref100


def centroid_seq_given_identity(seq_record_dict, identity_res_df):
    """
    Function just like "centroid_seq", but the pairwise identity has been calculated
    :param seq_record_dict: {uniref100_id: seq_record_by_Bio}
    :param identity_res_df: intermediate result by Blast
    :return: the single UniRef100 ID of the centroid sequence
    """
    # B: Get the centroid seq, highest identity to other sequences
    uniref100_mean_iden_dict = {}       # the average identity to other seq within the cluster
    for uniref100_id in seq_record_dict.keys():
        uniref100_mean_iden_dict.update({uniref100_id: identity_res_df.loc[identity_res_df['Seq1'] == uniref100_id]['Identity'].mean(skipna=True)})
    highest_iden = max(uniref100_mean_iden_dict.values())
    uniref100_id_highest_iden_list = [key for key, value in uniref100_mean_iden_dict.items() if value == highest_iden]
    # C: if there multiple seq with the same highest identity to other seqs, choose the longest
    if len(uniref100_id_highest_iden_list) == 1:
        centroid_uniref100 = uniref100_id_highest_iden_list[0]
    elif len(uniref100_id_highest_iden_list) > 1:
        uniref100_id_length_dict = {i: len(seq_record_dict[i].seq) for i in uniref100_id_highest_iden_list}
        centroid_uniref100 = max(uniref100_id_length_dict, key=uniref100_id_length_dict.get)
    else:
        print(f"In centroid_seq_given_identity, the identity_res_df is {identity_res_df}")
        raise ValueError("Invalid value provided")
    return centroid_uniref100

