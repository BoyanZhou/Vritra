"""
module add python/cpu/3.6.5

Script Name: clarify_species_different_uniref90.py
Author: Boyan Zhou
Date: 2025-02-02
Description: [For MU_MS, determine whether the seqs annotated as the same species from different UniRef90s
              should be merged.]

Function Overview:
    1. clarify_same_species_different_uniref90():
        - Purpose: main function
        - Dependencies: merge_overlapping_sets
    2. preprocess_data(data):
        - Purpose: Merges overlapping sets from a list of sets. For merging different UniRef100's UniRef90
        - Dependencies: None
"""

import os
import networkx as nx
import pandas as pd
import recluster_species_UniRef.identity_within_group as iwg
import recluster_species_UniRef.affinity_propagation as aff_p
import recluster_species_UniRef.get_centroid as get_centroid


###############################
# for situation of MU_MS_list #
###############################
def clarify_same_species_different_uniref90(corresponding_df_temp, seq_record_total, intermediate_res_dir,
                                            merge_iden_threshold=95, species_iden_threshold=97):
    """
    Identify whether species in different UniRef90 should be merged, after clarifying UniRef90 should be merged,
    MU_MS_list can be transformed into OneU_MS_list or other groups
    :param corresponding_df_temp: corresponding_df of current multiple uniref90
    :param seq_record_total: seq_record_total = UniRefSpecies_frc.record_dict
    :return: {species_name: within_species_centroid_UniRef100s_dict}
    """
    total_res_dict = {}
    for species_name in corresponding_df_temp["species"].unique():
        print(f"\nclarify_same_species_different_uniref90 for species {species_name} ... ...")
        within_species_df = corresponding_df_temp.loc[corresponding_df_temp["species"] == species_name]
        # debug start
        # duplicate_ids = within_species_df[within_species_df.duplicated(subset="UniRef100_ID", keep=False)]
        # Print the rows with duplicate UniRef100_ID
        # print(duplicate_ids)
        # debug end

        # which uniref90s this species covers
        unique_uniref90_list = within_species_df["UniRef90_ID"].unique()
        print(f"Unique_uniref90_list is {unique_uniref90_list}")
        if len(unique_uniref90_list) > 1:
            # same species name across different UniRef90 with n sequences
            seq_n_within_species = within_species_df.shape[0]
            # only show the first 2 uniref90 names in the temp fas filename, the species name in the  temp fas filename should not contain " ", "(", ")"
            species_name_standard = species_name.replace(' ', '_').replace('(', '_').replace(')', '_')
            output_temp_fas_path = os.path.join(intermediate_res_dir, f"{species_name_standard}_{'_'.join(unique_uniref90_list[:2])}_etal.fas")
            target_uniref_id_list = within_species_df["UniRef100_ID"]
            print(f"For species {species_name}, the length of target_uniref_id_list is {len(target_uniref_id_list)}.")
            # pairwise_iden_df.columns = ['Seq1', 'Seq2', 'Identity'], then add ['Seq1_UniRef90', 'Seq2_UniRef90']
            # 2025/02/19, target_uniref_id_list 需要在 seq_record_total 中
            pairwise_iden_df = iwg.pairwise_identity_by_two_method(seq_record_total, output_temp_fas_path, target_uniref_id_list, 20, True)
            # # annotate the UniRef90 for each 'Seq1', 'Seq2' (UniRef100)
            pairwise_iden_df['Seq1_UniRef90'] = pairwise_iden_df["Seq1"].map(within_species_df.set_index("UniRef100_ID")["UniRef90_ID"])
            pairwise_iden_df['Seq2_UniRef90'] = pairwise_iden_df["Seq2"].map(within_species_df.set_index("UniRef100_ID")["UniRef90_ID"])
            # # check the identity between seqs from different UniRef90
            pairwise_iden_df_diff_uniref90 = pairwise_iden_df.loc[(pairwise_iden_df["Seq1_UniRef90"] != pairwise_iden_df["Seq2_UniRef90"]) & (pairwise_iden_df["Identity"] >=merge_iden_threshold)]
            if pairwise_iden_df_diff_uniref90.shape[0] > 0:
                """
                pairwise_iden_df_diff_uniref90 = {'Seq1': ['A', 'B', 'C', 'B', 'E'], 'Seq2': ['B', 'C', 'D', 'D', 'F'], 'Identity': [97, 95, 96, 97, 99]}
                pairwise_iden_df_diff_uniref90 = pd.DataFrame(pairwise_iden_df_diff_uniref90)
                """
                print(f"For species {species_name}, there are seqs assigned to different UniRef90s should be combined. \n {pairwise_iden_df_diff_uniref90}")
                # uniref100_90_dict = {i: j for i, j in zip(within_species_df["UniRef100_ID"], within_species_df["UniRef90_ID"])}
                # Seqs assigned to different UniRef90s need to be analyzed together
                G1 = nx.Graph()
                # G1.add_edges_from(pairwise_iden_df_diff_uniref90.values)
                G1.add_edges_from(pairwise_iden_df_diff_uniref90.iloc[:, :2].itertuples(index=False, name=None))
                # Find connected components, each element is a set of UniRef100 {'A', 'B', 'C'}
                # each connected_component means there are close UniRef100s in different UniRef90s
                connected_components = list(nx.connected_components(G1))
                uniref90_in_component_list = []     # list of UniRef90 sets
                for connected_component in connected_components:
                    # for each connected_component, find the corresponding UniRef90s of each seq
                    uniref90_in_component = set(pd.Series(list(connected_component)).map(within_species_df.set_index("UniRef100_ID")["UniRef90_ID"]))
                    uniref90_in_component_list.append(uniref90_in_component)
                # Need to determine whether some component should be merged (i.e. with shared UniRef90) in analysis
                # For example, {2, 3, 4} and {4, 5} should be merged
                uniref90_in_component_merged_list = merge_overlapping_sets(uniref90_in_component_list)
                # For each separated component, affinity propagation to find core UniRef100, using precalculated indentity
                print(f"Total uniref90_in_component_merged_list is {uniref90_in_component_merged_list}")

                within_species_centroid_UniRef100s_dict = {}    # the results to be recoreded, {centroid_uniref100: list of represented uniref100s}
                for uniref90_in_component_merged in uniref90_in_component_merged_list:
                    # uniref90_in_component_merged is a list of uniref90s
                    print(f"uniref90_in_component_merged is {uniref90_in_component_merged}")
                    # the seqs from given species from these uniref90s should analyzed together
                    pairwise_iden_df_component_merged = pairwise_iden_df.loc[(pairwise_iden_df["Identity"]>merge_iden_threshold) & pairwise_iden_df['Seq1_UniRef90'].isin(uniref90_in_component_merged) & pairwise_iden_df['Seq2_UniRef90'].isin(uniref90_in_component_merged)]
                    print(f"The number of uniref100s for aff is {len(set(pairwise_iden_df_component_merged['Seq1']) | set(pairwise_iden_df_component_merged['Seq2']))}")
                    #######################################################
                    # get the clustered results from affinity propagation #
                    #######################################################
                    # {"cluster_0": [UniRef100_1, UniRef100_2]},
                    aff_p_cluster_dict = aff_p.affinity_propagation_clustering(pairwise_iden_df_component_merged)
                    centroid_cluster_dict = {}      # {centroid_uniref100: "cluster_0"}
                    centroid_n_seq_dict = {}        # {centroid_uniref100: {"n_seq": 20}}
                    # get the centroid seq for each cluster by affinity propagation
                    for cluster_name, clustered_UniRef100 in aff_p_cluster_dict.items():
                        cluster_seq_record_dict = {i: seq_record_total[i] for i in clustered_UniRef100}
                        cluster_name_temp = f"{cluster_name}_{clustered_UniRef100[0]}"
                        # GET centroid seq for the cluster
                        cluster_centroid_uniref100 = get_centroid.centroid_seq(cluster_seq_record_dict, cluster_name_temp, intermediate_res_dir)
                        centroid_cluster_dict.update({cluster_centroid_uniref100: cluster_name})
                        centroid_n_seq_dict.update({cluster_centroid_uniref100: {"n_seq": len(clustered_UniRef100)}})
                    print(f"centroid_n_seq_dict is {centroid_n_seq_dict}")
                    #####################################
                    # merge can be merged centroid seqs #
                    #####################################
                    # Merge if two centroids are too close to each other
                    final_centroid_UniRef100s_dict = {}     # {centroid_uniref100: all represented UniRef100s}
                    if len(centroid_n_seq_dict) > 1:
                        # get the pairwise identity df only for centroid
                        centroid_pairwise_identity_df = pairwise_iden_df.loc[pairwise_iden_df['Seq1_UniRef90'].isin(centroid_n_seq_dict.keys())
                                                                             & pairwise_iden_df['Seq2_UniRef90'].isin(centroid_n_seq_dict.keys())]
                        # {0: [unref100_A, unref100_B], 1: [unref100_C]}, the first is the representative one
                        group_centroid_list_dict = aff_p.greedy_affinity_clustering(centroid_pairwise_identity_df, centroid_n_seq_dict, species_iden_threshold)
                        print(f"The group_centroid_list_dict is {group_centroid_list_dict}")
                        for key, value in group_centroid_list_dict.items():
                            # value = [unref100_A, unref100_B] or [unref100_C]
                            if len(value) == 1:
                                final_centroid_UniRef100s_dict.update({value[0]: aff_p_cluster_dict[centroid_cluster_dict[value[0]]]})
                            else:
                                uniref100_list_temp = []
                                for centroid_temp in value:
                                    uniref100_list_temp.extend(aff_p_cluster_dict[centroid_cluster_dict[centroid_temp]])
                                final_centroid_UniRef100s_dict.update({value[0]: uniref100_list_temp})
                    elif len(centroid_n_seq_dict) == 1:
                        print(f"Only one centroid, no need to merge for {centroid_n_seq_dict}")
                        centroid_temp = list(centroid_cluster_dict.keys())[0]
                        final_centroid_UniRef100s_dict.update({centroid_temp: aff_p_cluster_dict[centroid_cluster_dict[centroid_temp]]})
                    ##################################################
                    # Extract cluster overlapping multiple UniRef90s #
                    ##################################################
                    for key, value in final_centroid_UniRef100s_dict.items():
                        # key is the centroid seq | value is the represented seqs
                        overlapped_uniref90 = set(pd.Series(value).map(within_species_df.set_index("UniRef100_ID")["UniRef90_ID"]))
                        if len(overlapped_uniref90) > 1:
                            within_species_centroid_UniRef100s_dict.update({key: value})
                # store the result to total_res_dict
                print(f"within_species_centroid_UniRef100s_dict is {within_species_centroid_UniRef100s_dict}")
                if len(within_species_centroid_UniRef100s_dict) > 0:
                    total_res_dict.update({species_name: within_species_centroid_UniRef100s_dict})
    return total_res_dict


def merge_overlapping_sets(set_list):
    """
    Merges overlapping sets from a list of sets. For merging different UniRef100's UniRef90
    :param set_list: A list containing multiple UniRef90s
    :return: A list of merged sets with all overlaps resolved.
    """
    merged = []  # List to store merged sets
    while set_list:
        current_set = set_list.pop(0)
        last_loop_set_len = -1
        while len(current_set) > last_loop_set_len:
            print(current_set)
            last_loop_set_len = len(current_set)
            # while the current_set still take in last loop, search for new overlapped set
            del_index_list = [set_index for set_index, set_temp in enumerate(set_list) if not set_temp.isdisjoint(current_set)]
            if len(del_index_list) > 0:
                print(f"Delete index is {sorted(del_index_list, reverse=True)}")
                for del_index in sorted(del_index_list, reverse=True):
                    current_set.update(set_list[del_index])
                    del set_list[del_index]
        merged.append(current_set)
    print(merged)
    return merged
