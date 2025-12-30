"""
Re-cluster the UniRef100s given the fas file & corresponding taxonomy csv file

module add python/cpu/3.6.5

Script Name: example_script.py
Author: Boyan Zhou
Date: 2025-02-02
Description: [Brief description of what the script does.]

Function Overview:
    1. load_data(file_path):
        - Purpose: Loads a dataset from a specified file path.
        - Dependencies: None
    2. preprocess_data(data):
        - Purpose: Cleans and preprocesses the input data.
        - Dependencies: None
    3. analyze_data(data):
        - Purpose: Performs analysis on the cleaned data.
        - Dependencies: preprocess_data
"""


import os
import json
# import networkx as nx
import pandas as pd
from Bio import SeqIO
import log_setup
import argparse

import recluster_species_UniRef
# import identity_within_group as iwg
# import affinity_propagation as aff_p
# import get_centroid
# import clarify_species_different_uniref90 as csdu
# import uniref100_represent_info as uni_rep_info


class UniRefSpecies:
    def __init__(self, uniref_taxonomy_path, uniref100_fas_path, logger):
        self.UniRef100_taxonomy_path = uniref_taxonomy_path
        self.record_dict = SeqIO.to_dict(SeqIO.parse(uniref100_fas_path, "fasta"))
        # ATTENTION!!!: make sure self.corresponding_df overlap with csv file of uniref_taxonomy
        corresponding_df = read_my_csv(uniref_taxonomy_path)
        """
        {'UniRef90_ID': 'UniRef90_A0A059DQ18', 'UniRef100_ID': 'UniRef100_A0A059DQ18', 'Protein_length': 489, 'Taxon_ID': '1280944', 'Taxon_ranks': 'd__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Hyphomonadales|f__Hyphomonadaceae|g__Hyphomonas|s__Hyphomonas sp. CY54-11-8',
         'superkingdom': 'Bacteria', 'phylum': 'Pseudomonadota', 'class': 'Alphaproteobacteria', 'order': 'Hyphomonadales', 'family': 'Hyphomonadaceae', 'genus': 'Hyphomonas', 'species': 'Hyphomonas sp. CY54-11-8'}
        """
        # drop NaN in "Protein_length" and convert it into int
        corresponding_df = corresponding_df.dropna(subset=["Protein_length"])
        corresponding_df["Protein_length"] = corresponding_df["Protein_length"].astype(int)
        # Keep only the row with the largest "Protein_length" for each "UniRef100_ID"
        corresponding_df_filtered = corresponding_df.loc[corresponding_df.groupby("UniRef100_ID")["Protein_length"].idxmax()]
        # Reset index (optional)
        corresponding_df_filtered = corresponding_df_filtered.reset_index(drop=True)
        self.corresponding_df = corresponding_df_filtered.loc[corresponding_df_filtered["UniRef100_ID"].isin(self.record_dict.keys())]

        self.logger = logger

        # at species OR no info at species level
        self.corresponding_df_species = self.corresponding_df.loc[self.corresponding_df['species'] != ""]
        self.corresponding_df_no_species = self.corresponding_df.loc[self.corresponding_df['species'] == ""]
        self.logger.info(f"The number of UniRef100 with species rank is {self.corresponding_df_species.shape[0]}; "
                         f"the number of UniRef100 without species rank is {self.corresponding_df_no_species.shape[0]}")
        # get the number of taxon ranks != "" for each species
        species_unique_df = self.corresponding_df_species.drop_duplicates(subset=['species'], keep='first')
        species_unique_df = species_unique_df.copy()
        species_unique_df['ranks_n'] = species_unique_df[['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']].apply(lambda row: sum(row != ""), axis=1)
        self.species_ranks_count_dict = species_unique_df.set_index('species')['ranks_n'].to_dict()
        # categorize_four_groups
        self.logger.info(f"Preliminary categorization of UniRef100s into four groups.")
        self.OneU_OneS_list, self.MU_OneS_list, self.OneU_MS_list, self.MU_MS_list = (
            categorize_uniref_species(self.corresponding_df_species))
        ##########################
        # result storage in dict #
        ##########################
        # {representative_uniref100: [all represented uniref100s]}
        self.uniref100_from_single_species = {}     # one uniref100 represents uniref100s from one species
        self.uniref100_from_multiple_species = {}   # one uniref100 represents uniref100s from multiple species

    def complete_process(self, gene_prefix, intermediate_folder, final_database_folder):
        # 1. Must process MU_MS first, after removing connected species from different UniRef90s
        #       process the left as OneU_OneS & MU_OneS (built in MU_MS_process())
        self.logger.info(f"STEP1: MU_MS_process ... ...")
        self.MU_MS_process(intermediate_folder)
        # 2. One UniRef90 contains Multiple Species
        self.logger.info(f"STEP2: OneU_MS_process ... ...")
        for uniref90_id_temp in self.OneU_MS_list:
            corresponding_df_temp = self.corresponding_df_species.loc[self.corresponding_df_species["UniRef90_ID"] == uniref90_id_temp]
            centroid_represent_uniref100_dict_temp = self.__OneU_MS_process(uniref90_id_temp, corresponding_df_temp, intermediate_folder, 16)
            self.uniref100_from_multiple_species.update(centroid_represent_uniref100_dict_temp)
        # 3. Multiple UniRef90 from One Species
        self.logger.info(f"STEP3: MU_OneS_process ... ...")
        self.MU_OneS_process(intermediate_folder)
        # 4. One UniRef90 to One Species
        self.logger.info(f"STEP4: OneU_OneS_process ... ...")
        self.OneU_OneS_process(intermediate_folder)
        # 5. Output the core target representative sequences database
        core_target_uniref100_set = self.uniref100_from_single_species.keys() | self.uniref100_from_multiple_species.keys()
        if not os.path.exists(final_database_folder):
            os.system(f"mkdir -p {final_database_folder}")
        core_target_uniref100_fas_path = os.path.join(final_database_folder, f"{gene_prefix}_core_representative_UniRef100.fas")
        core_target_uniref100_fas_f = open(core_target_uniref100_fas_path, "w")
        for uniref100_temp in core_target_uniref100_set:
            SeqIO.write(self.record_dict[uniref100_temp], core_target_uniref100_fas_f, "fasta")
        core_target_uniref100_fas_f.close()
        # 6. Output UniRef100 correspondence file (one UniRef100 for a cluster, with represented UniRef information)
        # core_uniref100_df added four columns: 'rep_species' (string); 'represented_uniref100' (list);
        # 'corresponding_uniref90' (list); 'corresponding_species' (list)
        core_uniref100_df = self.add_uniref100_rep_info()       # main results to output
        # output CSV core uniref100 correspondence
        output_core_uniref100_path = os.path.join(final_database_folder, f"{gene_prefix}_core_uniref100_correspondence.csv")
        core_uniref100_df.to_csv(output_core_uniref100_path, index=False, header=True)
        # Convert DataFrame to JSON (as a list of records)
        output_json_path = os.path.join(final_database_folder, f"{gene_prefix}_UniRef100_correspondence.json")
        corresponding_df_json = self.corresponding_df.to_dict(orient="records")
        core_uniref100_df_json = core_uniref100_df.to_dict(orient="records")

        # Combine both DataFrame JSON and dictionary into a single structure
        combined_data = {
            "corresponding_df": corresponding_df_json, "core_uniref100_df": core_uniref100_df_json
        }
        # Save to a JSON file
        with open(output_json_path, "w") as json_file:
            json.dump(combined_data, json_file, indent=4)
        # Output final results to json
        # self.output_representative_and_represented(output_json_path)

    def add_uniref100_rep_info(self):
        """
        :return: Add four columns to the original self.corresponding_df to get new core_uniref_df
        """
        # 1. merge all rep seqs (注意这里都是species有值的，不是"")
        uniref100_rep_represented_dict = self.uniref100_from_single_species.copy()
        uniref100_rep_represented_dict.update(self.uniref100_from_multiple_species)
        # 2. build the correspondence between {UniRef100: UniRef90}, and {UniRef100, species}
        uniref100_90_dict = self.corresponding_df.set_index('UniRef100_ID')['UniRef90_ID'].to_dict()
        uniref100_species_dict = self.corresponding_df.set_index('UniRef100_ID')['species'].to_dict()
        # 3. get represented species
        uniref100_rep_info_dict = {}
        for uniref100_temp, uniref100_represented_list in uniref100_rep_represented_dict.items():
            # contain info: rep_species; represented_uniref100_list; corresponding_uniref90_list; corresponding_species_list
            uniref100_info = recluster_species_UniRef.uniref100_represent_info.UniRepInfo(uniref100_temp, uniref100_represented_list, uniref100_90_dict, uniref100_species_dict)
            uniref100_info.get_represented_species(self.species_ranks_count_dict)
            uniref100_rep_info_dict.update({uniref100_temp: uniref100_info})
        # 4. add represented info: rep_species, represented_uniref100_list, corresponding_uniref90_list, corresponding_species_list
        core_uniref_df = self.corresponding_df.loc[self.corresponding_df['UniRef100_ID'].isin(uniref100_rep_info_dict.keys())]
        core_uniref_df = core_uniref_df.copy()
        core_uniref_df['rep_species'] = core_uniref_df['UniRef100_ID'].map(lambda x: uniref100_rep_info_dict[x].rep_species)
        core_uniref_df['rep_species_summary'] = core_uniref_df['UniRef100_ID'].map(lambda x: uniref100_rep_info_dict[x].corresponding_species_count_summary)
        core_uniref_df['represented_uniref100'] = core_uniref_df['UniRef100_ID'].map(lambda x: ";".join(uniref100_rep_info_dict[x].represented_uniref100_list))
        core_uniref_df['corresponding_uniref90'] = core_uniref_df['UniRef100_ID'].map(lambda x: ";".join(uniref100_rep_info_dict[x].corresponding_uniref90_list))
        core_uniref_df['corresponding_species'] = core_uniref_df['UniRef100_ID'].map(lambda x: ";".join(uniref100_rep_info_dict[x].corresponding_species_list))
        return core_uniref_df

    # !!! currently deprecated !!!, replaced by "add_uniref100_rep_info"
    def output_representative_and_represented(self, output_json_path):
        output_json_dir = os.path.split(output_json_path)[0]
        if not os.path.exists(output_json_dir):
            os.system(f"mkdir -p {output_json_dir}")
        # Convert DataFrame to JSON (as a list of records)
        corresponding_df_species_json = self.corresponding_df_species.to_dict(orient="records")

        # Combine both DataFrame JSON and dictionary into a single structure
        combined_data = {
            "corresponding_df_species": corresponding_df_species_json,
            "uniref100_from_single_species": self.uniref100_from_single_species,
            "uniref100_from_multiple_species": self.uniref100_from_multiple_species
        }
        # Save to a JSON file
        with open(output_json_path, "w") as json_file:
            json.dump(combined_data, json_file, indent=4)

    def MU_MS_process(self, intermediate_folder):
        if not os.path.exists(intermediate_folder):
            os.system(f"mkdir -p {intermediate_folder}")
        # processing record for MU_MS, finally recorded in logger, [[dig_out_species_uniref_dict, uniref90_list_temp], [...]]
        merged_species_from_diff_uniref90_in_MUMS = []
        # In each group, multiple UniRef90s assigned to multiple groups, but this may be caused by incorrectly
        # assignment of species; this is, seqs from one species split into two UniRef90s due to incorrect representative
        # Thus, in the function, we should break the UniRef90s connected in this way
        for uniref90_list_temp in self.MU_MS_list:
            self.logger.info(f"\nProcessing MutipleUniRef90-MultipleSpecies for uniref90_list {uniref90_list_temp}.")
            print(f"\nProcessing MutipleUniRef90-MultipleSpecies for uniref90_list {uniref90_list_temp}.")
            # uniref90_list_temp is a group of MU-to-MS isolated with others
            corresponding_df_temp = self.corresponding_df_species.loc[self.corresponding_df_species['UniRef90_ID'].isin(uniref90_list_temp)]
            seq_record_within_MU = {i: self.record_dict[i] for i in corresponding_df_temp["UniRef100_ID"] if i in self.record_dict}
            ###############################################################
            # Dig out UniRef100s from different groups that can be merged #
            ###############################################################
            # {species: {centroid: [uniref100, ... ...]}}
            dig_out_species_uniref_dict = recluster_species_UniRef.clarify_species_different_uniref90.clarify_same_species_different_uniref90(corresponding_df_temp, seq_record_within_MU, intermediate_folder, merge_iden_threshold=95)
            self.logger.info(f"The UniRef100s need to be dig out is {dig_out_species_uniref_dict}.")
            uniref100_to_remove = []
            if len(dig_out_species_uniref_dict) > 0:
                merged_species_from_diff_uniref90_in_MUMS.append([dig_out_species_uniref_dict, uniref90_list_temp])
            for key, value in dig_out_species_uniref_dict.items():
                # record the dig out results
                self.uniref100_from_single_species.update(value)
                for key1, value1 in value.items():
                    uniref100_to_remove.extend(value1)
            # get the data after digging out some merged uniref100s
            ################################################
            # the remaining should be OneU_MS or OneU_OneS #
            ################################################
            corresponding_df_temp_remain = corresponding_df_temp[~corresponding_df_temp["UniRef100_ID"].isin(set(uniref100_to_remove))]
            for uniref90_temp in corresponding_df_temp_remain["UniRef90_ID"].unique():
                corresponding_df_temp_remain_uniref90 = corresponding_df_temp_remain[corresponding_df_temp_remain["UniRef90_ID"] == uniref90_temp]
                # IF ONLY ONE SPECIES CONTAINED (OneU_OneS)
                if len(set(corresponding_df_temp_remain_uniref90['species'])) == 1:
                    self.logger.info(f"OneU_OneS processing after digging out, for {uniref90_temp} ... ...")
                    seq_record_dict_uniref90_temp = {i: self.record_dict[i] for i in corresponding_df_temp_remain_uniref90["UniRef100_ID"]}
                    centroid_seq_uniref100_id = recluster_species_UniRef.get_centroid.centroid_seq(seq_record_dict_uniref90_temp, uniref90_temp, intermediate_folder)
                    self.uniref100_from_single_species.update({centroid_seq_uniref100_id: list(seq_record_dict_uniref90_temp.keys())})
                # IF MULTIPLE SPECIES (OneU_MS)
                elif len(set(corresponding_df_temp_remain_uniref90['species'])) > 1:
                    self.logger.info(f"OneU_MultipleSpecies processing after digging out, for {uniref90_temp} ... ...")
                    centroid_represent_uniref100_dict_temp = self.__OneU_MS_process(uniref90_temp,
                                                                                    corresponding_df_temp_remain_uniref90,
                                                                                    intermediate_folder)
                    self.uniref100_from_multiple_species.update(centroid_represent_uniref100_dict_temp)
        self.logger.info(f"\nFor the summary of MU_MS_process, there are {len(merged_species_from_diff_uniref90_in_MUMS)} "
                         f"examples of merged species from different UniRef90s.\n")
        for i_temp in merged_species_from_diff_uniref90_in_MUMS:
            self.logger.info(f"Merged species in MU_MS_process, {i_temp}.")

    def __OneU_MS_process(self, uniref90_id, corresponding_df_temp, intermediate_folder, n_seq_threshold=16):
        """
        For single UniRef90 with target corresponding_df_temp (has been filtered)
        Attention! uniref90_id may not be real UniRef90, may be modified UniRef90 (some sequences may be dig out)
        """
        if not os.path.exists(intermediate_folder):
            os.system(f"mkdir -p {intermediate_folder}")
        self.logger.info(f"Processing OneUniRef90-MultipleSpecies for UniRef90 (not mean perfect match) {uniref90_id}.")
        # self.logger.info(f"For debug, It's corresponding_df_temp is {corresponding_df_temp}")
        ###############################
        # calculate pairwise identity #
        ###############################
        subtracted_uniref90_fas_path = os.path.join(intermediate_folder, f"{uniref90_id}.fas")
        target_uniref_id_list = corresponding_df_temp["UniRef100_ID"]
        identity_res_df = recluster_species_UniRef.identity_within_group.pairwise_identity_by_two_method(self.record_dict, subtracted_uniref90_fas_path, target_uniref_id_list, n_seq_threshold, True)
        ########################
        # Affinity Propagation #
        ########################
        # cluster_res_dict is like {'cluster_0': [UniRef100_a, ...]}
        cluster_res_dict = recluster_species_UniRef.affinity_propagation.affinity_propagation_clustering(identity_res_df)
        centroid_N_dict = {}        # {UniRef100_a: {'n_seq': number_represented_seqs}}
        centroid_cluster_dict = {}  # {UniRef100_a: 'cluster_0'}
        # self.logger.info(f"For debug, cluster_res_dict is {cluster_res_dict}; identity_res_df is {identity_res_df}")
        for key, uniref100_list in cluster_res_dict.items():
            if len(uniref100_list) == 1:
                # Attention! If only one UniRef100 in the cluster, will get empty for cluster_identity_res_df
                centroid_cluster_dict.update({uniref100_list[0]: key})
                centroid_N_dict.update({uniref100_list[0]: {'n_seq': 1}})
            else:
                # more than one UniRef100 in the cluster, find the centroid UniRef100
                cluster_identity_res_df = identity_res_df.loc[identity_res_df["Seq1"].isin(uniref100_list) & identity_res_df["Seq2"].isin(uniref100_list)]
                cluster_seq_record_dict = {i: self.record_dict[i] for i in uniref100_list}
                cluster_centroid = recluster_species_UniRef.get_centroid.centroid_seq_given_identity(cluster_seq_record_dict, cluster_identity_res_df)
                centroid_cluster_dict.update({cluster_centroid: key})
                centroid_N_dict.update({cluster_centroid: {'n_seq': len(uniref100_list)}})

        ##########################
        # Greedy merge centroids #
        ##########################
        # centroids can also be merged if they are too close to each other
        centroid_list = list(centroid_cluster_dict.keys())
        centroid_uniref100_n_seq_dict = {i: {"length_seq": len(self.record_dict[i].seq)} for i in centroid_list}
        self.logger.info(f"After affinity propagation, centroid_cluster_dict is {centroid_cluster_dict}.\n"
                         f"The sequence length of centroid in each cluster is {centroid_uniref100_n_seq_dict}.\n"
                         f"The number of represented seqs by each centroid is {centroid_N_dict}")
        cluster_centroid_identity_res_df = identity_res_df.loc[
            identity_res_df["Seq1"].isin(centroid_list) & identity_res_df["Seq2"].isin(centroid_list)]
        # {0: [unref100_A, unref100_B], 1: [unref100_C]}
        grouped_centroid_dict = recluster_species_UniRef.affinity_propagation.greedy_affinity_clustering(cluster_centroid_identity_res_df, centroid_N_dict, 95)
        self.logger.info(f"After greedy merge centroids, grouped_centroid_dict is {grouped_centroid_dict}.")

        #######################################
        # Centroid Represent Uniref100 Series #
        #######################################
        centroid_represent_uniref100_dict = {}      # {representative_centroid: pd.Series([UniRef100_a, ...])}
        for key, centroids_temp in grouped_centroid_dict.items():
            uniref100_list = []
            for uniref100_temp in centroids_temp:
                uniref100_list.extend(cluster_res_dict[centroid_cluster_dict[uniref100_temp]])
            # for each group by "Greedy merge centroids", choose the centroid with largest number of seqs
            centroid_max_seq = max(centroids_temp, key=lambda k: centroid_N_dict[k]['n_seq'])
            centroid_represent_uniref100_dict.update({centroid_max_seq: uniref100_list})
        self.logger.info(f"The final number of seqs in each representative centroid is "
                         f"{[[i, len(j)] for i, j in centroid_represent_uniref100_dict.items()]}")
        return centroid_represent_uniref100_dict

    def MU_OneS_process(self, intermediate_res_dir):
        print(f"Processing multiple UniRef90 within one species ... ...")
        merged_species_from_diff_uniref90_in_MUOneS = []
        species_grouped_centroid_dict = {}
        for uniref90_list_temp in self.MU_OneS_list:
            # uniref90_list_temp is the list of uniref90 for ONE SPECIES
            corresponding_df_temp = self.corresponding_df_species.loc[self.corresponding_df_species["UniRef90_ID"].isin(uniref90_list_temp)]
            current_species = list(corresponding_df_temp["species"])[0]
            centroid_uniref100_uniref90_dict = {}   # record the centroid seq and number of seq within the uniref90, {centroid_uniref100: {"UniRef90": UniRef90_ID, "n_seq": 20}}
            # 1. centroid seq for each uniref90
            for uniref90_temp in uniref90_list_temp:
                # the corresponding uniref100 series of current uniref90
                uniref100_series = corresponding_df_temp.loc[corresponding_df_temp["UniRef90_ID"] == uniref90_temp]["UniRef100_ID"]
                seq_record_dict_uniref90_temp = {i: self.record_dict[i] for i in uniref100_series}
                centroid_seq_uniref100_id = recluster_species_UniRef.get_centroid.centroid_seq(seq_record_dict_uniref90_temp, uniref90_temp, intermediate_res_dir)
                centroid_uniref100_uniref90_dict.update({centroid_seq_uniref100_id: {"UniRef90": uniref90_temp, "n_seq": len(seq_record_dict_uniref90_temp),
                                                                                     "UniRef100_list": list(uniref100_series)}})
            # 2. check whether these centroid seq are actually close to each other
            if len(centroid_uniref100_uniref90_dict) < 11:
                seq_record_of_centroid = {i: self.record_dict[i] for i in centroid_uniref100_uniref90_dict.keys() if i in self.record_dict}
                identity_res_df = recluster_species_UniRef.identity_within_group.seqs_pairwise_identity(seq_record_of_centroid)
            else:
                # there should not be " " in the species name of fas path
                current_species_unified_name = current_species.replace(' ', '_').replace("(", "_").replace(")", "_")
                subtracted_species_fas_path = os.path.join(intermediate_res_dir, f"{current_species_unified_name}.fas")
                # # Subtract all centroid sequences within species
                recluster_species_UniRef.identity_within_group.subset_fas_given_uniref100(self.record_dict, subtracted_species_fas_path, list(centroid_uniref100_uniref90_dict.keys()))
                # # Blast to it self to get pairwise identity file
                species_blast_res_path = recluster_species_UniRef.identity_within_group.blast_to_self(subtracted_species_fas_path, num_alignments=len(centroid_uniref100_uniref90_dict))
                identity_res_df = pd.read_table(species_blast_res_path, header=None)
                identity_res_df = identity_res_df.iloc[:, 0:3]
                identity_res_df.columns = ['Seq1', 'Seq2', 'Identity']
                # Attention, A-B and B-A may be different in blast result
                identity_res_df = identity_res_df.loc[identity_res_df['Seq1'] != identity_res_df['Seq2']]
            # # 2.1 greedy algorithm to merge possible
            # {0: [unref100_A, unref100_B], 1: [unref100_C]},  the 1st one in the list is the selected
            grouped_centroid_dict = recluster_species_UniRef.affinity_propagation.greedy_affinity_clustering(identity_res_df, centroid_uniref100_uniref90_dict, 95)
            for key, value in grouped_centroid_dict.items():
                if len(value) > 1:
                    # Merge of centroid from different UniRef90 exists
                    self.logger.info(f"In MU_OneS_process, centroid seqs {value} from different UniRef90 were merged to {value[0]}.")
                    merged_species_from_diff_uniref90_in_MUOneS.append(value)   # for logger record
                uniref100_list_temp = []
                for centroid_temp in value:
                    uniref100_list_temp.extend(centroid_uniref100_uniref90_dict[centroid_temp]["UniRef100_list"])
                self.uniref100_from_single_species.update({value[0]: uniref100_list_temp})
            # species_grouped_centroid_dict.update({current_species: grouped_centroid_dict})
        # return species_grouped_centroid_dict
        self.logger.info(
            f"\nFor the summary of MU_OneS_process, there are {len(merged_species_from_diff_uniref90_in_MUOneS)} "
            f"examples of merged UniRef100 (each represent one UniRef90) from different UniRef90s: {merged_species_from_diff_uniref90_in_MUOneS}\n")

    def OneU_OneS_process(self, intermediate_res_dir):
        # for this group, we only need to find the centroid sequence UniRef100 ID (representative) for each UniRef90
        for uniref90_temp in self.OneU_OneS_list:
            uniref100_series = self.corresponding_df_species.loc[self.corresponding_df_species["UniRef90_ID"] == uniref90_temp]["UniRef100_ID"]
            seq_record_dict_uniref90_temp = {i: self.record_dict[i] for i in uniref100_series}
            centroid_seq_uniref100_id = recluster_species_UniRef.get_centroid.centroid_seq(seq_record_dict_uniref90_temp, uniref90_temp, intermediate_res_dir)
            self.uniref100_from_single_species.update({centroid_seq_uniref100_id: list(uniref100_series)})


"""
UniRef90_O87838 和 UniRef90_A0A917RK06 中的 Streptomyces rochei 应该是一起的！！！
corresponding_df_temp.loc[corresponding_df_temp["UniRef90_ID"] == "UniRef90_O87838"]
corresponding_df_temp.loc[corresponding_df_temp["UniRef90_ID"] == "UniRef90_A0A917RK06"]

UniRef90_A0A917RK06 和 UniRef90_A0A640UML5 中的 Streptomyces sparsogenes 应该是一起的！！！
corresponding_df_temp.loc[corresponding_df_temp["UniRef90_ID"] == "UniRef90_A0A917RK06"]
corresponding_df_temp.loc[corresponding_df_temp["UniRef90_ID"] == "UniRef90_A0A640UML5"]
"""


#################################
# read uniref taxonomy csv file #
#################################
def read_my_csv(csv_path):
    table_by_list = []
    ncol = -1
    col_name = []
    with open(csv_path) as csv_f:
        for line in csv_f:
            cols = line.strip().split(",")
            if ncol < 0:
                ncol = len(cols)
            if len(col_name) == 0:
                col_name = cols
                continue
            if len(cols) == ncol:
                table_by_list.append(cols)
    data_df = pd.DataFrame(table_by_list, columns=col_name)
    return data_df


##################
# To Four Groups #
##################
def categorize_uniref_species(input_df):
    """
    Categorize the relationship between UniRef90 and Species into four groups
    1. 1UniRef90-1Species; 2. MultiUniRef90-1Species  3. 1UniRef90_MultiSpecies; 4. MultiUniRef90-MultiSpecies
    UniRef90_ID             UniRef100_ID Protein_length  ...               family               genus                        species
    UniRef90_A0A418WVG2     UniRef100_A0A418WVG2            415  ...     Oxalobacteraceae  Noviherbaspirillum    Noviherbaspirillum cavernae
    UniRef90_A0A418WVG2     UniRef100_A0A3D3PNG4            415  ...     Oxalobacteraceae                         Oxalobacteraceae bacterium
    UniRef90_A0A969IS53     UniRef100_A0A969IS53            417  ...                                              Rhodospirillales bacterium
    :param input_df=corresponding_df_species
    :return:
    """
    uniref_id_used = set()
    OneUniRef90_OneSpecies_list = []        # list of single UniRef90
    MultiUniRef90_OneSpecies_list = []      # each element is the UniRef90 list for OneSpecies
    OneUniRef90_MultiSpecies_list = []      # list of single UniRef90
    MultiUniRef90_MultiSpecies_list = []    # each element is the UniRef90 list for a MultiSpecies group
    # input_df=corresponding_df_species
    for uniref90_id in input_df["UniRef90_ID"].unique():
        # uniref90_id = "UniRef90_A0A418WVG2"
        # df_temp is the data with same uniref90 id
        if uniref90_id in uniref_id_used:
            continue
        df_temp = input_df.loc[input_df["UniRef90_ID"] == uniref90_id]
        df_temp2 = input_df.loc[input_df["species"].isin(df_temp["species"])]
        # 1. if this uniref90 only contain 1 species
        if len(df_temp["species"].unique()) == 1:
            # 1.1 if the species only exists in one uniref90
            if len(df_temp2["UniRef90_ID"].unique()) == 1:
                OneUniRef90_OneSpecies_list.append(uniref90_id)
                uniref_id_used.add(uniref90_id)
                continue
            # 1.2 if the species exists in multiple uniref90
            else:
                df_temp3 = input_df.loc[input_df["UniRef90_ID"].isin(df_temp2["UniRef90_ID"].unique())]
                # 1.2.1 if multiple uniref90 only in one species
                if len(df_temp3["species"].unique()) == 1:
                    MultiUniRef90_OneSpecies_list.append(list(df_temp3["UniRef90_ID"].unique()))
                    uniref_id_used.update(list(df_temp3["UniRef90_ID"]))
                    continue
        # 2. if this uniref90 contain multiple species
        else:
            # 2.1 if multiple species only exists in one uniref90
            if len(df_temp2["UniRef90_ID"].unique()) == 1:
                OneUniRef90_MultiSpecies_list.append(uniref90_id)
                uniref_id_used.add(uniref90_id)
                continue
        # 3. Deal with MultiUniRef90_MultiSpecies
        while df_temp2.shape[0] > df_temp.shape[0]:
            df_temp = input_df.loc[input_df["UniRef90_ID"].isin(df_temp2["UniRef90_ID"])]
            df_temp2 = input_df.loc[input_df["species"].isin(df_temp["species"])]
        MultiUniRef90_MultiSpecies_list.append(list(df_temp["UniRef90_ID"].unique()))
        uniref_id_used.update(list(df_temp["UniRef90_ID"].unique()))
    print(f"The number of used uniref90 is {len(uniref_id_used)}. The total number of unique UniRef90 is {len(input_df['UniRef90_ID'].unique())}")
    return OneUniRef90_OneSpecies_list, MultiUniRef90_OneSpecies_list, OneUniRef90_MultiSpecies_list, MultiUniRef90_MultiSpecies_list


#####################################
# Download corresponding UniRef100s #
#####################################
def expand_uniref90_to_100(gene_prefix, core_uniref90_path, out_dir, logger):
    """

    :param gene_prefix:
    :param core_uniref90_path: previous gene_blast.pre_labeled_new_labeled_merge_path
    :param out_dir:
    :return:
    """
    # out_dir = os.path.split(gene_blast.pre_labeled_new_labeled_merge_path)[0]   # "PATH/dataset_construction_by_annotated_seqs"
    if not os.path.exists(out_dir):
        os.system(f"mkdir -p {out_dir}")

    output_json_path = os.path.join(out_dir, f"{gene_prefix}_core_uniref90_100.json")
    output_taxon_ranks_path = os.path.join(out_dir, f"{gene_prefix}_core_uniref90_100_taxon_ranks.csv")
    out_uniref100_fas_path = os.path.join(out_dir, f"{gene_prefix}_core_uniref100.fas")

    uniref90_api_response_dict, uniref100_api_response_dict = recluster_species_UniRef.download_uniref100_from_uniref90.uniref90_to_100_fas(
        core_uniref90_path, output_json_path, output_taxon_ranks_path, out_uniref100_fas_path)
    logger.info(f"For {len(uniref90_api_response_dict)} UniRef90s, {len(uniref90_api_response_dict['Success 200: Success'])} are successfully fetched.")
    logger.info(f"For {len(uniref100_api_response_dict)} UniRef100s, {len(uniref100_api_response_dict['Success 200: Success'])} are successfully fetched.")
    logger.info(uniref90_api_response_dict)
    logger.info(uniref100_api_response_dict)


################################################
# Read the uniref100 sequences using Biopython #
################################################
"""
result_file = open('test_biopython_record_dict.fasta', "w")
record_dict = SeqIO.to_dict(SeqIO.parse("/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/FRC/frc_blast_iteration_final_set/frc_target_UniRef100.fas", "fasta"))

# 看看之前第一个cluster的在不在uniref100的fasta中
sum(input_df.loc[input_df["UniRef90_ID"].isin(pd.Series(d[0]))]["UniRef100_ID"].isin(record_dict.keys()))
len(input_df.loc[input_df["UniRef90_ID"].isin(pd.Series(d[0]))]["UniRef100_ID"])
"""


def get_args():
    """Parses command-line arguments and returns them."""

    # Define parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-g", "--gene-prefix", type=str, required=True, help="Prefix for gene files")
    parent_parser.add_argument("-o", "--output-dir", type=str, required=True, help="Directory of output ")
    parent_parser.add_argument("-L", "--logger-path", type=str, required=True, help="Path to save log files")

    # Main parser
    parser = argparse.ArgumentParser(description="Expand and Re-cluster all labeled sequences (UniRef90 to UniRef100).")
    subparsers = parser.add_subparsers(dest="mode", help="Mode of operation")

    # Subparser for 'uniref90_to_100' mode
    parser_expand = subparsers.add_parser("uniref90_to_100", parents=[parent_parser],
                                         help="Expand the refined Uniref90 dataset to corresponding UniRef100s")
    parser_expand.add_argument("-c", "--core_uniref90", type=str, required=True, help="Path to the core UniRef90 constructed in last step.")

    # Subparser for 're_cluster' mode
    parser_cluster = subparsers.add_parser("re_cluster", parents=[parent_parser], help="Re-cluster")
    parser_cluster.add_argument("-t", "--uniref-taxonomy-path", type=str, required=True,
                                help="Path to the taxonomy of UniRef100, generated by mode 'uniref90_to_100'")
    parser_cluster.add_argument("-f", "--uniref100-fas-path", type=str, required=True,
                                help="Path to the FASTA file of UniRef100, generated by mode 'uniref90_to_100'")
    parser_cluster.add_argument("-n", "--intermediate-folder", type=str, required=True,
                                help="Folder to store intermediate files")
    parser_cluster.add_argument("-b", "--blast-bin-dir", type=str, required=True,
                                help="The directory path of bin of downloaded ncbi-blast. Specify this parameter if "
                                     "ncbi-blast has not been added to the environment.")

    args1 = parser.parse_args()

    if args1.mode is None:
        parser.error("A subcommand is required. Use -h for help.")
    return args1


if __name__ == "__main__":
    args = get_args()
    Logger = log_setup.setup_log_file(args.logger_path)
    if args.mode == "uniref90_to_100":
        """
        Example:
        python ${vritra_dir}recluster_species_UniRef.py ${mode} -g ${gene_name} 
        -o ${parent_folder}/dataset_construction_refined 
        -c ${parent_folder}/dataset_construction_refined/${gene_name}_refined_core_UniRef90_filtered.fas 
        -L uniref90_to_uniref100_0526.log
        Output:
        xdhA_core_uniref90_100.json
        xdhA_core_uniref90_100_taxon_ranks.csv
        xdhA_core_uniref100.fas
        """
        expand_uniref90_to_100(args.gene_prefix, args.core_uniref90, args.output_dir, Logger)

    elif args.mode == "re_cluster":
        """
        Example:
        python ${vritra_dir}recluster_species_UniRef.py ${mode} -g ${gene_name} 
        -o ${parent_folder}/finalized_database 
        -b /gpfs/data/lilab/home/zhoub03/software/blast/ncbi-blast-2.16.0+/bin
        -t ${parent_folder}/dataset_construction_refined/${gene_name}_core_uniref90_100_taxon_ranks.csv 
        -f ${parent_folder}/dataset_construction_refined/${gene_name}_core_uniref100.fas 
        -n ${parent_folder}/recluster_intermediate 
        -L recluster_0526.log
        Output:
        xdhA_core_representative_UniRef100.fas
        xdhA_core_uniref100_correspondence.csv
        xdhA_UniRef100_correspondence.json
        """
        # add bin folder path of ncbi-blast to the linux environment
        if args.blast_bin_dir is not None:
            if args.blast_bin_dir not in os.environ["PATH"]:
                os.environ["PATH"] = os.environ["PATH"] + ":" + args.blast_bin_dir
        uniref_species = UniRefSpecies(args.uniref_taxonomy_path, args.uniref100_fas_path, Logger)
        uniref_species.complete_process(args.gene_prefix, args.intermediate_folder, args.output_dir)

    # Example of how the parameters can be used
    print(f"gene-prefix: {args.gene_prefix}")

