"""
给定一个种子fas文件，从一个较大的seq subset（the complete UniRef90 data set）中 blastp
找到和种子序列similarity相近的所有序列。再给定一个pre-labeled的sequence数据集合，通过我们自己的
算法从上面的集合中找到所有应该被重新标记为目标功能的序列（new-labeled）
下载pre-labeled和new-labeled sequences所对应的UniRef100


module add diamond/0.9.18
module add python/cpu/3.6.5

class GeneBlast: main class Given a gene name and seed sequence, build data set
    ----build_initial_uniref_local: blast seed sequence to get initial seq set
"""

import os

import sys
import pandas as pd
from Bio import SeqIO
import log_setup
import argparse
# import make_UniRefDB as make_uniref_db
# import make_CompleteDB as make_comDB

import construct_database
# import get_seq_by_id_from_uniref as gsbifu
# import detect_community as detect_com
# import download_uniref100_from_uniref90 as download_100_from_90


class GeneBlast:
    def __init__(self, output_prefix, seed_fas_path, seed_type, total_output_dir, logger):
        self.output_prefix = output_prefix      # like gene name, ygeX
        self.seed_path = seed_fas_path
        if seed_type not in ['nucl', 'prot']:
            raise ValueError("Invalid option. seed_type must be 'nucl' or 'prot'.")
        self.seed_type = seed_type
        self.seed_length_dict = construct_database.get_seq_by_id_from_uniref.get_seq_length_from_fas(seed_fas_path, seed_type)    # {"sp|P66899": 398}
        logger.info(f"The length of seed sequences is {self.seed_length_dict}")

        self.total_dir = total_output_dir
        # self.blast_bin_path = blast_bin_path
        # os.environ["PATH"] = os.environ["PATH"] + ":" + self.blast_bin_path
        self.initial_blast_dir = os.path.join(total_output_dir, f"1_blast_initial_set")
        if not os.path.exists(self.initial_blast_dir):
            os.system(f"mkdir -p {self.initial_blast_dir}")
        self.initial_blasted_seq_path = ""              # like PATH/uniref90_blast_by_P66899.fas

        self.initial_blast_res_path = ""
        self.initial_blast_prot_id = []
        self.initial_blast_prot_length_dict = {}        # {"UniRef90_UPI000FDD77E9": 398}
        """
        self.first_query_seq_id = []
        self.first_query_seq_fas_path = ""              # PATH/new_sequences_iteration_0.fas

        self.iteration_final_dir = os.path.join(self.total_dir, f"{self.output_prefix}_blast_iteration_final_set")
        self.target_uniref90_list = []
        """
        self.uniprotkb_incorporated_seq_path = ""   # if labeled seq in UniProtKB is provided
        self.labeled_to_uniref90_blast_res_path = ""
        self.labeled_connected_to_uniref90_blast_res_path = ""
        self.labeled_seq_set = set()
        self.labeled_seq_path = ""          # pre-labeled
        self.labeled_seq_no_uniprotkb = ""  # the self.labeled_seq_path before incorporating uniprotkb
        self.refined_labeled_seq_path = ""  # refined pre-labeled, identity > 35%
        self.refined_labeled_uniref90_set = set()   # uniref90 id set in the "self.refined_labeled_seq_path"
        self.quasi_labeled_seq_path = ""
        self.unlabeled_seq_path = ""    # directly (t1) or indirectly (t2) connected seq, pre-labeled excluded
        self.t1_t2_merged_only_set = set()     # set of IDs in self.unlabeled_seq_path
        self.labeled_and_t1_t2_blast_res_path = ""
        self.pre_labeled_new_labeled_merge_path = ""    # final core UniRef90 set, after length filter
        self.sponge_seqs_path = ""       # rest of t1 (remove new labeled), and all t2 seqs
        self.logger = logger

    def build_initial_uniref_local(self, uniref_path, uniref_name, num_alignments=10000):
        """
        :param uniref_path: default uniref should be UniRef90, like path/uniref90.fasta
        :param uniref_name: should be uniref90 for original design, or uniref50, or uniref100
        :param num_alignments: number of hits for blast
        :return:
        """
        ##################################################################################
        # STEP1: blast the seed fas to UniRef database, to get total potential sequences #
        ##################################################################################
        initial_blast_res_path = os.path.join(self.initial_blast_dir, f"{self.output_prefix}_{uniref_name}_initial_blast.txt")
        blast_method = ""
        if self.seed_type == "nucl":
            blast_method = "blastx"
        elif self.seed_type == "prot":
            blast_method = "blastp"

        if os.path.exists(initial_blast_res_path):
            self.logger.info(f"The results of blast of seed sequences existed in {initial_blast_res_path}. Skip blast of seed sequences!")
        else:
            blast_command = f"{blast_method} -query {self.seed_path} -db {uniref_path} -outfmt 6 " \
                            f"-num_alignments {num_alignments} > {initial_blast_res_path}"
            self.logger.info(f"Blast the seed sequences to {uniref_name} ...\n{blast_command}")
            os.system(blast_command)
        ################################################################################
        # STEP2: build blast database for the initial subset-- "initial_seq_set_path " #
        ################################################################################
        initial_blast_prot_id = construct_database.get_seq_by_id_from_uniref.get_prot_id_from_blast_res(initial_blast_res_path)
        self.initial_blast_res_path = initial_blast_res_path
        self.initial_blast_prot_id = initial_blast_prot_id

        # # output to initial_seq_set_path, get protein length dict
        initial_seq_set = f"{self.output_prefix}_initial_sequence_set.fas"
        initial_seq_set_path = os.path.join(self.initial_blast_dir, initial_seq_set)
        if os.path.exists(initial_seq_set_path):
            self.logger.info(f"The initial-blasted protein sequences existed in {initial_seq_set_path}. Skip extracting the initial-blasted protein sequences!")
        else:
            construct_database.get_seq_by_id_from_uniref.get_seq_from_fas_given_id(uniref_path, initial_seq_set_path, initial_blast_prot_id)
            self.logger.info(f"Extract the initial-blasted protein sequences to {initial_seq_set}")
        self.initial_blast_prot_length_dict = construct_database.get_seq_by_id_from_uniref.get_seq_length_from_fas(initial_seq_set_path, "prot")

        # # make db: the db path is the same with the initial_seq_set_path
        os.chdir(self.initial_blast_dir)
        if os.path.exists(f"{initial_seq_set}.pdb"):
            self.logger.info(f"The initial protein database existed in {initial_seq_set} under {self.initial_blast_dir}. Skip making initial blastdb!")
        else:
            os.system(f"makeblastdb -in {initial_seq_set} -dbtype prot -parse_seqids")
            self.logger.info(f"Made initial blastdb {os.path.join(self.initial_blast_dir, initial_seq_set)}")
        self.initial_blasted_seq_path = os.path.join(self.initial_blast_dir, initial_seq_set)

    # build the refined database by PSI_BLAST (之后需要改成原先使用的iterative search)
    def psi_blast_label_new_seq(self, new_labeled_seq_func_summary_path, new_labeled_seq_fas_path):
        psi_blast_command = ""
        pass

    def incorporate_uniprotkb(self, uniprotkb_fas_path):
        uniprotkb_dir = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs_uniprotkb")
        uniprotkb_intermediate_dir = os.path.join(self.total_dir, "2_dataset_construction_uniprotkb_intermediate")
        if not os.path.exists(uniprotkb_dir):
            os.system(f"mkdir -p {uniprotkb_dir}")
        if not os.path.exists(uniprotkb_intermediate_dir):
            os.system(f"mkdir -p {uniprotkb_intermediate_dir}")
        # 1. get the uniprotkb filtered by length
        uniprotkb_len_filtered_path = os.path.join(uniprotkb_dir, f"{self.output_prefix}_uniprotkb_len_filtered_seqs.fas")
        if not os.path.exists(uniprotkb_len_filtered_path):
            filter_fas_by_len(uniprotkb_fas_path, uniprotkb_len_filtered_path, self.logger, len_upper=0.25, len_lower=0.25)
        # --------------------------------
        # 2. get species refined fas and correspondence csv
        """
        'tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM': SeqRecord(seq=Seq('MEIGRSRHRVDAWSKVTGEAKYTADLFPDNCLTAKVIRSTIANGRVLSMDTREA...
        AYV'), id='tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM', name='tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM', 
        description='tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM Xanthine dehydrogenase molybdenum-binding subunit XdhA 
        OS=Candidatus Onthocola gallistercoris OX=2840876 GN=xdhA PE=4 SV=1', dbxrefs=[])
        """
        refined_species_fas_path = os.path.join(uniprotkb_dir, f"{self.output_prefix}_uniprotkb_species_refined.fas")
        correspondence_csv_path = os.path.join(uniprotkb_dir, f"{self.output_prefix}_uniprotkb_taxon_ranks.csv")
        if not os.path.exists(correspondence_csv_path):
            construct_database.uniprotkb_process.uniprotkb_simplified(uniprotkb_len_filtered_path, refined_species_fas_path,
                                                                  correspondence_csv_path, uniprotkb_intermediate_dir)
        # --------------------------------
        # 3. blast refined_species_fas_path to "self.labeled_seq_path (len filtered)"
        uniprotkb_to_labeled_blast_res_path = os.path.join(uniprotkb_dir, f"{self.output_prefix}_uniprotkb_to_labeled_blasted.txt")
        if os.path.exists(uniprotkb_to_labeled_blast_res_path):
            self.logger.info(
                f"{uniprotkb_to_labeled_blast_res_path} existed. Skip blast UnoProtKB to the labeled sequences.")
        else:
            blast_command = f"blastp -query {refined_species_fas_path} -db {self.labeled_seq_no_uniprotkb} -outfmt 6 " \
                            f"-num_alignments 1 > {uniprotkb_to_labeled_blast_res_path}"
            self.logger.info(f"Blast the filtered labeled sequences to initial uniref90 ...\n{blast_command}")
            os.system(blast_command)
        # --------------------------------
        # 4. get the uniprotkb that not fall into the identity-90 of labeled seqs
        uniprotkb_to_labeled_blast_res = pd.read_table(uniprotkb_to_labeled_blast_res_path, header=None)
        uniprotkb_shared_id = set(uniprotkb_to_labeled_blast_res.loc[uniprotkb_to_labeled_blast_res.iloc[:, 2] >= 90, ].iloc[:, 0])
        uniprotkb_incorporated_fas_path = os.path.join(uniprotkb_dir, f"{self.output_prefix}_uniprotkb_incorporated.fas")
        if not os.path.exists(uniprotkb_incorporated_fas_path):
            output_fas_given_seq(uniprotkb_shared_id, refined_species_fas_path, uniprotkb_incorporated_fas_path, reverse=True)
        self.logger.info(f"Seqs from uniprotkb that should be included are stored in {uniprotkb_incorporated_fas_path}")
        self.uniprotkb_incorporated_seq_path = uniprotkb_incorporated_fas_path

    def filter_labeled_by_length(self, labeled_seq_path, len_upper=0.25, len_lower=0.25):
        self.logger.info(f"Pairwise identity calculation for annotated sequence set to initial uniref90 set ... ...")
        dataset_construction_dir = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs")
        if not os.path.exists(dataset_construction_dir):
            os.system(f"mkdir -p {dataset_construction_dir}")
        self.labeled_seq_no_uniprotkb = os.path.join(dataset_construction_dir,
                                             f"{self.output_prefix}_labeled_filtered_seqs_no_uniprotkb.fas")
        # 1. filter the labeled according to sequence length
        if os.path.exists(self.labeled_seq_no_uniprotkb):
            self.logger.info(f"The filtered labeled sequences {self.labeled_seq_no_uniprotkb} existed. Skip filtering.")
        else:
            filter_fas_by_len(labeled_seq_path, self.labeled_seq_no_uniprotkb, self.logger, len_upper, len_lower)
        # 2. make blast db of length-filtered labeled seqs
        if not os.path.exists(f"{self.labeled_seq_no_uniprotkb}.pdb"):
            os.chdir(dataset_construction_dir)
            os.system(f"makeblastdb -in {self.output_prefix}_labeled_filtered_seqs_no_uniprotkb.fas -dbtype prot -parse_seqids")
            self.logger.info(f"Made blastdb for length-filtered labeled seqs {self.labeled_seq_no_uniprotkb}")

    # blast the labeled sequences to the sequences in UniRef90 initially blasted from seed sequences
    def calculate_labeled_to_uniref90_identity(self, num_alignments=500):
        """
        :param num_alignments:
        :return:
        """
        dataset_construction_dir = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs")
        self.labeled_seq_path = os.path.join(dataset_construction_dir, f"{self.output_prefix}_labeled_filtered_seqs.fas")
        self.labeled_to_uniref90_blast_res_path = os.path.join(dataset_construction_dir,
                                                               f"{self.output_prefix}_labeled_to_uniref90_blasted_res_alignments{num_alignments}.txt")
        # 1. consider whether uniprotKB is provided
        if len(self.uniprotkb_incorporated_seq_path) > 0:
            # merge previous labeled seqs with uniprotKB To "self.labeled_seq_path"
            if not os.path.exists(self.labeled_seq_path):
                os.system(f"cat {self.labeled_seq_no_uniprotkb} {self.uniprotkb_incorporated_seq_path} > {self.labeled_seq_path}")
        else:
            # uniprotKB not provided
            if not os.path.exists(self.labeled_seq_path):
                os.system(f"cp {self.labeled_seq_no_uniprotkb} {self.labeled_seq_path}")

        # 2. blast labeled seqs to initial raw dataset to get the pairwise identity
        if os.path.exists(self.labeled_to_uniref90_blast_res_path):
            self.logger.info(f"{self.labeled_to_uniref90_blast_res_path} existed. Skip blast the labeled sequences to initial uniref90.")
        else:
            blast_command = f"blastp -query {self.labeled_seq_path} -db {self.initial_blasted_seq_path} -outfmt 6 " \
                        f"-num_alignments {num_alignments} > {self.labeled_to_uniref90_blast_res_path}"
            self.logger.info(f"Blast the filtered labeled sequences to initial uniref90 ...\n{blast_command}")
            os.system(blast_command)
            self.logger.info(f"Blast the filtered labeled sequences to initial uniref90 End.")

    # get T-1 and T2, and identity to each other
    def calculate_labeled_connected_to_uniref90_identity(self, quasi_labeled_threshold=80,
                                                         connection_iden_threshold=50, num_alignments=500):
        """
        1. XXX_refined_labeled_filtered_seqs.fas (all qualified labeled), self.refined_labeled_seq_path
        2. XXX_quasi_labeled_filtered_seqs.fas
        3. XXX_labeled_connected_in_uniref90.fas (T1 only, from UniRef90 db)
        4. XXX_labeled_connected_and_its_connected_seq.fas (T1 and T2)
        5. XXX_labeled_search1_connected_and_its_connected_and_labeled_seq.fas (T1 and T2 + labeled seqs)
        """
        dataset_construction_dir = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs")

        labeled_to_uniref90_blast_res = pd.read_table(self.labeled_to_uniref90_blast_res_path, header=None)

        ##############################################################################
        # 1. define labeled seqs; quasi_labeled_seq; labeled_connected_only_seq (t1) #
        ##############################################################################
        # 1.1 all qualified labeled:
        # need to be closer to raw dataset (Identity >35%), this step ensure the labeled seq must be related
        labeled_seq_set = set(labeled_to_uniref90_blast_res.loc[labeled_to_uniref90_blast_res.iloc[:, 2] >= 35, ].iloc[:, 0])
        self.refined_labeled_uniref90_set = set([i for i in labeled_seq_set if i.startswith("UniRef90")])
        self.refined_labeled_seq_path = os.path.join(dataset_construction_dir,
                                                     f"{self.output_prefix}_refined_labeled_filtered_seqs.fas")
        if not os.path.exists(self.refined_labeled_seq_path):
            output_fas_given_seq(labeled_seq_set, self.labeled_seq_path, self.refined_labeled_seq_path)
        # 反过来通过 labeled_seq_set 优化 labeled_to_uniref90_blast_res (先暂时不加)
        # labeled_to_uniref90_blast_res = labeled_to_uniref90_blast_res.loc[labeled_to_uniref90_blast_res.iloc[:, 0] in labeled_seq_set, ]

        # 1.2 quasi: the seqs very close to original labeled seqs (default >80%), but not the same ID as labeled
        # from the target UniRef90 database
        quasi_labeled_seq_set = set(labeled_to_uniref90_blast_res.loc[labeled_to_uniref90_blast_res.iloc[:, 2] >= quasi_labeled_threshold, ].iloc[:, 1])
        quasi_labeled_seq_set = quasi_labeled_seq_set - labeled_seq_set
        self.quasi_labeled_seq_path = os.path.join(dataset_construction_dir,
                                                   f"{self.output_prefix}_quasi_labeled_filtered_seqs.fas")
        if not os.path.exists(self.quasi_labeled_seq_path):
            output_fas_given_seq(quasi_labeled_seq_set, self.initial_blasted_seq_path, self.quasi_labeled_seq_path)

        # 1.3 t1 (labeled_connected_only_seq_set): connected to labeled and remove "labeled" & "quasi-labeled" IDs
        # from the target UniRef90 database
        blast_res_identity_filtered = labeled_to_uniref90_blast_res.loc[labeled_to_uniref90_blast_res.iloc[:, 2] >= connection_iden_threshold]
        labeled_connected_seq_set = set(blast_res_identity_filtered.iloc[:, 1])
        labeled_connected_only_seq_set = labeled_connected_seq_set - labeled_seq_set - quasi_labeled_seq_set
        self.logger.info(f"{len(labeled_seq_set)} seqs in labeled_seq_set; "
                         f"{len(quasi_labeled_seq_set)} seqs in quasi_labeled_seq_set;"
                         f"{len(labeled_connected_only_seq_set)} seqs in labeled_connected_only_seq_set.")
        labeled_only_seq_set = labeled_seq_set - labeled_connected_seq_set
        overlapped_seq_set = labeled_seq_set & labeled_connected_seq_set
        ##########################################################################
        # 2. output sequence unique to the labeled_connected seq set (output T1) #
        ##########################################################################
        labeled_connected_only_seq_path = os.path.join(dataset_construction_dir,
                                                       f"{self.output_prefix}_labeled_connected_in_uniref90.fas")
        if not os.path.exists(labeled_connected_only_seq_path):
            output_fas_given_seq(labeled_connected_only_seq_set, self.initial_blasted_seq_path, labeled_connected_only_seq_path)
        # os.system(f"makeblastdb -in {labeled_connected_only_seq_path} -dbtype prot -parse_seqids")
        # self.logger.info(f"Made blastdb for {labeled_connected_only_seq_path}")

        #############################################################################
        # 3. get the seq identity from label-connected-seq (T1) to initial uniref90 #
        #############################################################################
        self.labeled_connected_to_uniref90_blast_res_path = \
            os.path.join(dataset_construction_dir,
                         f"{self.output_prefix}_labeled_connected_seq_to_uniref90_blasted_res_alignments{num_alignments}.txt")
        blast_command = f"blastp -query {labeled_connected_only_seq_path} " \
                        f"-db {self.initial_blasted_seq_path} " \
                        f"-outfmt 6 " \
                        f"-num_alignments {num_alignments} > " \
                        f"{self.labeled_connected_to_uniref90_blast_res_path}"
        if not os.path.exists(self.labeled_connected_to_uniref90_blast_res_path):
            self.logger.info(f"Blast the sequences connected to labeled sequences to initial uniref90 ...\n{blast_command}")
            os.system(blast_command)
            self.logger.info(f"Blast the sequences connected to labeled sequences to initial uniref90 End.")
        else:
            self.logger.info(f"Blast the sequences connected to labeled sequences, "
                             f"result {self.labeled_connected_to_uniref90_blast_res_path} existed. Skip this step.")
        ############################################
        # 4. get the seq connected to T1 (i.e. T2) #
        ############################################
        # here threshold is the same with "connection_iden_threshold"
        labeled_connected_to_uniref90_blast_res = pd.read_table(self.labeled_connected_to_uniref90_blast_res_path, header=None)
        blast_res_t2_identity_filtered = labeled_connected_to_uniref90_blast_res.loc[
            labeled_connected_to_uniref90_blast_res.iloc[:, 2] >= connection_iden_threshold]
        labeled_connected_seq_set2 = set(blast_res_t2_identity_filtered.iloc[:, 0])
        labeled_connected_t2_seq_set = set(blast_res_t2_identity_filtered.iloc[:, 1])
        t1_t2_merged_set = labeled_connected_seq_set2 | labeled_connected_t2_seq_set
        self.logger.info(f"Get t1 t2 merged set. {len(t1_t2_merged_set)} seqs.")

        # Output T2 and T1, in which "pre_labeled" & "quasi labeled" have been removed
        t1_t2_merged_only_set = t1_t2_merged_set - labeled_seq_set - quasi_labeled_seq_set  # set of UniRef90 IDs
        self.logger.info(f"Remove {len(labeled_seq_set)} labeled_seq_set and "
                         f"{len(quasi_labeled_seq_set)} quasi_labeled_seq_set."
                         f"Get t1 t2 merged only set with {len(t1_t2_merged_set)} seqs.")
        t1_t2_merged_only_seq_path = os.path.join(dataset_construction_dir,
                                                  f"{self.output_prefix}_labeled_connected_and_its_connected_seq.fas")
        if not os.path.exists(t1_t2_merged_only_seq_path):
            output_fas_given_seq(t1_t2_merged_only_set, self.initial_blasted_seq_path, t1_t2_merged_only_seq_path)
        self.logger.info(f"Output t1_t2_merged_only_seq (labeled excluded) to {t1_t2_merged_only_seq_path}")
        self.unlabeled_seq_path = t1_t2_merged_only_seq_path
        self.t1_t2_merged_only_set = t1_t2_merged_only_set

        #####################################################################################
        # 5. merge t1_t2_merged_only with original labeled sequences AND blast to them self #
        #####################################################################################
        # this step is to prepare for "Label Propagation algorithm"
        labeled_t1_t2_merged_seq_path = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs",
                                                  f"{self.output_prefix}_labeled_connected_and_its_connected_and_labeled_seq.fas")
        if not os.path.exists(labeled_t1_t2_merged_seq_path):
            os.system(f"cat {self.refined_labeled_seq_path} {self.quasi_labeled_seq_path} {self.unlabeled_seq_path} > {labeled_t1_t2_merged_seq_path}")
        self.logger.info(f"Merge original labeled and t1_t2_merged_only_seq (labeled excluded) to {labeled_t1_t2_merged_seq_path}")
        """
        # make blast db
        if os.path.exists(f"{labeled_t1_t2_merged_seq_path}.pdb"):
            self.logger.info(f"The initial protein database existed in {labeled_t1_t2_merged_seq_path}. Skip making initial blastdb!")
        else:
            os.system(f"makeblastdb -in {labeled_t1_t2_merged_seq_path} -dbtype prot -parse_seqids")
            self.logger.info(f"Made all merged blastdb {labeled_t1_t2_merged_seq_path}")
        # blast to it self
        labeled_and_t1_t2_blast_res_path = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs",
            f"{self.output_prefix}_labeled_connected_and_its_connected_and_labeled_seq_to_self_blasted_res_alignments{num_alignments*2}.txt")
        self.labeled_and_t1_t2_blast_res_path = labeled_and_t1_t2_blast_res_path
        blast_command = f"blastp -query {labeled_t1_t2_merged_seq_path} -db {labeled_t1_t2_merged_seq_path} -outfmt 6 " \
                        f"-num_alignments {num_alignments*2} > {labeled_and_t1_t2_blast_res_path}"
        if os.path.exists(labeled_and_t1_t2_blast_res_path):
            self.logger.info(f"Blast all merged to itself results existed! Skip this step.")
        else:
            self.logger.info(f"Blast all merged to itself ...\n{blast_command}")
            os.system(blast_command)
            self.logger.info(f"Blast all merged to itself End.")
        """
    def label_new_seq(self, new_labeled_seq_func_summary_path, new_labeled_seq_fas_path, db_name, connection_iden_threshold):
        """
        :param new_labeled_seq_func_summary_path:
        :param new_labeled_seq_fas_path:
        :param db_name: "uniprot" or "interpro"
        :param connection_iden_threshold: should be consistent with the parameter in "calculate_labeled_connected_to_uniref90_identity"
        :return:
        1. XXX_pre_labeled_and_new_labeled_raw.fas
        2. XXX_pre_labeled_and_new_labeled_filtered.fas
        3. XXX_sponge_UniRef90.fas
        """
        dataset_construction_dir = os.path.join(self.total_dir, "2_dataset_construction_by_annotated_seqs")

        # create the fas of all labeled UniRef90s; and the fas of sponge UniRef90s
        if len(self.labeled_connected_to_uniref90_blast_res_path) == 0:
            self.logger.info(f"There is not labeled blasted to initial uniref90 result {self.labeled_connected_to_uniref90_blast_res_path}. "
                             f"Need to run labeled to initial uniref90 first.")
        else:
            # closer_labeled_dict is like {seq_id: {"to_labeled": {"seq_id": "", "highest_identity": -1}, "to_unlabeled": {"seq_id": "", "highest_identity": -1}}}
            final_summary_table, closer_labeled_dict, closer_unlabeled_dict = construct_database.detect_community.detect_new_seq_in_community(
                self.refined_labeled_seq_path, self.quasi_labeled_seq_path, self.unlabeled_seq_path,
                self.initial_blasted_seq_path, self.labeled_to_uniref90_blast_res_path,
                self.labeled_connected_to_uniref90_blast_res_path)

            # Output function summary of pre-labeled and new-labeled sequences
            final_summary_table.to_csv(new_labeled_seq_func_summary_path, sep='\t', index=False)
            # Output new labeled seq to fas
            output_fas_given_seq(set(closer_labeled_dict.keys()), self.unlabeled_seq_path, new_labeled_seq_fas_path)

            # merge pre-labeled with new-labeled sequences (before filter)
            pre_labeled_new_labeled_merge_path = os.path.join(dataset_construction_dir,
                                                              f"{self.output_prefix}_pre_labeled_and_new_labeled_raw.fas")
            if db_name == "uniprot":
                """
                refined_labeled_seq_uniref90_path = os.path.join(dataset_construction_dir,
                                                                 f"{self.output_prefix}_refined_labeled_filtered_uniref90_seqs.fas")
                if not os.path.exists(refined_labeled_seq_uniref90_path):
                    output_fas_given_seq(self.refined_labeled_uniref90_set, self.refined_labeled_seq_path,
                                         refined_labeled_seq_uniref90_path)
                if len(self.uniprotkb_incorporated_seq_path) > 0 and os.path.exists(self.uniprotkb_incorporated_seq_path):
                    # if use UniProtKB; from "self.refined_labeled_seq_path" extract the seq from uniref90 only
                    os.system(f"cat {refined_labeled_seq_uniref90_path} {self.quasi_labeled_seq_path} {new_labeled_seq_fas_path} > {pre_labeled_new_labeled_merge_path}")
                else:
                # only UniRef90 without UniProtKB
                """

                os.system(f"cat {self.refined_labeled_seq_path} {self.quasi_labeled_seq_path} {new_labeled_seq_fas_path} > {pre_labeled_new_labeled_merge_path}")
            else:
                os.system(f"cat {self.quasi_labeled_seq_path} {new_labeled_seq_fas_path} > {pre_labeled_new_labeled_merge_path}")

            # filter the merged file by sequence length (refined UniRef90 dataset without sponge sequence)
            pre_labeled_new_labeled_merge_filtered_path = os.path.join(dataset_construction_dir,
                                                              f"{self.output_prefix}_pre_labeled_and_new_labeled_filtered.fas")
            filter_fas_by_len(pre_labeled_new_labeled_merge_path, pre_labeled_new_labeled_merge_filtered_path,
                              self.logger, len_upper=0.25, len_lower=0.25, by_sd=True)
            self.pre_labeled_new_labeled_merge_path = pre_labeled_new_labeled_merge_filtered_path

            # get pure sponge seqs (seqs in t1 closer to unlabeled, and all seqs in t2),
            # remove "closer to labeled" from t1_t2_merged_only_set
            sponge_id_set = self.t1_t2_merged_only_set - set(closer_labeled_dict.keys())
            sponge_seqs_path = os.path.join(dataset_construction_dir,
                                            f"{self.output_prefix}_sponge_UniRef90.fas")
            self.logger.info(f"Extract the rest {len(sponge_id_set)} sequences as sponge sequences to {sponge_seqs_path}.")
            output_fas_given_seq(sponge_id_set, self.unlabeled_seq_path, sponge_seqs_path)
            self.sponge_seqs_path = sponge_seqs_path


def output_fas_given_seq(seq_id_set, original_fas_path, extracted_fas_path, reverse=False):
    """
    Given the target seq id set, extract seqs to fas file from a fas file (ONLY FOR SMALL SAMPLE SIZE)
    :param seq_id_set:
    :param original_fas_path:
    :param extracted_fas_path:
    :param reverse: if True, output the seq not in the seq_id_set
    :return:
    """
    original_seq_record_dict = SeqIO.to_dict(SeqIO.parse(original_fas_path, "fasta"))
    original_seq_set = set(original_seq_record_dict.keys())
    extracted_fas_file = open(extracted_fas_path, "w")
    if reverse:
        # output the id not in given seq_id_set
        extracted_seq_id = original_seq_set - seq_id_set
        for seq_id_temp in extracted_seq_id:
            SeqIO.write(original_seq_record_dict[seq_id_temp], extracted_fas_file, "fasta")
        print(f"From {len(original_seq_set)} original seqs, extract {len(extracted_seq_id)} seqs not in {len(seq_id_set)} given ids.")
    else:
        # output the id in given seq_id_set
        extracted_seq_id = original_seq_set & seq_id_set
        for seq_id_temp in extracted_seq_id:
            SeqIO.write(original_seq_record_dict[seq_id_temp], extracted_fas_file, "fasta")
        print(f"From {len(original_seq_set)} original seqs, extract {len(extracted_seq_id)} seqs in {len(seq_id_set)} given ids.")
    extracted_fas_file.close()


def filter_fas_by_len(original_fas_path, extracted_fas_path, logger, len_upper=0.25, len_lower=0.25, by_sd=False):
    """"""
    original_seq_record_dict = SeqIO.to_dict(SeqIO.parse(original_fas_path, "fasta"))
    seq_len_dict = {}       # {seq_id: seq_len}
    for seq_id, seq_record in original_seq_record_dict.items():
        seq_len_dict.update({seq_id: len(seq_record.seq)})
    seq_len_series = pd.Series(list(seq_len_dict.values())) # convert the protein lengths to pandas series
    if by_sd:
        # filter the seq length by mean +- 2sd
        seq_len_mean = seq_len_series.mean()
        seq_len_std = seq_len_series.std()
        seq_len_upper = seq_len_mean + 2 * seq_len_std
        seq_len_lower = seq_len_mean - 2 * seq_len_std
        logger.info(f"The mean length of sequence is {seq_len_mean}, upper threshold is {seq_len_upper}, lower threshold is {seq_len_lower}.")
    else:
        # get the median length of sequences and upper/lower threshold
        seq_len_median = seq_len_series.median()
        seq_len_upper = seq_len_median * (1 + len_upper)
        seq_len_lower = seq_len_median * (1 - len_lower)
        logger.info(f"The median length of sequence is {seq_len_median}, upper threshold is {seq_len_upper}, lower threshold is {seq_len_lower}.")
    # output the filtered sequence with length fall into the thresholds
    extracted_fas_file = open(extracted_fas_path, "w")
    for seq_id, seq_record in original_seq_record_dict.items():
        if (seq_len_dict[seq_id] >= seq_len_lower) and (seq_len_dict[seq_id] <= seq_len_upper):
            SeqIO.write(seq_record, extracted_fas_file, "fasta")
    extracted_fas_file.close()
    logger.info(f"The filtered sequences are stored in {extracted_fas_path}")


############################
# Main function of model 1 #
############################
def construct_by_annotated_seqs(gene_prefix, seed_fas_path, seed_type, gene_total_output_dir, labeled_seq_path,
                                # labeled_prefix,
                                labeled_seq_uniprotkb_path,
                                db_name, quasi_thresh, search_thresh, uniref_path, num_align, logger_path):
    logger = log_setup.setup_log_file(logger_path)
    # 1. initialize the class "GeneBlast"
    gene_blast = GeneBlast(gene_prefix, seed_fas_path, seed_type, gene_total_output_dir, logger)

    # 2. blast seed to complete uniref90 to get initial set
    # # Uniref_path, Uniref_name = ["/gpfs/data/lilab/home/zhoub03/software/UniRef/uniref90.fasta", "uniref90"]
    gene_blast.build_initial_uniref_local(uniref_path, "uniref90", num_alignments=num_align)

    # 3.1 filter labeled seq by length (get self.labeled_seq_path)
    gene_blast.filter_labeled_by_length(labeled_seq_path)

    # 3. if "labeled_seq_uniprotkb_path" is not "None", refine and reduce the seqs from UniProtKB first
    if labeled_seq_uniprotkb_path:
        # get uniprotkb that should be incorporated (self.uniprotkb_incorporated_seq_path)
        gene_blast.incorporate_uniprotkb(labeled_seq_uniprotkb_path)

    # 3. blast labeled seq (has consider both uniref and uniprotkb) to initial set
    gene_blast.calculate_labeled_to_uniref90_identity(num_alignments=500)

    # 4. get seq connected to labeled seq (t1) and seq connected to t1 (t2); as a whole, get pairwise identity
    gene_blast.calculate_labeled_connected_to_uniref90_identity(quasi_thresh, search_thresh, num_alignments=500)

    # 5. get new labeled seq set (including pre labeled) and summary of sequence number by functions
    New_labeled_seq_func_summary_path = os.path.join(gene_total_output_dir, "2_dataset_construction_by_annotated_seqs", f"{gene_prefix}_pre_labeled_and_new_labeled_func_summary.txt")
    New_labeled_seq_fas_path = os.path.join(gene_total_output_dir, "2_dataset_construction_by_annotated_seqs", f"{gene_prefix}_new_labeled_seq.fas")
    # create the fas of all labeled UniRef90s; and the fas of sponge UniRef90s
    gene_blast.label_new_seq(New_labeled_seq_func_summary_path, New_labeled_seq_fas_path, db_name, search_thresh)



####################################################################
# 下面这个暂时还没改可能最后不用PSI-blast,改回先前的iterative search #
####################################################################
def construct_by_psi_blast(gene_prefix, seed_fas_path, seed_type, gene_total_output_dir, labeled_seq_path, labeled_prefix, uniref_path, num_align, logger_path):
    logger = log_setup.setup_log_file(logger_path)
    # 1. initialize
    # blast_bin_path = "/gpfs/data/lilab/home/zhoub03/software/blast/ncbi-blast-2.16.0+/bin"
    gene_blast = GeneBlast(gene_prefix, seed_fas_path, seed_type, gene_total_output_dir, logger)

    # 2. blast seed to complete uniref90 to get initial set
    # # Uniref_path, Uniref_name = ["/gpfs/data/lilab/home/zhoub03/software/UniRef/uniref90.fasta", "uniref90"]
    gene_blast.build_initial_uniref_local(uniref_path, "uniref90", num_alignments=num_align)

    # 3. blast labeled seq to initial set
    gene_blast.calculate_labeled_to_uniref90_identity(labeled_seq_path, labeled_prefix, num_alignments=500)

    # 4. get seq connected to labeled seq (t1) and seq connected to t1 (t1); as a whole, get pairwise identity
    gene_blast.calculate_labeled_connected_to_uniref90_identity(labeled_prefix,
                                                                connection_iden_threshold=50, num_alignments=500)
    # 5. get new labeled seq set (including pre labeled) and summary of sequence number by functions
    # New_labeled_seq_func_summary_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/P19409_baiB_202502/dataset_construction/pre_labeled_new_labeled_func_summary.txt"
    # New_labeled_seq_fas_path = "/gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/P19409_baiB_202502/dataset_construction/baiB_new_labeled_seq.fas"
    New_labeled_seq_func_summary_path = os.path.join(gene_total_output_dir, "2_dataset_construction_by_annotated_seqs", f"{gene_prefix}_{labeled_prefix}_pre_labeled_and_new_labeled_func_summary.txt")
    New_labeled_seq_fas_path = os.path.join(gene_total_output_dir, "2_dataset_construction_by_annotated_seqs", f"{gene_prefix}_{labeled_prefix}_new_labeled_seq.fas")
    # gene_blast.label_new_seq(New_labeled_seq_func_summary_path, New_labeled_seq_fas_path, db_name, )

    # 6. Download the corresponding UniRef100 of UniRef90; & generate correspondence between UniRef100 and taxonomy
    # 这部分挪到了后面recluster的程序中运行


def get_args():
    """Parses command-line arguments and returns them."""

    # Define parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-g", "--gene-prefix", type=str, required=True, help="Prefix for gene files")
    parent_parser.add_argument("-s", "--seed-fas-path", type=str, required=True, help="Path to the seed FASTA file")
    parent_parser.add_argument("-t", "--seed-type", type=str, required=True, choices=['prot', 'nucl'], default='prot',
                               help="Type of the seed")
    parent_parser.add_argument("-o", "--gene-total-output-dir", type=str, required=True,
                               help="Directory of total output for this gene")
    parent_parser.add_argument("-b", "--blast-bin-dir", type=str, required=True,
                               help="The directory path of bin of downloaded ncbi-blast. Specify this parameter if "
                                     "ncbi-blast has not been added to the environment.")
    parent_parser.add_argument("-L", "--logger-path", type=str, required=True, help="Path to save log files")

    # Main parser
    parser = argparse.ArgumentParser(description="Construct sequence database from seed seq using Annotated Seqs or PSI-BLAST.")
    subparsers = parser.add_subparsers(dest="mode", help="Mode of operation")

    # Subparser for 'pre_construct' mode (shared + specific args)

    "--download-uniref90"

    # Subparser for 'annotated_seqs' mode (shared + specific args)
    parser_annotated = subparsers.add_parser("annotated_seqs", parents=[parent_parser],
                                         help="Run in annotated sequence mode")
    parser_annotated.add_argument("-l", "--annotated-seq-path", type=str, required=True, help="Path to annotated sequence file")
    # parser_annotated.add_argument("-p", "--annotated-prefix", type=str, required=True, help="Prefix for annotated sequences")
    parser_annotated.add_argument("-d", "--annotated-database", type=str, required=True, choices=['uniprot', 'interpro'],
                               help="The name of the database where the annotated seqs are from.")
    parser_annotated.add_argument("-u", "--uniref90-db", type=str, required=True, help="Path to the UniRef90 database")
    parser_annotated.add_argument("--annotated-uniprotkb-path", type=str,
                                  help="Path to annotated sequence file from UniProtKB (Optional).")
    parser_annotated.add_argument("-q", "--quasi-threshold", type=float, default=80, help='Threshold value of quasi labeled seqs (default: 80)')
    parser_annotated.add_argument("-r", "--search-threshold", type=float, default=50, help="Threshold value of searching range (default: 50)")
    parser_annotated.add_argument("--seed-align-number", type=int, default=15000, help="Number of alignments by seed seq (default: 15000)")

    # Subparser for 'psi_blast' mode (only shared args)
    parser_psi = subparsers.add_parser("psi_blast", parents=[parent_parser], help="Run in PSI-BLAST mode")
    parser_psi.add_argument("-e", "--e_value", type=float, default=1e-3, help="E-value threshold for PSI-BLAST search (default: 1e-3)")
    args1 = parser.parse_args()
    if args1.mode is None:
        parser.error("A subcommand is required. Use -h for help.")
    return args1


if __name__ == "__main__":
    args = get_args()
    if args.mode == "annotated_seqs":
        """
        python construct_database.py annotated_seqs -g ygeX -s /gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/P66899_ygeX_202502/P66899_ygeX.fas 
        -t prot -o /gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/P66899_ygeX_202502 -b /gpfs/data/lilab/home/zhoub03/software/blast/ncbi-blast-2.16.0+/bin
        -l /gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/P66899_ygeX_202502/uniref_Diaminopropionate_ammonia_lya_2025_02_19.fasta 
        -p search_C -L /gpfs/data/lilab/home/zhoub03/generalized_pipeline_20240925/result/P66899_ygeX_202502/YGEX_search_C_contruct_database_0219.log
        """
        if args.blast_bin_dir not in os.environ["PATH"]:
            # blast_bin_path = "/gpfs/data/lilab/home/zhoub03/software/blast/ncbi-blast-2.16.0+/bin"
            os.environ["PATH"] = os.environ["PATH"] + ":" + args.blast_bin_dir
        construct_by_annotated_seqs(args.gene_prefix, args.seed_fas_path, args.seed_type, args.gene_total_output_dir,
                                    args.annotated_seq_path,
                                    args.annotated_uniprotkb_path, args.annotated_database,
                                    args.quasi_threshold, args.search_threshold, args.uniref90_db,
                                    args.seed_align_number, args.logger_path)
    elif args.mode == "psi_blast":
        pass

    # Example of how the parameters can be used
    print(f"Gene Prefix: {args.gene_prefix}")
    print(f"Seed FASTA Path: {args.seed_fas_path}")
    print(f"Seed Type: {args.seed_type}")
    print(f"Gene Total Output Directory: {args.gene_total_output_dir}")
    # print(f"Labeled Sequence Path: {args.labeled_seq_path}")
    # print(f"Labeled Prefix: {args.labeled_prefix}")
    print(f"Logger Path: {args.logger_path}")
    print("End of construct_procedure.")


