"""
module add python/cpu/3.7.2

Script Name: diamond_summarize.py
Author: Boyan Zhou
Date: 2025-05-29
Description: []

Function Overview:
    1. ():
        - Purpose:
        - Dependencies:
"""

import subprocess
import pandas as pd
import os


##############
# Main Class #
##############
def summarize_diamond_report(gene_prefix, read_count_filename, sample_fq_table_path, sample_report_table_path,
                             taxonomy_mapping_csv, output_dir, my_logger, iden_thresh=70, matched_len_thresh=25):
    """

    :param gene_prefix: gene symbol, e.g. FRC
    :param read_count_filename: precision_medicine_MGX_read_count.tsv (output read count file)
    :param sample_fq_table_path: mapping file of sample ID and fqs
    :param sample_report_table_path: mapping file of sample ID and diamond report
    :param taxonomy_mapping_csv: mapping file of uniref seq and taxonomy
    :param output_dir:
    :param my_logger:
    :param iden_thresh: identity threshold of aligned reads
    :param matched_len_thresh: length threshold of aligned reads
    :return:
    """
    ##################################################
    # get the sample with both fq and diamond report #
    ##################################################
    sample_fq_dict = sample_dict_from_csv(sample_fq_table_path, "fq", my_logger)
    sample_diamond_dict = sample_dict_from_csv(sample_report_table_path, "diamond_report", my_logger)
    sample_valid_list = list(sample_fq_dict.keys() & sample_diamond_dict.keys())
    my_logger.info(f"There are {len(sample_fq_dict)} fastq samples and {len(sample_diamond_dict)} samples.\n"
                   f"The number of overlapping samples are {len(sample_valid_list)}: {sample_valid_list}")
    # only include sample with both fq and diamond output
    sample_info_dict = {i: {"diamond_result_path_list": sample_diamond_dict[i], "fq_path_list": sample_fq_dict[i]}
                        for i in sample_valid_list}

    #############
    # summarize #
    #############
    diamond_sum = DiamondSum(gene_prefix, taxonomy_mapping_csv, sample_info_dict, output_dir, my_logger)
    read_count_report_path = read_count_filename
    if os.path.dirname(read_count_filename) == '':
        read_count_report_path = os.path.join(output_dir, read_count_filename)
    diamond_sum.get_reads_count(read_count_report_path)
    diamond_sum.diamond_result_summary(read_count_report_path, iden_thresh, matched_len_thresh)


def sample_dict_from_csv(sample_info_path, sample_type, my_logger):
    """
    :param sample_info_path:
    :param sample_type="fq" or "diamond_report"
    :param my_logger:
    :return:
    """
    if sample_type not in ("fq", "diamond_report"):
        raise ValueError(f"Current sample_type is {sample_type}, it must be 'fq' or 'diamond_report'")

    # 1. get sample fq-path dict
    sample_file_dict = {}     # {sample: ["PATH/XXX_R1.fq", "PATH/XXX_R2.fq"]}
    with open(sample_info_path, "r") as sample_info_f:
        for line in sample_info_f:
            cols = line.strip().split(",")
            fq_path_list = cols[1].split(":")   # ["PATH/XXX.fq"], ["PATH/XXX_R1.fq", "PATH/XXX_R2.fq"]
            if len(fq_path_list) > 2:
                print(f"Warning! The sample {cols[0]} has more than two {sample_type} samples {fq_path_list}. Skip this sample!")
                my_logger.info(f"Warning! The sample {cols[0]} has more than two fq samples {fq_path_list}. Skip this sample!")
            else:
                sample_file_dict.update({cols[0]: fq_path_list})
    return sample_file_dict


##############
# Main Class #
##############
class DiamondSum:
    def __init__(self, gene_name, seq_taxonomy_mapping_csv_path, sample_info_dict, output_folder, logger):
        """
        :param sample_info_dict: {sample_ID: {"fq_path_list":[R1.fq, ]}} or {sample_ID: {"fq_path_list":[path.fq]}}
        """
        self.seq_taxonomy_mapping_csv_path = seq_taxonomy_mapping_csv_path
        # initialize from a given core set of uniref100
        self.logger = logger
        self.gene_name = gene_name
        # self.diamond_db_path = core_uniref100_dmnd_path

        self.sample_info_dict = sample_info_dict        # {sample_ID: {"fq_path_list":[,]}}
        self.sample_uniref_count_dict = {}              # {sample_ID: uniref_count_table} summarized from diamond result
        self.output_folder = output_folder

    def get_reads_count(self, read_count_report_path):
        self.logger.info(f"Counting the number of reads for each sample using __get_reads_count ... ... ")
        if os.path.exists(read_count_report_path):
            print(f"The read count file already exists. Skip this step! {read_count_report_path}")
            self.logger.info(f"The read count file already exists. Skip this step! {read_count_report_path}")
            return 1
        read_count_report = open(read_count_report_path, "w")
        for sample_id, sample_info in self.sample_info_dict.items():
            if os.path.exists(sample_info["fq_path_list"][0]):
                read_count = count_reads_shell(sample_info["fq_path_list"][0])
                self.logger.info(f"The reads count of {sample_info['fq_path_list'][0]} is {read_count}")
                sample_info["reads_count"] = read_count
                # NOTICE: here read_count is only R1 for paired-end fq sample
                read_count_report.write(f"{sample_id}\t{sample_info['fq_path_list'][0]}\t{read_count}\n")
            else:
                self.logger.info(f"No file: {sample_info['fq_path_list'][0]} for sample {sample_id}!")
        read_count_report.close()

    def diamond_result_summary(self, read_count_report_path, identical_percent_threshold=70, matched_length_threshold=25):
        # 0. Read the correspondence table
        self.logger.info(f"In diamond_result_summary, step0: Read the correspondence table")
        # seq_taxonomy_df = pd.read_csv("C:\\Users\\boyan\\OneDrive - NYU Langone Health\\research\\uric_acid_to_xanthine_20240528\\generalized_pipeline_20240715\\code\\Vritra\\Vritra_v2_github_with_annotation_20260319\\summarize_contribution\\FRC_finalized_correspondence.csv")
        seq_taxonomy_df = pd.read_csv(self.seq_taxonomy_mapping_csv_path)


        # 1. Get the reads count for each fq/fq.gz sample (have considered paired end as one sample)
        self.logger.info(f"In diamond_result_summary, step1: Get the reads count for each fq/fq.gz sample")
        if not os.path.exists(read_count_report_path):
            self.get_reads_count(read_count_report_path)    # self.sample_info_dict[sample_id]["reads_count"]
        read_count_report_df = pd.read_csv(read_count_report_path, sep="\t", header=None)
        sample_read_count_dict = dict(zip(read_count_report_df.iloc[:, 0].astype(str), read_count_report_df.loc[:, 2].astype(int)))

        # 2. Collect the count of each UniRef100 for each sample-id
        self.logger.info(f"In diamond_result_summary, step2: Collect the count of each UniRef100 for each sample-id")
        for sample_id, sample_info in self.sample_info_dict.items():
            self.logger.info(f"Collect the sample {sample_id} {sample_info}")
            if 'diamond_result_path_list' not in sample_info:
                # 2.1 Check whether the diamond results exist
                continue
            else:
                # 2.2 Check Single-End or Paired-End
                if len(sample_info['diamond_result_path_list']) == 1:
                    ###################################
                    # if single end sequencing sample #
                    ###################################
                    diamond_result_path = sample_info['diamond_result_path_list'][0]
                    uniref_count_dict = blastx_to_taxon_count_unipro(diamond_result_path, identical_percent_threshold, matched_length_threshold)
                else:
                    ###################################################################
                    # if pair end sequencing sample, count in R1/R2 is counted as 0.5 #
                    ###################################################################
                    diamond_result_path1 = sample_info['diamond_result_path_list'][0]
                    diamond_result_path2 = sample_info['diamond_result_path_list'][1]
                    uniref_count_dict1 = blastx_to_taxon_count_unipro(diamond_result_path1, identical_percent_threshold, matched_length_threshold)
                    uniref_count_dict2 = blastx_to_taxon_count_unipro(diamond_result_path2, identical_percent_threshold, matched_length_threshold)
                    uniref_count_dict = {key: (uniref_count_dict1.get(key, 0) + uniref_count_dict2.get(key, 0))/2 for key in set(uniref_count_dict1) | set(uniref_count_dict2)}
                # sample_info['uniref_count_dict'] = uniref_count_dict
                self.sample_uniref_count_dict.update({sample_id: uniref_count_dict})

        # 3. Merge the count info to a whole data frame
        self.logger.info(f"In diamond_result_summary, step3: Merge the count info to a whole data frame")
        # # convert each dictionary into a DataFrame and merge them, row is uniref100, col is sample_id
        uniref_count_df = pd.DataFrame.from_dict(self.sample_uniref_count_dict, orient='columns').fillna(0)
        # # convert NaN to integer (optional, if all values are integers)
        uniref_count_df = uniref_count_df.astype(int)

        # 4. Calculate CPM and RPKM
        self.logger.info(f"In diamond_result_summary, step4: Calculate CPM and RPKM")
        # # reorder the sample-reads_count series, to align with uniref_count_df (col order is sample)
        total_reads_series = pd.Series(sample_read_count_dict).reindex(uniref_count_df.columns)
        # # CPM = (reads_N * 1e6)/total_N; RPKM = (reads_N * 1e9)/(total_N * gene_length)
        cpm_df = uniref_count_df.div(total_reads_series, axis=1) * 1e6
        uniref_length_dict = dict(zip(seq_taxonomy_df['UniRef100_ID'].astype(str), seq_taxonomy_df['Protein_length'].astype(int)))
        # # reorder the gene_length series, to align with uniref_count_df (row order is seqs)
        gene_length_series = pd.Series(uniref_length_dict).reindex(uniref_count_df.index)
        self.logger.info(f"uniref_count_df is {uniref_count_df}\ntotal_reads_series is {total_reads_series}\ngene_length_series is {gene_length_series}")
        rpkm_df = (uniref_count_df * 1e9).div(total_reads_series, axis=1).div(gene_length_series, axis=0)
        # # 4.1 Merge taxonomy info with count, cpm and rpkm
        uniref_count_with_info_df = seq_taxonomy_df.merge(uniref_count_df, left_on="UniRef100_ID", right_index=True, how="inner")
        cpm_with_info_df = seq_taxonomy_df.merge(cpm_df, left_on="UniRef100_ID", right_index=True, how="inner")
        rpkm_with_info_df = seq_taxonomy_df.merge(rpkm_df, left_on="UniRef100_ID", right_index=True, how="inner")
        # # output
        uniref_count_with_info_df.to_csv(os.path.join(self.output_folder, f"{self.gene_name}_uniref_count_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=False)
        cpm_with_info_df.to_csv(os.path.join(self.output_folder, f"{self.gene_name}_uniref_cpm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=False)
        rpkm_with_info_df.to_csv(os.path.join(self.output_folder, f"{self.gene_name}_uniref_rpkm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=False)
        # --------------------------------------------------------------------------------------------------------------
        # # 4.2 Merge the row of count, cpm and rpkm df with the same rep-species
        uniref_rep_species_dict = dict(zip(seq_taxonomy_df['UniRef100_ID'].astype(str), seq_taxonomy_df['rep_species'].astype(str)))
        uniref_count_df.index = uniref_count_df.index.map(uniref_rep_species_dict)
        uniref_count_df_species_level = uniref_count_df.groupby(uniref_count_df.index).sum()
        cpm_df.index = cpm_df.index.map(uniref_rep_species_dict)
        cpm_df_species_level = cpm_df.groupby(cpm_df.index).sum()
        rpkm_df.index = rpkm_df.index.map(uniref_rep_species_dict)
        rpkm_df_species_level = rpkm_df.groupby(rpkm_df.index).sum()
        # # output
        uniref_count_df_species_level.to_csv(os.path.join(self.output_folder, f"{self.gene_name}_species_count_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=True)
        cpm_df_species_level.to_csv(os.path.join(self.output_folder, f"{self.gene_name}_species_cpm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=True)
        rpkm_df_species_level.to_csv(os.path.join(self.output_folder, f"{self.gene_name}_species_rpkm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=True)


def count_reads_shell(fastq_path):
    # fastq_path = "path/example_R1.fq.gz"
    if fastq_path.endswith(".gz"):
        cmd = f"zcat {fastq_path} | wc -l"
    else:
        cmd = f"cat {fastq_path} | wc -l"
    # result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    line_count = int(result.stdout.strip())
    read_count = int(line_count/4)  # Each read has 4 lines
    return read_count


# (0712 updated)
def blastx_to_taxon_count_unipro(diamond_report_path, identical_percent_threshold=70, matched_length_threshold=25):
    """
    This is for the results of diamond using UniPro database (not InterPro)
    For each unique read, only count the first record
    :param diamond_report_path:
    :param identical_percent_threshold: threshold of identical aa percent, default 70
    :param matched_length_threshold: the aa length >= 25 (that is >=75 in bp)
    :return: reads count dict of each UniRef, {"UniRef100_A0A2D6BS33": 10, "UniRef100_UPI000C7BF04C": 16}
    """
    protein_count_dict = {}     # {"UniRef100_": 35, ...}
    diamond_raw_list = []
    with open(diamond_report_path, "r") as report_f:
        for line in report_f:
            if line.startswith("#"):
                continue
            else:
                diamond_raw_list.append(line.strip().split("\t"))
    diamond_raw_df = pd.DataFrame(diamond_raw_list)

    if diamond_raw_df.empty:
        print(f"No valid DIAMOND hits found (only comments or empty file) for {diamond_report_path}")
        return protein_count_dict

    # flexible assignment of colnames according to the real output column number
    diamond_raw_df.columns = ["Query ID", "Subject ID", "Percentage of identical matches", "Alignment length",
                              "Number of mismatches", "Number of gap openings", "Start of alignment in query",
                              "End of alignment in query", "Start of alignment in subject",
                              "End of alignment in subject", "Expected value", "Bit score", "Query length"][0:diamond_raw_df.shape[1]]

    # diamond_raw_df["Percentage of identical matches"][0]
    diamond_raw_df.iloc[:, 2:] = diamond_raw_df.iloc[:, 2:].apply(pd.to_numeric)

    diamond_best_match_df = diamond_raw_df[
        ~diamond_raw_df["Subject ID"].str.startswith("UniRef90", na=False)]

    # 1. Only the best match for each reads
    diamond_best_match_df = diamond_best_match_df.drop_duplicates(subset='Query ID', keep='first')

    # 2. Filter by alignment length, And output
    diamond_best_match_filter1_df = diamond_best_match_df[diamond_best_match_df["Alignment length"] >= matched_length_threshold]
    diamond_best_match_filter2_df = diamond_best_match_filter1_df[diamond_best_match_filter1_df["Percentage of identical matches"] >= identical_percent_threshold]
    # diamond_best_match_length_df.to_csv(filtered_diamond_report_path, sep='\t', index=False)

    # count of each uniref100
    protein_count_dict = diamond_best_match_filter2_df['Subject ID'].value_counts().to_dict()

    return protein_count_dict

