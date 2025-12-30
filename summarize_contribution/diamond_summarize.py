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


class DiamondSum:
    def __init__(self, gene_name, core_uniref100_dmnd_path, seq_taxonomy_mapping_csv_path, sample_info_dict, output_folder, logger):
        """
        :param sample_info_dict: {sample_ID: {"fq_path_list":[R1.fq, ]}} or {sample_ID: {"fq_path_list":[path.fq]}}
        """
        self.seq_taxonomy_mapping_csv_path = seq_taxonomy_mapping_csv_path
        # initialize from a given core set of uniref100
        self.logger = logger
        self.gene_name = gene_name
        self.diamond_db_path = core_uniref100_dmnd_path

        self.sample_info_dict = sample_info_dict        # {sample_ID: {"fq_path_list":[,]}}
        self.sample_uniref_count_dict = {}              # {sample_ID: uniref_count_table} summarized from diamond result
        self.output_folder = output_folder

    def get_reads_count(self, read_count_report_path):
        self.logger.info(f"Counting the number of reads for each sample using __get_reads_count ... ... ")
        read_count_report = open(read_count_report_path, "w")
        for sample_id, sample_info in self.sample_info_dict.items():
            if os.path.exists(sample_info["diamond_result_path_list"][0]):
                read_count = count_reads_shell(sample_info["diamond_result_path_list"][0])
                sample_info["reads_count"] = read_count
                # NOTICE: here read_count is only R1 for paired-end fq sample
                read_count_report.write(f"{sample_id}\t{sample_info['diamond_result_path_list'][0]}\t{read_count}\n")
            else:
                self.logger.info(f"No file: {sample_info['diamond_result_path_list'][0]} for sample {sample_id}!")
        read_count_report.close()

    def diamond_result_summary(self, read_count_report_path, output_folder, identical_percent_threshold=70, matched_length_threshold=25):
        # 0. Read the correspondence table
        self.logger.info(f"In diamond_result_summary, step0: Read the correspondence table")
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
            if 'diamond_result_path_list' not in sample_info:
                # 2.1 Check whether the diamond results exist
                continue
            else:
                # 2.2 Check Single-End or Paired-End
                if len(sample_info['diamond_result_path_list']) == 1:
                    # if single end sequencing sample
                    diamond_result_path = sample_info['diamond_result_path_list'][0]
                    uniref_count_dict = blastx_to_taxon_count_unipro(diamond_result_path, identical_percent_threshold, matched_length_threshold)
                else:
                    # if pair end sequencing sample, count in R1/R2 is counted as 0.5
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
        # # reorder the sample-reads_count series, to align with uniref_count_df (col order)
        total_reads_series = pd.Series(sample_read_count_dict).reindex(uniref_count_df.columns)
        # # CPM = (reads_N * 1e6)/total_N; RPKM = (reads_N * 1e9)/(total_N * gene_length)
        cpm_df = uniref_count_df.div(total_reads_series, axis=1) * 1e6
        uniref_length_dict = dict(zip(seq_taxonomy_df['UniRef100_ID'].astype(str), seq_taxonomy_df['Protein_length'].astype(int)))
        gene_length_series = pd.Series(uniref_length_dict).reindex(uniref_count_df.columns)
        rpkm_df = (uniref_count_df * 1e9).div(total_reads_series, axis=1).div(gene_length_series, axis=0)
        # # 4.1 Merge taxonomy info with count, cpm and rpkm
        uniref_count_with_info_df = seq_taxonomy_df.merge(uniref_count_df, left_on="UniRef100_ID", right_index=True, how="inner")
        cpm_with_info_df = seq_taxonomy_df.merge(cpm_df, left_on="UniRef100_ID", right_index=True, how="inner")
        rpkm_with_info_df = seq_taxonomy_df.merge(rpkm_df, left_on="UniRef100_ID", right_index=True, how="inner")
        # # output
        uniref_count_with_info_df.to_csv(os.path.join(output_folder, f"{self.gene_name}_uniref_count_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=False)
        cpm_with_info_df.to_csv(os.path.join(output_folder, f"{self.gene_name}_uniref_cpm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=False)
        rpkm_with_info_df.to_csv(os.path.join(output_folder, f"{self.gene_name}_uniref_rpkm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=False)
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
        uniref_count_df_species_level.to_csv(os.path.join(output_folder, f"{self.gene_name}_species_count_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=True)
        cpm_df_species_level.to_csv(os.path.join(output_folder, f"{self.gene_name}_species_cpm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=True)
        rpkm_df_species_level.to_csv(os.path.join(output_folder, f"{self.gene_name}_species_rpkm_q{identical_percent_threshold}_len{matched_length_threshold}.csv"), index=True)


def count_reads_shell(fastq_path):
    # fastq_path = "path/example_R1.fq.gz"
    if fastq_path.endswith(".gz"):
        cmd = f"zcat {fastq_path} | wc -l"
    else:
        cmd = f"cat {fastq_path} | wc -l"
    # result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
    line_count = int(result.stdout.strip())
    read_count = int(line_count/4)  # Each read has 4 lines
    return read_count



