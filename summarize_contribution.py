"""
# metagenomic requires, then diamond v2.0.15.153 will be loaded automatically
module add python/cpu/3.7.2

#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --job-name=kraken_boyan
#SBATCH --mem-per-cpu=5G
#SBATCH --time=72:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
"""

import os
import json
import pandas as pd
import subprocess
import log_setup
import argparse
import summarize_contribution


# main function
def diamond_summarize(gene_name, core_uniref100_fas_path, seq_taxonomy_mapping_csv_path, sample_info_dict,
                      output_folder, summary_folder, read_count_report_path, logger):
    diamond_summary = summarize_contribution.diamond_summarize.DiamondSum(gene_name, core_uniref100_fas_path, seq_taxonomy_mapping_csv_path, sample_info_dict, output_folder, logger)
    # diamond_summary.diamond_align(thread=4, max_target=5, whether_overwrite=False)
    diamond_summary.diamond_result_summary(read_count_report_path, summary_folder, identical_percent_threshold=70, matched_length_threshold=25)



# deprecated now
def build_dict_from_seq_taxonomy_mapping(mapping_csv_path):
    """
    :param mapping_csv_path = "FRC_core_uniref100_correspondence.csv"
    :return:
    """
    seq_taxonomy_df = pd.read_csv(mapping_csv_path)
    # seq_taxonomy_df.columns
    # ['UniRef90_ID','UniRef100_ID','Protein_length','Taxon_ID','Taxon_ranks','superkingdom','phylum','class',
    #  'order','family','genus','species','rep_species','rep_species_summary','represented_uniref100',
    #  'corresponding_uniref90', 'corresponding_species']
    pass


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


def blastx_to_taxon_count_unipro(blastx_output_path, identical_percent_threshold=70, matched_length_threshold=25):
    """
    This is for the results of diamond using UniPro database (not InterPro)
    For each unique read, only count the first record
    :param blastx_output_path:
    :param identical_percent_threshold: threshold of identical aa percent, default 70
    :param matched_length_threshold: the aa length >= 25 (that is >=75 in bp)
    :return: reads count dict of each UniRef, {"UniRef100_A0A2D6BS33": 10, "UniRef100_UPI000C7BF04C": 16}
    """
    # 1. read the unique protein id that the read mapped to, like "A0A6N2R7Z8"; get the count of each id appeared
    protein_count_dict = {}     # {"UniRef100_": 35, ...}
    last_read_id = ""
    with open(blastx_output_path, "r") as blastx_f:
        for line in blastx_f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            # if the read has been recorded, skip other records (only the best hit will be recorded)
            if cols[0] == last_read_id:
                continue
            last_read_id = cols[0]
            protein_id = cols[1]  # like "UniRef100_UPI002DDAB8EB"
            identical_percent = float(cols[2])
            aligned_length = int(cols[3])
            # only consider the sequences passing the threshold
            if identical_percent < identical_percent_threshold:
                continue
            if aligned_length < matched_length_threshold:
                continue
            if protein_id in protein_count_dict:
                protein_count_dict[protein_id] += 1
            else:
                protein_count_dict.update({protein_id: 1})
    return protein_count_dict


#################################
# read uniref taxonomy csv file #
#################################
# not used in this
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


def get_args():
    """Parses command-line arguments and returns them."""

    # Define parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-g", "--gene-prefix", type=str, required=True, help="Prefix for gene files")
    parent_parser.add_argument("-o", "--output-dir", type=str, required=True, help="Directory of output ")
    parent_parser.add_argument("-L", "--logger-path", type=str, required=True, help="Path to save log files")

    # Main parser
    parser = argparse.ArgumentParser(description="Build diamond_db, align sequence by diamond, and summarize.")
    subparsers = parser.add_subparsers(dest="mode", help="Mode of operation")

    # Subparser for build 'diamond_db' mode
    parser_diamond = subparsers.add_parser("build_diamond_db", parents=[parent_parser],
                                           help="Build the diamond database for target gene "
                                                "by core Uniref100 and sponge sequences.")
    parser_diamond.add_argument("-c", "--core-rep-uniref100", type=str, required=True,
                                help="Path to the core representative UniRef100 "
                                     "generated by 'recluster_species_UniRef.py'.")
    parser_diamond.add_argument("-s", "--sponge-fas-path", type=str,
                                help="Path to the FASTA file of sponge Sequence, "
                                     "generated by 'refine_core_UniRef90s.py'.")
    parser_diamond.add_argument("--uniprotkb-dir", type=str,
                                help="Path to the directory of UniProtKB generated by 'construct_database.py'.")

    # Subparser for 'align' sequence by diamond mode
    parser_align = subparsers.add_parser("align", parents=[parent_parser],
                                         help="Align raw sequence to finalized database by "
                                              "'build_diamond_db' in last step.")
    parser_align.add_argument("-d", "--diamond-db", type=str, required=True,
                              help="Path to the diamond database (.dmnd) built in last step by 'build_diamond_db'")
    parser_align.add_argument("-t", "--sample-table", type=str,
                              help="For multiple samples: not used together with --sample, --fastq \n"
                                   "Path to the csv table (separated by ',') containing sample information. "
                                    "The first column is the sample ID and the second column is the path to fq or "
                                    "fq.gz. For paired-end sample, the path to R1 and the path to R2 should be "
                                    "separated by ':' within the second column (PATH/XXX_R1.fq:PATH/XXX_R2.fq).")
    parser_align.add_argument("--sample", type=str, help="For single sample, sample ID.")
    parser_align.add_argument("--fastq", type=str,
                              help="For single sample, path to the fq or fq.gz. \n"
                                   "For paired-end sample, path to R1 and R2 should be separated by ':', "
                                   "like (PATH/XXX_R1.fq:PATH/XXX_R2.fq)")
    parser_align.add_argument("--thread", type=int, default=4, help="Number of thread for diamond align.")
    parser_align.add_argument("--max-hit", type=int, default=10, help="Number of max target seqs per query in diamond.")
    parser_align.add_argument("--identity-threshold", type=int, default=85,
                              help="Threshold of alignment identity (0~100).")
    parser_align.add_argument("--query-coverage", type=int, default=85,
                              help="Threshold of query covered percentage (0~100).")

    # Subparser for 'summarize' sequence by diamond mode
    parser_summarize = subparsers.add_parser("summarize", parents=[parent_parser],
                                             help="Align raw sequence to finalized database by "
                                                  "'build_diamond_db' in last step.")
    parser_summarize.add_argument("-t", "--sample-table", type=str, required=True,
                                  help="Path to the csv table (separated by ',') containing sample information. "
                                       "You can just used the same paramenter in last step ('align' -t). "
                                       "Or any csv file with the first column being the sample ID that matches the "
                                       "sample ID in last step ('align').")
    parser_summarize.add_argument("-f", "--uniref100-fas-path", type=str, required=True,
                                  help="Path to the FASTA file of UniRef100, generated by mode 'uniref90_to_100'")

    args1 = parser.parse_args()

    if args1.mode is None:
        parser.error("A subcommand is required. Use -h for help.")
    return args1


if __name__ == "__main__":
    args = get_args()
    logger = log_setup.setup_log_file(args.logger_path)
    # uniref_species = UniRefSpecies(args.uniref_taxonomy_path, args.uniref100_fas_path, logger)
    if args.mode == "build_diamond_db":
        """
        Example:
        echo "$gene"
        gene_dir=${parent_folder}/${gene}_2507
        python ${vritra_folder}summarize_contribution.py build_diamond_db \
        -g ${gene} \
        -o ${gene_dir}/4_finalized_database \
        -c ${gene_dir}/4_finalized_database/${gene}_core_representative_UniRef100.fas \
        -s ${gene_dir}/3_dataset_construction_refined/${gene}_sponge_UniRef90.fas \
        --uniprotkb-dir ${gene_dir}/2_dataset_construction_by_annotated_seqs_uniprotkb \
        -L ${gene_dir}/${gene}_construct_database_${date}.log
        """
        summarize_contribution.build_diamond_db.build_diamond_database(args.gene_prefix, args.core_rep_uniref100,
                                                                       args.sponge_fas_path, args.uniprotkb_dir,
                                                                       args.output_dir, logger)


    elif args.mode == "align":
        if args.sample_table is not None:
            """
            Example: Provide sample-ID and path to fastq in a table file, batch processing 
            """
            summarize_contribution.diamond_align.align_samples_from_csv(args.gene_prefix, args.sample_table,
                                                                        args.diamond_db, args.output_dir, args.thread,
                                                                        args.max_hit, args.identity_threshold,
                                                                        args.query_coverage, logger)
        elif (args.sample is not None) and (args.fastq is not None):
            """
            Example: Provide sample-ID and path to fastq, align one by one
            vritra_align = f"python {os.path.join(vritra_dir, 'summarize_contribution.py')} align " \
                       f"-g {gene_name} " \
                       f"--sample {srr_prefix} " \
                       f"--fastq {':'.join(input_fq_path_list)} " \
                       f"-d {db_path} " \
                       f"-o {gene_dir} " \
                       f"-L {logger_path}"            
            """
            fastq_path_list = args.fastq.split(":")
            summarize_contribution.diamond_align.align_sample_by_diamond(args.gene_prefix, args.sample, fastq_path_list,
                                                                         args.diamond_db, args.output_dir, args.thread,
                                                                         args.max_hit, args.identity_threshold,
                                                                         args.query_coverage, logger)
        else:
            print(f"Error! sample information is missing. Use -t or or '--sample and --fastq'")

    elif args.mode == "summarize":

        pass
        # diamond_summarize(args.gene_name, args.core_uniref100_fas_path, args.uniref_taxonomy_path, sample_info_dict, args.output_folder, args.summary_folder, read_count_report_path, logger)


