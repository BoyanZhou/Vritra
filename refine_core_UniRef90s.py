"""
Refine the constructed core UniRef90s database according to self domain knowledge, provided in txt file

module add diamond/0.9.18
module add python/cpu/3.6.5


"""

import os
import re
import pandas as pd
from Bio import SeqIO
import log_setup
import argparse

# import get_seq_by_id_from_uniref as gsbifu
# import detect_community as detect_com
# import download_uniref100_from_uniref90 as download_100_from_90


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


# extract the function of protein
def extract_description(text):
    """
    can deal with both uniref90 and UniProtKB
    'UniRef90_A0A117EBN7 Bile acid-coenzyme A ligase n=24 Tax=Streptomyces TaxID=1883 RepID=A0A117EBN7_STRSC\n'
    "tr|A0A021XAR7|A0A021XAR7_9HYPH Xanthine dehydrogenase molybdenum-binding subunit OS=Shinella sp. DD12 OX=1410620 GN=xdhA PE=3 SV=1"
    :param text:
    :return:
    """
    match = re.match(r'^(\S+)\s+(.*?)(?=\s+\S+=)', text)
    if match:
        record_id, description = match.groups()
        return record_id, description.strip()
    return None, None


def extract_SeqIO_given_func(original_fas_path, func_list, db_name="uniprot"):
    """
    'UniRef90_A0A8T7LD91': SeqRecord(seq=Seq('MLITDGTLITNETPNRILSDHALYIEDSVIQAVGRTDELLARYPQAERLNARGQ...TLL'),
    id='UniRef90_A0A8T7LD91', name='UniRef90_A0A8T7LD91',
    description='UniRef90_A0A8T7LD91 Putative aminohydrolase SsnA n=1 Tax=Chloroflexota bacterium TaxID=2026724 RepID=A0A8T7LD91_UNCCH', dbxrefs=[])
    :param original_fas_path:
    :param func_list:
    :param db_name:
    :return:
    """
    func_set = set(func_list)
    func_count = {i: 0 for i in func_list}
    # extract the seq records which match the given function list
    extracted_seq_record_dict = {}  # standard seq record from SeqIO

    original_seq_record_dict = SeqIO.to_dict(SeqIO.parse(original_fas_path, "fasta"))
    for seq_id, seq_record in original_seq_record_dict.items():
        seq_id_temp, seq_func = extract_description(seq_record.description)
        # if the function of seq match with the function target set
        if seq_func in func_set:
            func_count[seq_func] += 1
            extracted_seq_record_dict.update({seq_id: seq_record})
    return extracted_seq_record_dict, func_count


def get_func_name_from_file(func_name_file_path):
    """
    :param func_name_file_path: --expand or --identified or --revert
    :return:
    """
    if func_name_file_path is None:
        print(f"The file path of one of the function names is not specified. It is okay.")
        return []
    else:
        func_name_list = []
        with open(func_name_file_path, "r") as func_f:
            for line in func_f:
                func_name_list.append(line.strip().split("\t")[0])
        return func_name_list


def get_seq_fas_path(search_dir, input_pattern):
    """
    Raise error if no match or more than one match
    :param search_dir:
    :param input_pattern:
    :return: path of file matching the pattern
    """
    # Find matching files
    matching_files = [f for f in os.listdir(search_dir) if input_pattern.match(f)]

    # Raise error if none or more than one file matches
    if len(matching_files) == 0:
        raise FileNotFoundError(f"Did not find one file matching pattern {input_pattern.pattern}")
    elif len(matching_files) > 1:
        raise FileNotFoundError(f"Found more than one file matching pattern {input_pattern.pattern},"
                                f"{matching_files}")
    matched_file_path = os.path.join(search_dir, matching_files[0])
    return matched_file_path


def refine_core_uniref90_dataset(gene_prefix, dataset_construction_dir, dataset_refined_dir, func_revert_path, func_identified_path,
                                 func_expand_path, logger_path):
    logger = log_setup.setup_log_file(logger_path)
    if not os.path.exists(dataset_refined_dir):
        os.system(f"mkdir -p {dataset_refined_dir}")

    core_uniref90_raw_path = os.path.join(dataset_refined_dir, f"{gene_prefix}_refined_core_UniRef90_raw.fas")
    core_uniref90_filtered_path = os.path.join(dataset_refined_dir, f"{gene_prefix}_refined_core_UniRef90_filtered.fas")    # after length filter
    core_uniref90_set = set()   # not recorded in this set are put into refined_sponge
    refined_sponge_path = os.path.join(dataset_refined_dir, f"{gene_prefix}_sponge_UniRef90.fas")

    core_uniref90_raw_file = open(core_uniref90_raw_path, "w")  # write seq record to fas path
    refined_func_count = {}

    ####################
    # func_revert_path #
    ####################
    # xdhA_refined_labeled_filtered_seqs.fas + ygeY_quasi_labeled_search1_filtered_seqs.fas
    pattern1 = re.compile(rf"^{re.escape(gene_prefix)}_refined_labeled.*?_filtered_seqs\.fas$")
    refined_labeled_fas_path = get_seq_fas_path(dataset_construction_dir, pattern1)
    pattern2 = re.compile(rf"^{re.escape(gene_prefix)}_quasi_labeled.*?_filtered_seqs\.fas$")
    quasi_labeled_fas_path = get_seq_fas_path(dataset_construction_dir, pattern2)
    func_to_search_in_pre_labeled = get_func_name_from_file(func_revert_path)
    if len(func_to_search_in_pre_labeled) > 0:
        refined_labeled_seq_record_dict, refined_labeled_func_count = extract_SeqIO_given_func(refined_labeled_fas_path, func_to_search_in_pre_labeled, db_name="uniprot")
        quasi_labeled_seq_record_dict, quasi_labeled_func_count = extract_SeqIO_given_func(quasi_labeled_fas_path, func_to_search_in_pre_labeled, db_name="uniprot")
        core_uniref90_set.update(set(refined_labeled_seq_record_dict.keys()))
        core_uniref90_set.update(set(quasi_labeled_seq_record_dict.keys()))
        for seq_id_temp, seq_record_temp in refined_labeled_seq_record_dict.items():
            SeqIO.write(seq_record_temp, core_uniref90_raw_file, "fasta")
        for k, v in refined_labeled_func_count.items():
            refined_func_count[k] = refined_func_count.get(k, 0) + v

        for seq_id_temp, seq_record_temp in quasi_labeled_seq_record_dict.items():
            SeqIO.write(seq_record_temp, core_uniref90_raw_file, "fasta")
        for k, v in quasi_labeled_func_count.items():
            refined_func_count[k] = refined_func_count.get(k, 0) + v

    else:
        logger.info(f"No function given to retrieve from pre-labeled sequences.")

    ########################
    # func_identified_path #
    ########################
    # ygeY_pre_labeled_and_new_labeled_raw.fas
    pattern3 = re.compile(rf"^{re.escape(gene_prefix)}_pre_labeled_and_new_labeled_filtered\.fas$")
    identified_fas_path = get_seq_fas_path(dataset_construction_dir, pattern3)
    func_to_search_in_identified = get_func_name_from_file(func_identified_path)
    if len(func_to_search_in_identified) > 0:
        identified_seq_record_dict, identified_func_count = extract_SeqIO_given_func(identified_fas_path, func_to_search_in_identified, db_name="uniprot")
        core_uniref90_set.update(set(identified_seq_record_dict.keys()))
        for seq_id_temp, seq_record_temp in identified_seq_record_dict.items():
            SeqIO.write(seq_record_temp, core_uniref90_raw_file, "fasta")
        for k, v in identified_func_count.items():
            refined_func_count[k] = refined_func_count.get(k, 0) + v
    else:
        logger.info(f"No function given to retrieve from identified sequences.")

    ####################
    # func_expand_path #
    ####################
    # ygeY_labeled_search1_connected_and_its_connected_and_labeled_seq.fas
    pattern4 = re.compile(rf"^{re.escape(gene_prefix)}_labeled.*?_connected_and_its_connected_and_labeled_seq\.fas$")
    expand_fas_path = get_seq_fas_path(dataset_construction_dir, pattern4)
    func_to_search_in_expand = get_func_name_from_file(func_expand_path)
    if len(func_to_search_in_expand) > 0:
        expand_seq_record_dict, expand_func_count = extract_SeqIO_given_func(expand_fas_path, func_to_search_in_expand, db_name="uniprot")
        core_uniref90_set.update(set(expand_seq_record_dict.keys()))
        for seq_id_temp, seq_record_temp in expand_seq_record_dict.items():
            SeqIO.write(seq_record_temp, core_uniref90_raw_file, "fasta")
        for k, v in expand_func_count.items():
            refined_func_count[k] = refined_func_count.get(k, 0) + v

    else:
        logger.info(f"No function given to retrieve from expand sequences.")

    core_uniref90_raw_file.close()

    # filter core_uniref90_fas by length
    filter_fas_by_len(core_uniref90_raw_path, core_uniref90_filtered_path, logger, len_upper=0.2, len_lower=0.2, by_sd=False)

    #########################
    # output refined sponge #
    #########################
    # expand_fas_path contains pre-labeled, t1, and t2
    refined_sponge_f = open(refined_sponge_path, "w")
    seq_record_dict = SeqIO.to_dict(SeqIO.parse(expand_fas_path, "fasta"))
    for seq_id_temp, seq_record_temp in seq_record_dict.items():
        if seq_id_temp not in core_uniref90_set:
            SeqIO.write(seq_record_temp, refined_sponge_f, "fasta")
    refined_sponge_f.close()

    print(f"The path of sponge sequences is {refined_sponge_path}\n"
          f"The path of finally refined core UniRef90 sequences is {core_uniref90_filtered_path}")
    logger.info(f"The path of sponge sequences is {refined_sponge_path}\n"
                f"The path of finally refined core UniRef90 sequences is {core_uniref90_filtered_path}")


def get_args():
    """Parses command-line arguments and returns them."""

    # Define parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-g", "--gene-prefix", type=str, required=True, help="Prefix for gene files")
    parent_parser.add_argument("-e", "--expand", type=str, help="Path to the file. For seqs annotated as functions in this file,"
                                                                "pre-labeled, new labled, and new unlabeled seqs should be included.")
    parent_parser.add_argument("-i", "--identified", type=str, help="Path to the file. For seqs annotated as functions in this file,"
                                                                    "only pre-labeled and new labled seqs should be included (discard new unlabeled seqs).")
    parent_parser.add_argument("-r", "--revert", type=str, help="Path to the file. For seqs annotated as functions in this file,"
                                                                "only pre-labeled seqs should be included (discard new labeled and new unlabeled seqs).")
    parent_parser.add_argument("-c", "--construct-dir", type=str, required=True,
                               help="Directory of dataset construction in last step, like PATH/dataset_construction_by_annotated_seqs")
    parent_parser.add_argument("-o", "--output-dir", type=str, required=True,
                               help="Directory of total output for processing")
    parent_parser.add_argument("-L", "--logger-path", type=str, required=True, help="Path to save log files")

    # Main parser
    parser = argparse.ArgumentParser(description="Curate the labeled results by 'reapply' or 'refine'")
    subparsers = parser.add_subparsers(dest="mode", help="Mode of manually curate")

    # Subparser for 'Reapply' mode (shared + specific args)
    parser_reapply = subparsers.add_parser("reapply", parents=[parent_parser],
                                             help="Repeat the label extension using current labeled results")
    # parser_reapply.add_argument("-l", "--annotated-seq-path", type=str, required=True, help="Path to annotated sequence file")

    # Subparser for 'Refine' mode (only refine the current labeled results, not repeat label processing)
    parser_refine = subparsers.add_parser("refine", parents=[parent_parser], help="Refine the current labeled results")
    # parser_refine.add_argument("-e", "--e_value", type=float, default=1e-3, help="E-value threshold for PSI-BLAST search (default: 1e-3)")

    args1 = parser.parse_args()

    if args1.mode is None:
        parser.error("A subcommand is required. Use -h for help.")
    return args1


if __name__ == "__main__":
    args = get_args()

    if args.mode == "reapply":
        pass

    elif args.mode == "refine":
        """
        Example:
        python /gpfs/data/lilab/home/zhoub03/Liisa_Beta_glucuronidase/Vritra_v2_20250528/refine_core_UniRef90s.py 
        refine 
        -g ${gene_name} 
        -c ${gene_fir}/dataset_construction_by_annotated_seqs 
        -o ${gene_fir}/dataset_construction_refined 
        -e ${gene_fir}/dataset_construction_by_annotated_seqs/${gene_name}_search1_expand.txt 
        -i ${gene_fir}/dataset_construction_by_annotated_seqs/${gene_name}_search1_identified.txt 
        -L core_uniref90_function_refined_0523.log
        Output:
        3_dataset_construction_refined/ssnA_sponge_UniRef90.fas
        3_dataset_construction_refined/ssnA_refined_core_UniRef90_raw.fas
        3_dataset_construction_refined/ssnA_refined_core_UniRef90_filtered.fas
        """
        refine_core_uniref90_dataset(args.gene_prefix, args.construct_dir, args.output_dir, args.revert,
                                     args.identified, args.expand, args.logger_path)

    print("End of core UniRef90 refinement.")
