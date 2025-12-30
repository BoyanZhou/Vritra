"""
module add python/cpu/3.7.2

Script Name: build_diamond_db.py
Author: Boyan Zhou
Date: 2025-05-29
Description: []

Function Overview:
    1. align_samples_from_csv():
        - Purpose: multi samples in a csv file
        - Dependencies:

    1. align_samples_by_diamond():
        - Purpose: per sample
        - Dependencies:

"""

import os


def align_samples_from_csv(gene_prefix, sample_info_path, diamond_db_path, output_dir, thread, max_hit, id_thresh, query_cover, my_logger):
    # {sample: ["PATH/XXX_R1.fq", "PATH/XXX_R2.fq"]}
    sample_fq_dict = sample_fq_dict_from_csv(sample_info_path, my_logger)

    for sample_id, fq_list in sample_fq_dict.items():
        align_sample_by_diamond(gene_prefix, sample_id, fq_list, diamond_db_path, output_dir, thread, max_hit,
                                 id_thresh, query_cover, my_logger)


def sample_fq_dict_from_csv(sample_info_path, my_logger):
    # 1. get sample fq-path dict
    sample_fq_dict = {}     # {sample: ["PATH/XXX_R1.fq", "PATH/XXX_R2.fq"]}
    with open(sample_info_path, "r") as sample_info_f:
        for line in sample_info_f:
            cols = line.strip().split(",")
            fq_path_list = cols[1].split(":")   # ["PATH/XXX.fq"], ["PATH/XXX_R1.fq", "PATH/XXX_R2.fq"]
            if len(fq_path_list) > 2:
                my_logger.info(f"The sample {cols[0]} has more than two fq samples {fq_path_list}. Skip this sample!")
            else:
                sample_fq_dict.update({cols[0]: fq_path_list})
    return sample_fq_dict


def align_sample_by_diamond(gene_prefix, sample_id, fq_list, diamond_db_path, output_dir, thread, max_hit, id_thresh, query_cover, my_logger):
    """
    --log of diamond This will create a file named diamond.log in your current working directory
    :param gene_prefix:
    :param sample_id:
    :param fq_list:
    :param diamond_db_path:
    :param output_dir:
    :param thread:
    :param max_hit:
    :param id_thresh:
    :param query_cover:
    :param my_logger:
    :return:
    """
    if not os.path.exists(output_dir):
        os.system(f"mkdir -p {output_dir}")
    os.chdir(output_dir)

    ##############
    # SINGLE END #
    ##############
    if len(fq_list) == 1:
        # fq_list = ["PATH/XXX.fq"]
        if not os.path.exists(fq_list[0]):
            my_logger.info(f"Warning! For sample {sample_id}, fq file {fq_list[0]} not exist! Skip this sample.")
        else:
            log_filename = f"{sample_id}_{gene_prefix}_diamond.log"
            diamond_report_path = os.path.join(output_dir, f"{sample_id}_{gene_prefix}_diamond.output")
            diamond_command = f"diamond blastx -d {diamond_db_path} -p {thread} -q {fq_list[0]} -o {diamond_report_path}" \
                            f" --header --max-target-seqs {max_hit} " \
                            f"--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen " \
                            f"--id {id_thresh} --query-cover {query_cover} --log 2>{log_filename}"

            my_logger.info(diamond_command)
            os.system(diamond_command)

    ##############
    # PAIRED END #
    ##############
    elif len(fq_list) == 2:
        # fq_list = ["PATH/XXX_R1.fq", "PATH/XXX_R2.fq"]
        if not (os.path.exists(fq_list[0]) and os.path.exists(fq_list[1])) :
            my_logger.info(f"Warning! For sample {sample_id}, fq files {fq_list} not exist! Skip this sample.")
        else:
            for i, fq_path in enumerate(fq_list):
                # temp folder for log file per run (remove after getting the log file)
                log_filename = f"{sample_id}_{gene_prefix}_R{i+1}_diamond.log"
                diamond_report_path = os.path.join(output_dir, f"{sample_id}_{gene_prefix}_diamond_R{i+1}.output")
                diamond_command = f"diamond blastx -d {diamond_db_path} -p {thread} -q {fq_path} -o {diamond_report_path}" \
                                f" --header --max-target-seqs {max_hit} " \
                                f"--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen " \
                                f"--id {id_thresh} --query-cover {query_cover} --log 2>{log_filename}"

                my_logger.info(diamond_command)
                os.system(diamond_command)

    else:
        my_logger.info(f"Warning! For sample {sample_id}, fq files are {fq_list}. Neither single nor paired.")
