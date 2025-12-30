"""
module add python/cpu/3.7.2

Script Name: build_diamond_db.py
Author: Boyan Zhou
Date: 2025-07-28
Description: []

Function Overview:
    1. build_diamond_database():
        - Purpose:
        - Dependencies:
# 2025-07-28
compatible to UniProtKB
"""

import os
import pandas as pd
from Bio import SeqIO


#############################
# for mode build_diamond_db #
#############################
def build_diamond_database(gene_prefix, core_rep_uniref100_path, sponge_seq_path, uniprotkb_dir, output_dir, my_logger):
    if not os.path.exists(output_dir):
        os.system(f"mkdir -p {output_dir}")
    finalized_db_fas_path = os.path.join(output_dir, f"{gene_prefix}_finalized_database_no_sponge.fas")
    finalized_db_with_sponge_fas_path = os.path.join(output_dir, f"{gene_prefix}_finalized_database_with_sponge.fas")
    if uniprotkb_dir:
        # if UniProtKB was used, compare the correspondence obtained from uniprotkb and uniref
        finalized_db_dir = os.path.split(core_rep_uniref100_path)[0]
        uniprotkb_correspondence_path = os.path.join(uniprotkb_dir, f"{gene_prefix}_uniprotkb_taxon_ranks.csv")
        uniref_correspondence_path = os.path.join(finalized_db_dir, f"{gene_prefix}_core_uniref100_correspondence.csv")
        # core_uniref100_correspondence = pd.read_csv("xdhA_core_uniref100_correspondence.csv")
        # uniprotkb_correspondence = pd.read_csv("xdhA_uniprotkb_taxon_ranks.csv")
        uniprotkb_correspondence = pd.read_csv(uniprotkb_correspondence_path)
        core_uniref100_correspondence = pd.read_csv(uniref_correspondence_path)
        my_logger.info(f"Compare the correspondence files: \n"
                       f"{uniprotkb_correspondence_path}\n{uniref_correspondence_path}")

        # get the set of species in core_uniref100, core_uniref100_correspondence.columns
        species_core_uniref100 = set(core_uniref100_correspondence['species']) | set(
            core_uniref100_correspondence['rep_species'])
        # extract species not included in core_uniref100
        uniprotkb_correspondence_unique = uniprotkb_correspondence.loc[~uniprotkb_correspondence['species'].isin(species_core_uniref100)]
        uniprotkb_original_fas_path = os.path.join(uniprotkb_dir, f"{gene_prefix}_uniprotkb_species_refined.fas")
        uniprotkb_extracted_fas_path = os.path.join(finalized_db_dir, f"{gene_prefix}_core_representative_uniprotkb.fas")
        output_fas_given_seq(set(uniprotkb_correspondence_unique['UniRef90_ID']), uniprotkb_original_fas_path, uniprotkb_extracted_fas_path, reverse=False)
        # merge the correspondence of uniprotkb and uniref
        merged_correspondence_df = pd.concat([core_uniref100_correspondence, uniprotkb_correspondence_unique], ignore_index=True)
        merged_correspondence_path = os.path.join(output_dir, f"{gene_prefix}_finalized_correspondence.csv")
        merged_correspondence_df.to_csv(merged_correspondence_path, index=False, header=True)
        # output combined fas, with sponge and without sponge
        my_logger.info(f"UniProtKB incorporated ... ...\n"
                       f"cat {core_rep_uniref100_path} {uniprotkb_extracted_fas_path} > {finalized_db_fas_path}\n"
                       f"cat {core_rep_uniref100_path} {uniprotkb_extracted_fas_path} {sponge_seq_path} > {finalized_db_with_sponge_fas_path}")
        os.system(f"cat {core_rep_uniref100_path} {uniprotkb_extracted_fas_path} > {finalized_db_fas_path}")
        os.system(f"cat {core_rep_uniref100_path} {uniprotkb_extracted_fas_path} {sponge_seq_path} > {finalized_db_with_sponge_fas_path}")
    else:
        # if UniProtKB was not used
        my_logger.info(f"UniProtKB not used ... ...\n"
                       f"cat {core_rep_uniref100_path} {sponge_seq_path} > {finalized_db_with_sponge_fas_path}")
        os.system(f"cat {core_rep_uniref100_path} {sponge_seq_path} > {finalized_db_with_sponge_fas_path}")

    # make diamond db
    os.chdir(output_dir)
    mk_db_command = f"diamond makedb --in {finalized_db_with_sponge_fas_path} " \
                    f"-d {gene_prefix}_finalized_database_with_sponge"
    my_logger.info(mk_db_command)
    os.system(mk_db_command)


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