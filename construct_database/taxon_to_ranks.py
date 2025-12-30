import os
from ete3 import NCBITaxa


"""
# metagenomic requires, then diamond v2.0.15.153 will be loaded automatically
module add python/cpu/3.7.2
# !!! but for this script, we must use python 3.6.5, because only under this python, ete3 is available 
module add python/cpu/3.6.5

"""

###########################
# Output of blastx result #
###########################
# Notice that: some have taxID, some do not. This taxID is from InterPro database, should be NCBI taxon ID
"""
# Fields: Query ID, Subject ID, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query,
End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
SRR5155786.36178        A0A0L6ZX88|unreviewed|D-phenylhydantoinase|taxID:562    100     33      0       0       1       99      358     390     1.76e-19    71.6
SRR5155786.81704        A0A0H3EP60|unreviewed|D-phenylhydantoinase|taxID:685038 100     30      0       0       1       90      90      119     3.92e-19    70.5

# DIAMOND v2.0.15. http://github.com/bbuchfink/diamond
# Invocation: diamond blastx -d /gpfs/data/lilab/home/zhoub03/Uric_Acid_Lama/InterPro_7_genes/IPR050028-xdhAC_nonredundant.dmnd -p 8 -q /gpfs/data/lilab/home
/zhoub03/Uric_Acid_Lama/datasets/USman/10042055/dna/SRR5155786.fastq -o /gpfs/data/lilab/home/zhoub03/Uric_Acid_Lama/USman_diamond_result/10042055/dna/SRR515
5786_xdhAC_diamond.output --header --max-target-seqs 1
# Fields: Query ID, Subject ID, Percentage of identical matches, Alignment length, Number of mismatches, Number of gap openings, Start of alignment in query,
 End of alignment in query, Start of alignment in subject, End of alignment in subject, Expected value, Bit score
SRR5155786.9917 A0A2T5CF16|unreviewed|Aldehyde  62.1    29      11      0       14      100     639     667     1.52e-08        43.1
SRR5155786.12584        A0A6N2R7Z8|unreviewed|Aldehyde  57.1    21      9       0       81      19      15      35      5.46e-05        32.7
SRR5155786.17121        A0A5R9LIU7|unreviewed|Aldehyde  100     18      0       0       3       56      643     660     1.11e-08        43.5
"""


# function3: get taxon-id to taxon ranks
def taxon_to_ranks(taxon_id, desired_ranks):
    """
    :param taxon_id: taxon_id = 562
    :param ncbi: ncbi = NCBITaxa()
    :param desired_ranks: ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    :return: "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Lachnoclostridium|s__Lachnoclostridium phytofermentans"
    """
    ncbi = NCBITaxa()
    rank_abbr_dict = {'superkingdom': "d", 'phylum': "p", 'class': "c", 'order': "o", 'family': "f",
                      'genus': "g", 'species': "s"}
    try:
        lineage = ncbi.get_lineage(taxon_id)     # lineage to taxon_id, [1, 131567, 2, 1224, 1236, 91347, 543, 561, 562]
    except ValueError:
        print(f"{taxon_id} taxid not found")
        return None, None
    names = ncbi.get_taxid_translator(lineage)  # {1: 'root', 2: 'Bacteria', 543: 'Enterobacteriaceae', 561: 'Escherichia', 562: 'Escherichia coli'}
    lineage2ranks = ncbi.get_rank(names)    # {1: 'no rank', 2: 'superkingdom', 543: 'family'}
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    # {'superkingdom': 'Bacteria', 'phylum': 'Proteobacteria', 'class': 'Deltaproteobacteria', 'order': 'Myxococcales',
    # 'family': 'Nannocystaceae', 'genus': 'Nannocystis', 'species': None}
    rank_name_dict = {rank: names.get(ranks2lineage.get(rank, None), None) for rank in desired_ranks}
    full_ranks = "|".join([f"{rank_abbr_dict[rank]}__{rank_name_dict[rank]}" for rank in desired_ranks if rank_name_dict[rank]])   # "d__Bacteria|p__Firmicutes"
    return full_ranks, rank_name_dict


# main function1:
def blastx_to_taxon_multi_samples(sample_id_list, blastx_output_path_list, protein_seq_fas, output_path, logger, identical_percent_threshold=70):
    # 1. get the protein to taxon dict, protein: taxon
    # # A0A3A5JXE0|unreviewed|D-phenylhydantoinase|taxID:82991, {"A0A3A5JXE0": 562}
    protein_taxon_dict = protein_to_taxon(protein_seq_fas)

    # 2. for each sample, get the count of each taxon
    taxon_id_unique_set = set()        # unique taxon id appears in at least one sample
    sample_taxon_count_dict = {}    # {sample_id: {562:35, ...}}
    sample_id_with_result_list = []
    for sample_id, blastx_output_path in zip(sample_id_list, blastx_output_path_list):
        logger.info(f"Read sample {sample_id}'s file {blastx_output_path}")
        if not os.path.exists(blastx_output_path):
            logger.info(f"sample {sample_id}'s file {blastx_output_path} not exists")
            continue
        taxon_count_dict = blastx_to_taxon_count(blastx_output_path, protein_taxon_dict, identical_percent_threshold)
        taxon_id_unique_set.update(set(taxon_count_dict.keys()))
        sample_taxon_count_dict.update({sample_id: taxon_count_dict})
        sample_id_with_result_list.append(sample_id)
    logger.info(f"All unique taxon id are {taxon_id_unique_set}")

    # 3. get the full ranks of each taxon id (int)
    # "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Lachnoclostridium|s__Lachnoclostridium phytofermentans"
    desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxon_id_ranks_dict = {taxon_id: taxon_to_ranks(taxon_id, desired_ranks)[0] for taxon_id in list(taxon_id_unique_set) if taxon_to_ranks(taxon_id, desired_ranks)}
    ranks_taxon_id_dict = {value: key for key, value in taxon_id_ranks_dict.items()}

    # 4. output the counts of all samples
    with open(output_path, "w") as output_f:
        # write the header: sample ids
        output_f.write("Taxon_ranks" + "\t" + "\t".join(sample_id_with_result_list) + "\n")
        for full_rank in sorted(ranks_taxon_id_dict.keys()):
            # "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|..."
            taxon_id_temp = ranks_taxon_id_dict[full_rank]      # int, 562
            # if this taxon in this sample, get its count number, otherwise get 0
            line_of_count = [str(sample_taxon_count_dict[sample_id].get(taxon_id_temp, 0)) for sample_id in sample_id_with_result_list]
            output_f.write(full_rank + "\t" + "\t".join(line_of_count) + "\n")


# main function1:
def blastx_to_taxon_count(blastx_output_path, protein_taxon_dict, identical_percent_threshold=70):
    """
    :param blastx_output_path:
    :param protein_taxon_dict:
    :param identical_percent_threshold: threshold of identical aa percent, default 70
    :return: reads count of each taxon
    """
    # 1. read the unique protein id that the read mapped to, like "A0A6N2R7Z8"; get the count of each id appeared
    protein_count_dict = {}     # {"A0A6N2R7Z8": 35, ...}
    with open(blastx_output_path, "r") as blastx_f:
        for line in blastx_f:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            protein_id = cols[1].split("|")[0]  # like "A0A6N2R7Z8"
            identical_percent = float(cols[2])
            # only consider the sequences passing the threshold
            if identical_percent < identical_percent_threshold:
                continue
            if protein_id in protein_count_dict:
                protein_count_dict[protein_id] += 1
            else:
                protein_count_dict.update({protein_id: 1})

    # 3. summarize the count under each existing taxon id (different protein may be from same taxon ID)
    taxon_count_dict = {}   # {562: 35}
    for protein_id, count in protein_count_dict.items():
        if protein_id not in protein_taxon_dict:
            print(f"Warning! {protein_id} is not in the protein_taxon_dict!")
        else:
            taxon_id = protein_taxon_dict[protein_id]   # int, 562
            if taxon_id not in taxon_count_dict:
                taxon_count_dict.update({taxon_id: protein_count_dict[protein_id]})
            else:
                taxon_count_dict[taxon_id] += protein_count_dict[protein_id]
    return taxon_count_dict


# function2: get protein-id to taxon-id dict
def protein_to_taxon(protein_seq_fas):
    protein_taxon_dict = {}     # # A0A3A5JXE0|unreviewed|D-phenylhydantoinase|taxID:82991, {"A0A3A5JXE0", 562}
    with open(protein_seq_fas, "r") as protein_f:
        for line in protein_f:
            if line.startswith(">"):
                cols = line.strip().split("|")
                protein_id = cols[0][1:].strip()
                if "taxID:" not in cols[-1]:
                    print(f"There is no taxID for line {line}!")
                else:
                    taxID = int(cols[-1].split(":")[-1].strip())
                    protein_taxon_dict.update({protein_id: taxID})
    return protein_taxon_dict


if __name__ == "__main__":
    pass

