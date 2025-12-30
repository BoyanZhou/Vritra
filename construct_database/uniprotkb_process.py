from Bio import SeqIO
import construct_database.taxon_to_ranks as taxon_to_ranks
import re
import recluster_species_UniRef
import os
import pandas as pd


def uniprotkb_simplified(uniprotkb_fas_path, refined_species_fas_path, correspondence_csv_path, uniprotkb_intermediate_dir):
    """
    :param uniprotkb_fas_path: has been filtered by length
    :param refined_species_fas_path: output for downstream analysis
    :param correspondence_csv_path: correspondence between uniprotkb id and taxon rank, "XXXX_uniprot_correspondence.csv"
    :param uniprotkb_intermediate_dir:
    :return:
    """
    """
    'tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM': SeqRecord(seq=Seq('MEIGRSRHRVDAWSKVTGEAKYTADLFPDNCLTAKVIRSTIANGRVLSMDTREA...
    AYV'), id='tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM', name='tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM', 
    description='tr|A0A9D1HGB9|A0A9D1HGB9_9FIRM Xanthine dehydrogenase molybdenum-binding subunit XdhA 
    OS=Candidatus Onthocola gallistercoris OX=2840876 GN=xdhA PE=4 SV=1', dbxrefs=[])
    """
    uniprotkb_seq_record_dict = SeqIO.to_dict(SeqIO.parse(uniprotkb_fas_path, "fasta"))
    # 1. extract taxon rank for each record to "uniprotkb_ranks_dict" and "species_uniref_id_dict"
    print(f"UniProtKB processing step1 ...")
    pattern = re.compile(r'OX=(\d+)')
    desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    uniprotkb_ranks_dict = {}   # {uniref_id: {"species": , "taxon_id", "full_ranks": , "rank_name_dict": }}
    species_uniref_id_dict = {}     # {species: [uniref_id1, uniref_id2,]}
    for key, value in uniprotkb_seq_record_dict.items():
        match = pattern.search(value.description)
        if match:
            taxon_id = match.group(1)     # 2840876
            if not taxon_id:
                continue
            taxon_id = int(taxon_id)
            # full_ranks: "d__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Lachnoclostridium|s__Lachnoclostridium phytofermentans"
            full_ranks, rank_name_dict = taxon_to_ranks.taxon_to_ranks(taxon_id, desired_ranks)
            if rank_name_dict:
                if rank_name_dict['species']:
                    species_name_temp = rank_name_dict['species']
                    uniprotkb_ranks_dict.update({key: {"species": species_name_temp, "taxon_id": taxon_id,
                                                       "full_ranks": full_ranks, "rank_name_dict": rank_name_dict}})
                    if species_name_temp in species_uniref_id_dict:
                        species_uniref_id_dict[species_name_temp].append(key)
                    else:
                        species_uniref_id_dict.update({species_name_temp: [key]})
    # ------------------------------------------------------------------------------------------------------------------
    # 2. centroid sequences for each recorded species, stored in "centroid_uniref_id_list"
    if not os.path.exists(uniprotkb_intermediate_dir):
        os.system(f"mkdir -p {uniprotkb_intermediate_dir}")
    centroid_uniref_id_list = []
    for species_name_temp, uniref_id_list in species_uniref_id_dict.items():
        print(f"UniProtKB processing step2, species: {species_name_temp}")
        seq_record_dict_temp = {i: uniprotkb_seq_record_dict[i] for i in uniref_id_list}
        species_name_temp = species_name_temp.replace(" ", "")
        centroid_seq_uniref_id = recluster_species_UniRef.get_centroid.centroid_seq(seq_record_dict_temp,
                                                                                    species_name_temp,
                                                                                    uniprotkb_intermediate_dir)
        centroid_uniref_id_list.append(centroid_seq_uniref_id)
    # ----------------------------------------
    # 3. output refined species fas
    extracted_fas_file = open(refined_species_fas_path, "w")
    print(f"All selected {len(centroid_uniref_id_list)} uniref_ids are {centroid_uniref_id_list}")
    for seq_id_temp in centroid_uniref_id_list:
        SeqIO.write(uniprotkb_seq_record_dict[seq_id_temp], extracted_fas_file, "fasta")
    extracted_fas_file.close()
    # ----------------------------------------
    # 4. output correspondence_csv file for UniProtKB
    correspondence_info_list = []
    # # extract the info of each represented species
    for seq_id_temp in centroid_uniref_id_list:
        protein_len = len(uniprotkb_seq_record_dict[seq_id_temp].seq)
        # "UniRef90_ID,UniRef100_ID,Protein_length,Taxon_ID,Taxon_ranks"
        species_info_list = [seq_id_temp, seq_id_temp, protein_len, uniprotkb_ranks_dict[seq_id_temp]["taxon_id"],
                             uniprotkb_ranks_dict[seq_id_temp]["full_ranks"]]
        # "superkingdom,phylum,class,order,family,genus,species"
        ranks_list = [uniprotkb_ranks_dict[seq_id_temp]["rank_name_dict"].get(i, "") for i in desired_ranks]
        species_info_list.extend(ranks_list)
        species_name_temp = uniprotkb_ranks_dict[seq_id_temp]["species"]
        rep_species_summary = f"{species_name_temp}:{len(species_uniref_id_dict[species_name_temp])}"
        # 'rep_species, rep_species_summary, represented_uniref100, corresponding_uniref90, corresponding_species'
        species_info_list.extend([species_name_temp, rep_species_summary, seq_id_temp, seq_id_temp, species_name_temp])
        correspondence_info_list.append(species_info_list)

    correspondence_info_df = pd.DataFrame(
            correspondence_info_list,
            columns=["UniRef90_ID","UniRef100_ID","Protein_length","Taxon_ID","Taxon_ranks", "superkingdom", "phylum",
                    "class", "order", "family","genus", "species", "rep_species", "rep_species_summary",
                     "represented_uniref100", "corresponding_uniref90", "corresponding_species"]
        )
    correspondence_info_df.to_csv(correspondence_csv_path, index=False, header=True)


if __name__ == "__main__":
    pass

