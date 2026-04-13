import construct_database.seqid_from_fas as sff
import construct_database.uniprot_api as u_api
import os
import json
from construct_database.taxon_to_ranks import taxon_to_ranks

def uniref90_to_100_taxon(uniref90_fas_path, output_json_path, output_taxon_ranks_path, fas_type='UniProt'):
    """
    :param uniref90_fas_path:
    :param fas_type: fas_type must be 'InterPro' or 'UniProt'.
    :param output_json_path:
    :param output_taxon_ranks_path:
    :return: json file of UniRef90 to UniRef100; taxonomy file
    """

    # step1
    uniref90_list_temp = sff.protein_id_from_fas(uniref90_fas_path, fas_type)

    # step2:

    # # uniref_api_response_dict is like {"Success 200: Success": [], "error 404: Resource not found": [],
    # #                               "error 403: Access forbidden": [], "error 500: Internal server error": [],
    # #                               "unexpected error": []}
    uniref90_full_dict, uniref_api_response_dict = u_api.uniref_full_dict(uniref90_list_temp, in_uniref_type='UniRef90')
    # Optionally save the JSON data to a file
    with open(output_json_path, "w") as json_file:
        json.dump(uniref90_full_dict, json_file, indent=4)

    # step3
    # uniref90_uniref100_dict = {}    # {uniref90_id: [uniref100_id1, uniref100_id2]}
    taxon_ranks_file = open(output_taxon_ranks_path, "w")
    desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxon_ranks_file.write("UniRef90_ID" + "," + "UniRef100_ID" + "," + "Protein_length" + "," + "Taxon_ID" + "," + "Taxon_ranks" + "," + ",".join(desired_ranks) + "\n")

    uniref90_no_query_result_list = []
    for key, value in uniref90_full_dict.items():
        # key is the uniref90 ID
        if value is None:
            # if no searched result for this UniRef90 ID (not found on UniProt)
            uniref90_no_query_result_list.append(key)
        else:
            " Record the main sequence first"
            # value['commonTaxon'] = {'scientificName': 'Pseudomonas frederiksbergensis', 'taxonId': 104087}
            taxon_id = value['representativeMember']['organismTaxId']
            protein_length = value['representativeMember']['sequenceLength']
            full_ranks, rank_name_dict = taxon_to_ranks(taxon_id, desired_ranks)

            # first check whether this taxon id is recognized in "from ete3 import NCBITaxa"
            if rank_name_dict is None:
                uniref90_no_query_result_list.append(taxon_id)
                print(f"The taxon {taxon_id} has no result in the NCBITaxa imported from ete3")
                continue

            full_ranks_list = []
            for rank in desired_ranks:
                if rank_name_dict[rank]:
                    full_ranks_list.append(rank_name_dict[rank])
                else:
                    full_ranks_list.append("")

            taxon_ranks_file.write(key + "," + value['representativeMember']["uniref100Id"] + "," + str(protein_length) + "," +
                                   str(taxon_id) + "," + full_ranks + "," + ",".join(full_ranks_list) + "\n")

            if "members" in value:
                # [{'memberIdType': 'UniParc', 'memberId': 'UPI002180F047', 'organismName': 'Klebsiella pneumoniae',
                #   'organismTaxId': 573, 'sequenceLength': 383, 'proteinName': 'diaminopropionate ammonia-lyase',
                #   'uniref90Id': 'UniRef90_W9B689', 'uniref100Id': 'UniRef100_UPI002180F047'}, {...}]
                """ OTHER SEQUENCES within this uniref90 ID """
                for protein_record in value['members']:
                    # {'superkingdom': 'Bacteria', 'phylum': 'Proteobacteria', 'class': 'Deltaproteobacteria', 'order': 'Myxococcales',
                    #  'family': 'Nannocystaceae', 'genus': 'Nannocystis', 'species': None}
                    taxon_id = protein_record['organismTaxId']          # 573
                    protein_length = protein_record['sequenceLength']   # 383
                    full_ranks, rank_name_dict = taxon_to_ranks(taxon_id, desired_ranks)

                    if rank_name_dict is None:
                        uniref90_no_query_result_list.append(taxon_id)
                        print(f"The taxon {taxon_id} has no result in the NCBITaxa imported from ete3")
                        continue

                    full_ranks_list = []
                    for rank in desired_ranks:
                        if rank_name_dict[rank]:
                            full_ranks_list.append(rank_name_dict[rank])
                        else:
                            full_ranks_list.append("")

                    taxon_ranks_file.write(key + "," + protein_record["uniref100Id"] + "," + str(protein_length) + "," +
                                           str(taxon_id) + "," + full_ranks + "," + ",".join(full_ranks_list) + "\n")
    taxon_ranks_file.close()
    print(f"In uniref90_to_100_taxon, uniref90_no_query_result_list is {uniref90_no_query_result_list}")
    return uniref_api_response_dict


if __name__ == "__main__":
    pass
