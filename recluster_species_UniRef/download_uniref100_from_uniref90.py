# import seqid_from_fas as sff

import construct_database.uniprot_api as u_api
import os
# import sys
import json
import sys
import construct_database.uniref90_to_uniref100_taxon_ranks as uniref_taxon


#################
# main function #
#################
def uniref90_to_100_fas(in_uniref90_fas_path, output_json_path, output_taxon_ranks_path, out_uniref100_fas_path):
    # get the json file and taxon ranks file
    uniref90_api_response_dict = uniref_taxon.uniref90_to_100_taxon(in_uniref90_fas_path, output_json_path, output_taxon_ranks_path, fas_type='UniProt')
    # download UniRef100 sequence
    uniref100_api_response_dict = download_seq_from_json(output_json_path, out_uniref100_fas_path)
    return uniref90_api_response_dict, uniref100_api_response_dict


def download_seq_from_json(json_path, download_seq_path):
    """
    This json is the output of script1_uniref90_fas_to_uniref100_fas_taxon_ranks.py
    :param json_path: dict of UniRef90 ID containing UniRef100
    :param download_seq_path: output the downloaded sequences
    :return:
    """
    # READ JSON
    with open(json_path, 'r') as file:
        uniref_data = json.load(file)
    # OUTPUT FILE
    output_seq_f = open(download_seq_path, 'w')

    # Record the api response for list of UniRef100s
    uniref100_api_response_dict = {"Success 200: Success": [], "error 404: Resource not found": [],
                                   "error 403: Access forbidden": [], "error 500: Internal server error": [],
                                   "unexpected error": []}

    for uniref90_id, uniref90_dict in uniref_data.items():
        # 1. Get the main sequence, representative sequence
        # uniref90_id="UniRef90_UPI001BA6BB4E"
        if uniref90_dict is None:
            print(f"{uniref90_id} has None for uniref90_dict. Skip it!")
            continue
        member_count = uniref90_dict['memberCount']
        uniref100_id = uniref90_dict['representativeMember']['uniref100Id']
        protein_name = uniref90_dict['representativeMember']['proteinName']
        taxon_id = uniref90_dict['representativeMember']['organismTaxId']
        taxon_name = uniref90_dict['representativeMember']['organismName']
        seq_text = uniref90_dict['representativeMember']['sequence']['value']
        seq_header = f">{uniref100_id} {protein_name} n={member_count} Tax={taxon_name} TaxID={taxon_id} RepID={uniref90_id.split('_')[-1]}"
        # >UniRef100_A0A4S5C724 Diaminopropionate ammonia-lyase n=1 Tax=Aeromonas veronii TaxID=654 RepID=A0A4S5C724_AERVE\nMSQFSLKMDIADNRFFTGDPSPLFSR\n
        output_seq_f.write(seq_header + "\n" + seq_text + "\n")
        uniref100_api_response_dict["Success 200: Success"].append(uniref100_id)

        # 2. Check whether there are more 'members', a list containing other sequences
        if 'members' in uniref90_dict:
            current_uniref100_list = [uniref100_id]  # avoid duplicate uniref100 id
            for member_dict in uniref90_dict['members']:
                # {'memberIdType': 'UniParc', 'memberId': 'UPI002602F8B1', 'organismName': 'Desulfovibrio sp.', 'organismTaxId': 885,
                # 'sequenceLength': 405, 'proteinName': 'diaminopropionate ammonia-lyase', 'uniref50Id': 'UniRef50_A0A743TVN6',
                # 'uniref100Id': 'UniRef100_UPI002602F8B1'}
                uniref100_id_temp = member_dict['uniref100Id']
                if uniref100_id_temp in current_uniref100_list:
                    continue
                else:
                    current_uniref100_list.append(uniref100_id_temp)
                    uniref100_seq, response_name = u_api.get_seq_by_uniref100(uniref100_id_temp)
                    uniref100_api_response_dict[response_name].append(uniref100_id_temp)
                    if uniref100_seq is not None:
                        output_seq_f.write(uniref100_seq)

    output_seq_f.close()
    return uniref100_api_response_dict


if __name__ == "__main__":
    pass
