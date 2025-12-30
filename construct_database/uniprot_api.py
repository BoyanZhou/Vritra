import requests


"""
1. 给定uniref的ID，提取出其对应的另一个uniref的ID。例如，给定uniref90，提取出其每个id包含的uniref100

# metagenomic requires, then diamond v2.0.15.153 will be loaded automatically
module add python/cpu/3.7.2
# !!! but for this script, we must use python 3.6.5, because only under this python, ete3 is available 
module add python/cpu/3.6.5
"""


# 给定UniRef ID list， 返还full result dict
def uniref_full_dict(input_uniref_id, in_uniref_type):
    """
    :param input_uniref_id: list of uniref id ["UniRef50_P40817" or "P40817", ...]
    :param in_uniref_type: must be 'UniRef100' or 'UniRef90' or 'UniRef50'
    :return: full dict of uniref_query
    """
    if in_uniref_type not in ['UniRef100', 'UniRef90', 'UniRef50']:
        raise ValueError("Invalid option. in_uniref_type must be 'UniRef100' or 'UniRef90' or 'UniRef50'.")

    # 1. add "UniRef50_" to the uniref_id_list, if the id does not contain the prefix
    uniref_id_list = []     # transform to the format of "UniRef90_A0A371S335"
    for id_temp in input_uniref_id:
        if id_temp.startswith("UniRef"):
            if id_temp.startswith(in_uniref_type):
                uniref_id_list.append(id_temp)   # format correct, like "UniRef90_A0A371S335"
            else:
                raise ValueError(f"Invalid protein ID. {id_temp} should start with {in_uniref_type}!")
        else:
            uniref_id_list.append(f"{in_uniref_type}_{id_temp}")      # append prefix, like UniRef90_ + "A0A371S335"
    uniref_id_list = list(set(uniref_id_list))  # in the same format "UniRef90_A0A371S335"

    # 2. full dict result of a uniref id
    uniref_full_dict = {}       # {uniref_id: result_dict or None}
    uniref_api_response_dict = {"Success 200: Success": [], "error 404: Resource not found": [],
                                "error 403: Access forbidden": [], "error 500: Internal server error": [],
                                "unexpected error": []}
    for id_temp in uniref_id_list:
        response_res, response_name = uniref_id_query(id_temp, out_format="json")
        uniref_full_dict.update({id_temp: response_res})
        uniref_api_response_dict[response_name].append(id_temp)
    return uniref_full_dict, uniref_api_response_dict


# 单个query的结果，主要被函数uniref_full_dict调用
def uniref_id_query(uniref_id, out_format):
    """
    :param uniref_id="UniRef90_P40817"
    :param out_format="json" or "tsv", "fasta", ... ...
    :return: single query result of one uniref id
    """
    # Define the URL
    url = f"https://rest.uniprot.org/uniref/{uniref_id}.{out_format}"
    # Send a GET request to the URL
    response = requests.get(url)
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response
        print(f"Success: for {uniref_id}")
        if out_format == "json":
            response_dict = response.json()
            response_name = "Success 200: Success"
            return response_dict, response_name
    else:
        if response.status_code == 404:
            print(f"For {uniref_id}, error 404: Resource not found.")
            response_name = "error 404: Resource not found"
        elif response.status_code == 403:
            print(f"For {uniref_id}, error 403: Access forbidden.")
            response_name = "error 403: Access forbidden"
        elif response.status_code == 500:
            print(f"For {uniref_id}, error 500: Internal server error.")
            response_name = "error 500: Internal server error"
        else:
            print(f"For {uniref_id}, unexpected error: {response.status_code}")
            response_name = "unexpected error"
        return None, response_name


# 得到单个uniref100 的fasta文本序列
def get_seq_by_uniref100(uniref100_id):
    """
    get the text sequence from a single uniref100_id
    :param uniref100_id:
    :return: '>UniRef100_A0A4S5C724 Diaminopropionate ammonia-lyase n=1 Tax=Aeromonas veronii TaxID=654 RepID=A0A4S5C724_AERVE\nMSQFSLKMDIADNRFFTGDPSPLFSR\n'
    """
    # Define the URL
    url = f"https://rest.uniprot.org/uniref/{uniref100_id}.fasta"
    # Send a GET request to the URL
    response = requests.get(url)
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response
        print(f"Success: for {uniref100_id}")
        out_seq = response.text
        response_name = "Success 200: Success"
        return out_seq, response_name
    else:
        if response.status_code == 404:
            print(f"For {uniref100_id}, error 404: Resource not found.")
            response_name = "error 404: Resource not found"
        elif response.status_code == 403:
            print(f"For {uniref100_id}, error 403: Access forbidden.")
            response_name = "error 403: Access forbidden"
        elif response.status_code == 500:
            print(f"For {uniref100_id}, error 500: Internal server error.")
            response_name = "error 500: Internal server error"
        else:
            print(f"For {uniref100_id}, unexpected error: {response.status_code}")
            response_name = "unexpected error"
        return None, response_name


"""
# 根据UniRef100 ID, 下载单个蛋白质序列
url = "https://rest.uniprot.org/uniref/UniRef100_A0A4S5C724.fasta"
response = requests.get(url)
response.text   # >UniRef100_A0A4S5C724 Diaminopropionate ammonia-lyase n=1 Tax=Aeromonas veronii TaxID=654 RepID=A0A4S5C724_AERVE\nMSQFSLKMDIADNRFFTGDPSPLFSR\n
"""

"""
# Define the URL
url = "https://rest.uniprot.org/uniref/UniRef50_P40817.json"

# Send a GET request to the URL
response = requests.get(url)

# Check if the request was successful
if response.status_code == 200:
    # Parse the JSON response
    data = response.json()
    print(data)
    data.keys() # dict_keys(['id', 'name', 'memberCount', 'updated', 'entryType', 'commonTaxon', 'seedId', 'goTerms', 'representativeMember', 'members'])
    data["id"]  # 'UniRef50_P40817'
    data["name"]    # 'Cluster: Diaminopropionate ammonia-lyase'
    data["entryType"]   # 'UniRef50'
    data["commonTaxon"] # {'scientificName': 'root', 'taxonId': 1}
    data["goTerms"]     # [{'goId': 'GO:0016829', 'aspect': 'GO Molecular Function'}, {'goId': 'GO:0016020', 'aspect': 'GO Cellular Component'}]
    data['representativeMember'] # {'memberIdType': 'UniProtKB ID',
 'memberId': 'DPAL_SALTY',
 'organismName': 'Salmonella typhimurium (strain LT2 / SGSC1412 / ATCC 700720)',
 'organismTaxId': 99287,
 'sequenceLength': 404,
 'proteinName': 'Diaminopropionate ammonia-lyase',
 'accessions': ['P40817'],
 'uniref90Id': 'UniRef90_P40817',
 'uniref100Id': 'UniRef100_P40817',
 'uniparcId': 'UPI0000129737',
 'sequence': {'value': 'MHELIKYQFNTRRKKYGTGAALSLLNGNVGHEVLAFHKKLPNYAVTPLHNLAHLSQRLGLGSIHIKDESWRFGLNAFKGLGGSYAVGKYLADKLQCDINSLSFAALNTPEIKEKIKDCVFVTATDGNHGRGVAWAAEQLGLKAVVYMPKGSSLIRAENIRHHGAECTITDLNYDDAVRLAHRMAQTKGWVLLQDTAWTGYEEIPTWIMQGYMTLAVEAYEQLAETNSPLPTHLILQAGVGSFAGSVMGYFVEKMQENIPNIIVVEPHQANCLYQSAVMDDGQPHCVTGDMATIMAGLACGEPNIISWPIIRDNTSCFISADDCLAAKGMRISAAPRPGTDTPFISGESGAIGVGLLYELMNNMHYQDLANRLQLDASAHVLLISTEGDTSPDIYEDIVWNGRSA',
  'length': 404,
  'molWeight': 44152,
  'crc64': '9D605DFEAA691F00',
  'md5': 'F31BA5E5BCDB09A5E08C4780C8CBF949'}}
    
    
    data['members']     # [{'memberIdType': 'UniParc',
  'memberId': 'UPI002180F047',
  'organismName': 'Klebsiella pneumoniae',
  'organismTaxId': 573,
  'sequenceLength': 383,
  'proteinName': 'diaminopropionate ammonia-lyase',
  'uniref90Id': 'UniRef90_W9B689',
  'uniref100Id': 'UniRef100_UPI002180F047'},
 {'memberIdType': 'UniParc',
  'memberId': 'UPI00129D25AA',
  'organismName': 'Klebsiella pneumoniae',
  'organismTaxId': 573,
  'sequenceLength': 383,
  'proteinName': 'diaminopropionate ammonia-lyase',
  'uniref90Id': 'UniRef90_W9B689',
  'uniref100Id': 'UniRef100_UPI00129D25AA'}]

    # Optionally save the JSON data to a file
    with open("UniRef90_P40817.json", "w") as json_file:
        json.dump(data, json_file, indent=4)

    print("Data downloaded and saved successfully!")
else:
    print(f"Failed to retrieve data. Status code: {response.status_code}")
"""

