import os
import sys


def protein_id_from_fas(fas_path, fas_type):
    """
    :param fas_path:
    :param fas_type: fas_type must be 'InterPro' or 'UniProt'.
    :return: protein_id_list like, ["A0A371S335", "A0A089PYJ7"]
    """
    if fas_type not in ['InterPro', 'UniProt']:
        raise ValueError("Invalid option. fas_type must be 'InterPro' or 'UniProt'.")

    # 1. get the headers of all sequences, startwith(">")

    protein_id_list = []    # like, ["A0A371S335", "A0A089PYJ7"]

    with open(fas_path, "r") as fas_f:
        for line in fas_f:
            if line.startswith(">"):
                if fas_type == "InterPro":
                    protein_id = line.split("|")[0][1:].strip()
                else:
                    protein_id = line.split(" ")[0]     # >UniRef90_A0A371S335
                    if protein_id.startswith(">UniRef"):
                        protein_id = protein_id.split("_")[-1]
                    else:
                        continue
                protein_id_list.append(protein_id)
    return protein_id_list


if __name__ == "__main__":
    pass


