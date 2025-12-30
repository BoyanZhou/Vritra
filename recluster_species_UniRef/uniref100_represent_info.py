"""
Build a class, contain all information represented by the final chosen UniRef100
"""

from collections import Counter


class UniRepInfo:
    def __init__(self, rep_uniref100, represented_uniref100_list, uniref100_90_dict, uniref100_species_dict):
        self.rep_uniref100 = rep_uniref100
        self.represented_uniref100_list = represented_uniref100_list
        self.corresponding_uniref90_list = [uniref100_90_dict[i] for i in represented_uniref100_list]
        # ['Micromonospora rifamycinica', 'Micromonospora sp. WMMD714', 'Micromonospora rifamycinica']
        self.corresponding_species_list = [uniref100_species_dict[i] for i in represented_uniref100_list]
        species_count_list = [f"{species_name}:{count_n}" for species_name, count_n in Counter(self.corresponding_species_list).items()]
        # "Micromonospora rifamycinica:2;Micromonospora sp. WMMD714:1"
        self.corresponding_species_count_summary = ";".join(species_count_list)
        self.rep_species = ""

    def get_represented_species(self, species_ranks_count_dict):
        """
        :param species_ranks_count_dict: for each species name, the number of slots != "" in seven taxon ranks
        :return:
        """
        # species_count_dict = Counter(['Micromonospora rifamycinica', 'Micromonospora sp. WMMD714', 'Micromonospora rifamycinica'])
        species_count_dict = Counter(self.corresponding_species_list)
        species_name_list = list(species_count_dict.keys())
        # 1. First choose the species names with taxon ranks as full as possible (i.e. full info in taxon ranks)
        n_ranks_max = max(species_ranks_count_dict[i] for i in species_name_list)
        # Find all species that have the max number of ranks
        species_name_list1 = [i for i in species_name_list if species_ranks_count_dict[i] == n_ranks_max]
        if len(species_name_list1) == 1:
            self.rep_species = species_name_list1[0]

        # 2. choose the species not contain "sp.", "bacterium", "uncultured"
        # species_name_list1 = ['Micromonospora rifamycinica', 'Micromonospora sp. WMMD714', 'Micromonospora bacterium']
        species_name_list2 = [e for e in species_name_list1 if
                              not any(substring in e.lower() for substring in ["sp.", "bacterium", "uncultured"])]
        # 3. choose the species with max count
        if len(species_name_list2) == 0:
            # if no prior species name need to be considered, use all rest of species in species_name_list1
            species_name_list2 = species_name_list1
        species_name_count_dict = Counter(species_name_list2)
        # species_name_list3 = ["sp.", "bacterium", "uncultured"]
        species_name_list3 = [i for i, j in species_name_count_dict.items() if j == max(species_name_count_dict.values())]
        species_name_list3_sorted = sorted(species_name_list3)
        self.rep_species = species_name_list3_sorted[0]
