"""
Using Community Detection Method, replacing previous iterative search
Given pre-labeled sequences and other directly/indirectly connected sequences

module add diamond/0.9.18
module add python/cpu/3.6.5

"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import re


# suitable for UniRef90 and UniProtK
def get_seq_id_function_from_fas(fas_path, db_name="uniprot"):
    """
    :param fas_path:
    :param db_name: "uniprot" or "interpro"
    :return:
    """
    # fas_path = "uniref_name_Bile_acid_coenzyme_A_lig_Search_C_2025_02_04.fasta"
    # fas_path = "baiB_labeled_Search_C_connected_and_its_connected_seq.fas"
    seq_header_list = []
    with open(fas_path, "r") as fas_f:
        for line in fas_f:
            if line.startswith(">"):
                seq_header_list.append(line)
    # build seq_id-function dict: from '>UniRef90_A0A117EBN7 Bile acid-coenzyme A ligase n=24 Tax=Streptomyces TaxID=1883 RepID=A0A117EBN7_STRSC\n'
    # OR ">tr|A0A021XAR7|A0A021XAR7_9HYPH Xanthine dehydrogenase molybdenum-binding subunit OS=Shinella sp. DD12 OX=1410620 GN=xdhA PE=3 SV=1"
    seq_id_func_dict = {}
    for seq_info in seq_header_list:
        # >UniRef90_A0A117EBN7 Bile acid-coenzyme A ligase
        seq_id, seq_func = extract_description(seq_info)
        if seq_func and seq_id:
            seq_id_func_dict.update({seq_id: seq_func})
    func_count_dict = Counter(seq_id_func_dict.values())
    return seq_id_func_dict, func_count_dict


# extract the function of protein
def extract_description(text):
    """
    can deal with both uniref90 and UniProtKB
    '>UniRef90_A0A117EBN7 Bile acid-coenzyme A ligase n=24 Tax=Streptomyces TaxID=1883 RepID=A0A117EBN7_STRSC\n'
    ">tr|A0A021XAR7|A0A021XAR7_9HYPH Xanthine dehydrogenase molybdenum-binding subunit OS=Shinella sp. DD12 OX=1410620 GN=xdhA PE=3 SV=1"
    :param text:
    :return:
    """
    match = re.match(r'^>(\S+)\s+(.*?)(?=\s+\S+=)', text)
    if match:
        record_id, description = match.groups()
        return record_id, description.strip()
    return None, None


# ----------------------------------------------------------------------------------------------------------------------
# 下面这个是v4版本 (0721之后的版本)， v3见以前的文件
def label_propagation_given_edges_v4(labeled_to_initial_blast_res_df, t1_to_initial_blast_res_df, pre_labeled_id_set):
    """
    :param labeled_to_initial_blast_res_df:
    :param t1_to_initial_blast_res_df:
    :param pre_labeled_id_set: seq id set of pre-labeled seq
    :return: {seq_id: {"to_labeled": {"seq_id": "", "highest_identity": -1}, "to_unlabeled": {"seq_id": "", "highest_identity": -1}}}
    Returned "closer to labeled" and "closer to unlabeled" all from t1. t2 is not included here
    """
    # 1. remove the self-to-self sequences, t1 is contained in original initial
    labeled_to_initial_df = labeled_to_initial_blast_res_df.loc[labeled_to_initial_blast_res_df.iloc[:, 0] != labeled_to_initial_blast_res_df.iloc[:, 1]]
    labeled_to_initial_df.columns = ['Query', 'Target', 'Identity'] + labeled_to_initial_df.columns[3:].tolist()
    t1_to_initial_df = t1_to_initial_blast_res_df.loc[t1_to_initial_blast_res_df.iloc[:, 0] != t1_to_initial_blast_res_df.iloc[:, 1]]
    t1_to_initial_df.columns = ['Query', 'Target', 'Identity'] + t1_to_initial_df.columns[3:].tolist()

    # 2. divide t1_to_initial_df into two df: a. all pre_labeled; b. unlabeled
    t1_to_initial_labeled_df = t1_to_initial_df.loc[t1_to_initial_df['Target'].isin(pre_labeled_id_set)]
    t1_to_initial_unlabeled_df = t1_to_initial_df.loc[~t1_to_initial_df['Target'].isin(pre_labeled_id_set)]

    #
    t1_seq_id_set = set(t1_to_initial_df['Query'])

    # 对每个t1进行循环
    closer_labeled_dict = {}
    closer_unlabeled_dict = {}
    for seq_id_temp in t1_seq_id_set:
        # seq_id_temp = 'UniRef90_A0A1A2TS18'
        ###################################################################
        # t1 to labeled seq取最大的identity: stored in t1_to_labeled_dict #
        ###################################################################
        t1_to_labeled_dict = {"seq_id": "", "highest_identity": -1}
        edge_df_t1_current0 = labeled_to_initial_df.loc[labeled_to_initial_df['Target'] == seq_id_temp]
        edge_df_t1_current1 = t1_to_initial_labeled_df.loc[t1_to_initial_labeled_df['Query'] == seq_id_temp]
        if edge_df_t1_current0.shape[0] > 0:
            iden_max_temp = edge_df_t1_current0["Identity"].max()
            if iden_max_temp > t1_to_labeled_dict["highest_identity"]:
                t1_to_labeled_dict["seq_id"] = edge_df_t1_current0.loc[edge_df_t1_current0["Identity"].idxmax(), 'Query']
                t1_to_labeled_dict["highest_identity"] = iden_max_temp
        if edge_df_t1_current1.shape[0] > 0:
            iden_max_temp = edge_df_t1_current1["Identity"].max()
            if iden_max_temp > t1_to_labeled_dict["highest_identity"]:
                t1_to_labeled_dict["seq_id"] = edge_df_t1_current1.loc[edge_df_t1_current1["Identity"].idxmax(), 'Target']
                t1_to_labeled_dict["highest_identity"] = iden_max_temp
        ####################################
        # t1 to unlabeled 取最大的identity #
        ####################################
        edge_df_t1_unlabeled_current = t1_to_initial_unlabeled_df.loc[t1_to_initial_unlabeled_df["Query"] == seq_id_temp]
        t1_to_unlabeled_dict = {"seq_id": "", "highest_identity": -1}
        if edge_df_t1_unlabeled_current.shape[0] > 0:
            iden_max_temp = edge_df_t1_unlabeled_current["Identity"].max()
            if iden_max_temp > t1_to_unlabeled_dict["highest_identity"]:
                t1_to_unlabeled_dict["seq_id"] = edge_df_t1_unlabeled_current.loc[edge_df_t1_unlabeled_current["Identity"].idxmax(), 'Target']
                t1_to_unlabeled_dict["highest_identity"] = iden_max_temp

        #################################################
        # Current Seq Is Closer To Labeled Or Unlabeled #
        #################################################
        if t1_to_labeled_dict["highest_identity"] > t1_to_unlabeled_dict["highest_identity"]:
            closer_labeled_dict.update({seq_id_temp: {"to_labeled": t1_to_labeled_dict, "to_unlabeled": t1_to_unlabeled_dict}})
        elif t1_to_unlabeled_dict["highest_identity"] > t1_to_labeled_dict["highest_identity"]:
            closer_unlabeled_dict.update({seq_id_temp: {"to_labeled": t1_to_labeled_dict, "to_unlabeled": t1_to_unlabeled_dict}})
        else:
            # both of them are -1
            print(f"For {seq_id_temp}, no connection is detected in labeled_to_initial_blast_res_df, t1_to_initial_blast_res_df.")
    print(f"In detect_commnunity.py, {len(closer_labeled_dict)} seqs closer to labeled; {len(closer_unlabeled_dict)} seqs closer to unlabeled.")
    # here "closer to labeled" and "closer to unlabeled" all from t1. t2 is not included here
    return closer_labeled_dict, closer_unlabeled_dict


# 用来summarize新标记的seq
def sum_new_labeled_func(all_seq_func_dict, pre_labeled_seq, quasi_labeled_seq, new_labeled_seq, new_unlabeled_seq, initial_blast_remain_seq_set):
    # Get the count of each unique function in three set
    pre_labeled_func_count = Counter([all_seq_func_dict[i] for i in pre_labeled_seq])
    quasi_labeled_func_count = Counter([all_seq_func_dict[i] for i in quasi_labeled_seq])   # is a subset of pre_labeled
    new_labeled_func_count = Counter([all_seq_func_dict[i] for i in new_labeled_seq])
    new_unlabeled_func_count = Counter([all_seq_func_dict[i] for i in new_unlabeled_seq])
    initial_blast_remain_func_count = Counter([all_seq_func_dict[i] for i in initial_blast_remain_seq_set])

    # Get all unique keys
    all_func_set = set(pre_labeled_func_count.keys()) | set(new_labeled_func_count.keys()) | set(new_unlabeled_func_count.keys()) | set(initial_blast_remain_func_count.keys())

    # Construct DataFrame
    summary_df = pd.DataFrame({
        "Sequence": list(all_func_set),
        "All_Labeled_Number": [pre_labeled_func_count.get(k, 0) + new_labeled_func_count.get(k, 0) for k in all_func_set],
        "Pre_Labeled_Number": [pre_labeled_func_count.get(k, 0) for k in all_func_set],
        "Quasi_Labeled_Number": [quasi_labeled_func_count.get(k, 0) for k in all_func_set],
        "New_Labeled_Number": [new_labeled_func_count.get(k, 0) for k in all_func_set],
        "New_Unlabeled_Number": [new_unlabeled_func_count.get(k, 0) for k in all_func_set],
        "Raw_Dataset_Remained_Number": [initial_blast_remain_func_count.get(k, 0) for k in all_func_set]
    })

    # Sort by the second column in descending order
    summary_df = summary_df.sort_values(by=["Pre_Labeled_Number", "Quasi_Labeled_Number", "New_Labeled_Number"], ascending=[False, False, False]).reset_index(drop=True)
    return summary_df


#################
# main function #
#################
def detect_new_seq_in_community(pre_labeled_seq_fas_path, quasi_labeled_seq_fas_path,
                                pre_labeled_connected_and_indirect_fas_path,
                                initial_blast_seq_from_seed_path,
                                labeled_to_initial_blast_res_path, t1_to_initial_blast_res_path
                                # total_self_blast_res_path
                                ):
    """
    :param pre_labeled_seq_fas_path= "xdhA_refined_labeled_filtered_seqs.fas",
    uniref90 + uniprotkb & identity-to-initial > 35%
    :param quasi_labeled_seq_fas_path= "xdhA_quasi_labeled_filtered_seqs.fas"
    :param pre_labeled_connected_and_indirect_fas_path= "xdhA_labeled_connected_and_its_connected_seq.fas", t1 & t2
    t1 + t2 without pre-labeled and quasi-labeled
    :param initial_blast_seq_from_seed_path= "xdhA_initial_sequence_set.fas"
    :param total_self_blast_res_path= "self_blasted_res_alignments200.txt", include both pre-labeled, t1, and t2
    :return: closer_labeled_dict is like {}
    """
    # total_self_blast_res_df = pd.read_table(total_self_blast_res_path, header=None)
    # total_self_blast_res_df_iden_70 = total_self_blast_res_df.loc[total_self_blast_res_df.iloc[:, 2] >= iden_threshold]

    labeled_to_initial_blast_res_df = pd.read_table(labeled_to_initial_blast_res_path, header=None)
    t1_to_initial_blast_res_df = pd.read_table(t1_to_initial_blast_res_path, header=None)

    # 1. previously labeled uniref90 + uniprotkb (then, add quasi labeled to pre-labeled)
    pre_labeled_seq_id_func_dict, pre_labeled_func_count_dict = get_seq_id_function_from_fas(pre_labeled_seq_fas_path)
    quasi_labeled_seq_id_func_dict, quasi_labeled_func_count_dict = get_seq_id_function_from_fas(quasi_labeled_seq_fas_path)
    pre_labeled_seq_id_func_dict.update(quasi_labeled_seq_id_func_dict)
    # 2. previously unlabeled seqs (i.e. t1 + t2)
    other_seq_id_func_dict, other_func_count_dict = get_seq_id_function_from_fas(pre_labeled_connected_and_indirect_fas_path)
    # 3. the full uniref90 set, i.e. initial blast from seed); then remove seq id above
    initial_blast_seq_id_func_dict, initial_blast_func_count_dict = get_seq_id_function_from_fas(initial_blast_seq_from_seed_path)
    initial_remain_id_set = set(initial_blast_seq_id_func_dict.keys()) - set(pre_labeled_seq_id_func_dict.keys()) - set(other_seq_id_func_dict.keys())

    # 所有序列被注释的功能，合并上面三个文件中所有的seq id with corresponding function
    all_seq_id_func_dict = pre_labeled_seq_id_func_dict.copy()
    all_seq_id_func_dict.update(other_seq_id_func_dict)
    all_seq_id_func_dict.update({i: initial_blast_seq_id_func_dict[i] for i in initial_remain_id_set})
    # 检测直接连接的序列是否与pre-labeled更加接近， 类似Label Propagation Algorithm
    closer_labeled_dict, closer_unlabeled_dict = label_propagation_given_edges_v4(labeled_to_initial_blast_res_df,
                                                                                  t1_to_initial_blast_res_df, set(
                                                                                   pre_labeled_seq_id_func_dict.keys()))

    final_summary_table = sum_new_labeled_func(all_seq_id_func_dict, pre_labeled_seq_id_func_dict.keys(), quasi_labeled_seq_id_func_dict.keys(),
                                               closer_labeled_dict.keys(), closer_unlabeled_dict.keys(),
                                               initial_remain_id_set)
    return final_summary_table, closer_labeled_dict, closer_unlabeled_dict


"""
测试代码如下
Pre_labeled_seq_fas_path= "uniref_name_Bile_acid_coenzyme_A_lig_Search_C_2025_02_04.fasta"
Pre_labeled_connected_and_indirect_fas_path= "baiB_labeled_Search_C_connected_and_its_connected_seq.fas"
Total_self_blast_res_path= "self_blasted_res_alignments200.txt"
Final_summary_table = detect_new_seq_in_community(Pre_labeled_seq_fas_path, Pre_labeled_connected_and_indirect_fas_path, Total_self_blast_res_path, iden_threshold=70)

print(Final_summary_table.to_string())



# 这个函数暂时用不着
def visulize_labeled_seq(edge_df, pre_labeled):
    # pre_labeled = set(labeled_seq_id_func_dict.keys())
    # edge_df = total_self_blast_res_df_iden_70
    edge_df = edge_df.loc[edge_df.iloc[:, 0] != edge_df.iloc[:, 1]]
    edge_df_labeled_only = edge_df.loc[edge_df.iloc[:, 0].isin(pre_labeled) & edge_df.iloc[:, 1].isin(pre_labeled)]
    node1_col, node2_col, weight_col = edge_df_labeled_only.columns[:3]
    G = nx.from_pandas_edgelist(edge_df_labeled_only, source=node1_col, target=node2_col, edge_attr=[weight_col], create_using = nx.Graph())
    len(G.nodes())
    # len(set(edge_df_labeled_only.iloc[:, 0]))
    # len(pre_labeled & set(edge_df.iloc[:, 0]))
    # len(pre_labeled & set(edge_df.iloc[:, 1]))
    plt.figure(figsize=(10, 7))
    pos = nx.spring_layout(G, seed=42)  # Layout for visualization
    nx.draw(G, pos, with_labels=True, node_size=500,
            edge_color="gray", font_size=8)
    plt.title("Label nodes only")
    plt.show()
"""

