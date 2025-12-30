"""
Affinity propagation algorithm

module add diamond/0.9.18
module add python/cpu/3.6.5

"""

# import networkx as nx
import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation
# from sklearn.metrics.pairwise import pairwise_distances


def affinity_propagation_clustering(pairwise_identity_df):
    """
    Perform Affinity Propagation clustering on sequences.
    pairwise_identity_df = identity_res_df_97
    Args:
        import numpy as np (DataFrame): Input data with columns [query_sequence, target_sequence, identity].
    Returns:
        dict: A dictionary mapping cluster labels to lists of sequences.
    """
    # Create a list of unique sequences
    sequences = pd.unique(pairwise_identity_df[['Seq1', 'Seq2']].values.ravel('K'))
    n = len(sequences)
    # Create a similarity matrix
    sim_matrix = np.zeros((n, n))
    seq_to_index = {seq: i for i, seq in enumerate(sequences)}
    for _, row in pairwise_identity_df.iterrows():
        i, j = seq_to_index[row['Seq1']], seq_to_index[row['Seq2']]
        sim_matrix[i, j] = row['Identity']
        sim_matrix[j, i] = row['Identity']  # Symmetric
    # Affinity Propagation clustering
    ap = AffinityPropagation(affinity='precomputed', random_state=9527)
    labels = ap.fit_predict(sim_matrix)
    # Map clusters to sequences
    clusters = {f'cluster_{label}': [] for label in np.unique(labels)}
    for idx, label in enumerate(labels):
        clusters[f'cluster_{label}'].append(sequences[idx])
    return clusters


"""
# Example usage
data = pd.DataFrame({
    'query_sequence': ['seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'],
    'target_sequence': ['seq2', 'seq3', 'seq4', 'seq5', 'seq1', 'seq3'],
    'identity': [98, 99, 97, 100, 99, 98]
})

clusters_ap = affinity_propagation_clustering(data)
print(clusters_ap)
"""


def greedy_affinity_clustering(pairwise_identity_df, centroid_uniref100_uniref90_dict, iden_threshold=95):
    """
    Self-built greedy clustering base on identity threshold and number of sequences.
    pairwise_identity_df = identity_res_df_97
    :param pairwise_identity_df: pandas DF of three columns ['Seq1', 'Seq2', 'Identity'], two direction
    :param centroid_uniref100_uniref90_dict: {centroid_uniref100: {"UniRef90": UniRef90_ID, "n_seq": 20}}
    :param iden_threshold:
    :return: group assignment dict like {0: [unref100_A, unref100_B], 1: [unref100_C]}
    the first of output list is the representative one
    """
    processed_centroid_id = set()
    # sort centroid_uniref100 according to protein length in decreasing order
    sorted_centroid = sorted(centroid_uniref100_uniref90_dict.keys(), key=lambda k: centroid_uniref100_uniref90_dict[k]["n_seq"], reverse=True)
    pairwise_identity_connected_df = pairwise_identity_df.loc[pairwise_identity_df["Identity"] >= iden_threshold]
    group_index = 0
    centroid_group_dict = {}    # {0: [unref100_A, unref100_B], 1: [unref100_C]}
    for centroid_id in sorted_centroid:
        if centroid_id in processed_centroid_id:
            continue
        # if this centroid_id has not bee recorded
        connected_df = pairwise_identity_connected_df.loc[pairwise_identity_connected_df["Seq1"] == centroid_id]
        if connected_df.shape[0] == 0:
            # it is a separate uniref90, its centroid seq is not within the 95% of other uniref90
            centroid_group_dict.update({group_index: [centroid_id]})
        else:
            # it should be combined with a previous
            centroid_group_dict.update({group_index: [centroid_id] + list(connected_df["Seq2"])})
        processed_centroid_id.update(centroid_group_dict[group_index])
        group_index += 1
    return centroid_group_dict
