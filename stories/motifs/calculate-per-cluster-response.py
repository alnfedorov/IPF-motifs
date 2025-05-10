import pickle

import numpy as np
import pandas as pd

import ld

jaspar = pd.read_pickle(ld.jaspar.parsed_clusters)
responses = pd.read_pickle(ld.scoring.response_per_motif).set_index(['Gene ID', 'roi-type'])

# # Keep only protein coding genes PLS
# annotation = GRCh38.gencode.load()
# genes = {gene.ind.split('.')[0]: gene for gene in annotation.genes.values()}

# responses['biotype'] = [genes[x].attrs.type for x in responses['Gene ID']]
# mask = responses['biotype'] == 'protein_coding'
# responses = responses[mask].copy().drop(columns=['biotype'])

# Drop TF names
renaming = {x: x[0] for x in responses.columns}
assert len(set(renaming.values())) == len(renaming)
responses = responses.rename(columns=renaming).copy()

# Calculate the score per cluster by summing the z-scores of all motifs in the cluster
cluster_responses = {}
for cluster, ids in jaspar[['cluster', 'id']].itertuples(index=False, name=None):
    # Scores for each gene
    scores = np.zeros(len(responses), dtype=np.float32)

    # Sum Z-normalized scores for each motif in the cluster
    for id in ids:
        scores += (responses[id] - responses[id].mean()) / responses[id].std()
    # Re-standardize the scores
    scores = (scores - scores.mean()) / scores.std()

    assert cluster not in cluster_responses
    cluster_responses[cluster] = scores

cluster_responses = pd.DataFrame(cluster_responses).reset_index()
cluster_responses.to_pickle(ld.scoring.response_per_cluster, protocol=pickle.HIGHEST_PROTOCOL)
