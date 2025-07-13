from collections import defaultdict

import pandas as pd

import GRCh38
import ld
from stories.cCRE import ld as cCRE

# Load mapping between ROIs and transcripts
overlaps = pd.read_pickle(cCRE.cCRE.overlaps.pkl)
overlaps = overlaps.loc[overlaps['roi-type'] == 'PLS', ['Transcript ID', 'imputed', 'seqid', 'roi-start', 'roi-end']]

sequences = pd.read_pickle(cCRE.sequences.saveto)
sequences = sequences[sequences['roi-type'] == 'PLS']

overlaps = overlaps.merge(sequences, on=['seqid', 'roi-start', 'roi-end'], how='outer')
overlaps = overlaps[['Transcript ID', 'imputed', 'seqid', 'roi-norm-start', 'roi-norm-end']]
assert overlaps.isna().sum().sum() == 0, "There are NaN values in the overlaps!"

# Match ROI coordinates with Transcript IDs and imputed status
roi2transcripts, roi2imputed = defaultdict(set), {}
for tid, imputed, seqid, start, end in overlaps.itertuples(index=False, name=None):
    roi2transcripts[seqid, start, end].add(tid)
    if (seqid, start, end) not in roi2imputed:
        roi2imputed[seqid, start, end] = imputed
    assert roi2imputed[seqid, start, end] == imputed, (seqid, start, end, imputed)
roi2transcripts = {k: sorted(v) for k, v in roi2transcripts.items()}

# Load pre-calculated scores for all promoters and match them with Transcript IDs / Imputed status
responses = pd.read_pickle(ld.scoring.raw_signal)
responses['Transcript ID'] = [
    roi2transcripts[k] for k in zip(responses['seqid'], responses['roi-norm-start'], responses['roi-norm-end'])
]
assert responses['Transcript ID'].apply(len).min() > 0, "Some ROIs do not have any transcripts!"

responses['imputed'] = [
    roi2imputed[k] for k in zip(responses['seqid'], responses['roi-norm-start'], responses['roi-norm-end'])
]
assert responses['imputed'].isna().sum() == 0, "There are NaN values in the imputed scores!"

# Reference promoters are those that are not imputed and have direct matches with ENCODE cCREs
responses['is-reference'] = ~responses['imputed']

# Z-score the response scores for each promoter using reference promoters
motifs = [col for col in responses.columns if isinstance(col, tuple) and len(col) == 2]
references = responses.loc[responses['is-reference'], motifs]
mean, std = references.mean(), references.std()

for col in motifs:
    responses[col] = (responses[col] - mean[col]) / std[col]

# Clean up the DataFrame and aggregate to the level of clusters
responses = responses[['Transcript ID', *motifs]].copy()
jaspar = pd.read_pickle(ld.jaspar.parsed_clusters)

# Drop TF names
renaming = {x: x[0] for x in motifs}
assert len(set(renaming.values())) == len(renaming)
responses = responses.rename(columns=renaming)

for cluster, ids in jaspar[['cluster', 'id']].itertuples(index=False, name=None):
    # Check that all motifs in the cluster are present in the responses
    assert all(id in responses.columns for id in ids), \
        f"Cluster {cluster} has missing motifs: {set(ids) - set(responses.columns)}"

    responses[cluster] = responses[list(ids)].sum(axis=1)
responses = responses.drop(columns=list(renaming.values()))

# Aggregate to the level of genes
gencode = GRCh38.gencode.load()
tid2gene = {tid.split('.')[0]: rna.gene.split('.')[0] for tid, rna in gencode.rnas.items()}

responses['Gene ID'] = responses['Transcript ID'].map(lambda x: {tid2gene[tid] for tid in x})
responses = responses.drop(columns=['Transcript ID']).explode('Gene ID').groupby('Gene ID').max()

# Re-standardize the motif cluster scores
for motif in responses.columns:
    responses[motif] = (responses[motif] - responses[motif].mean()) / responses[motif].std()

responses.to_pickle(ld.scoring.response, protocol=-1)
