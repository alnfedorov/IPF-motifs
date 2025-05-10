import pickle
from collections import defaultdict
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

import ld
from stories.cCRE import ld as cCRE


def job(fname: Path):
    with open(fname, 'rb') as f:
        data = pickle.load(f)

    responses = {}
    for region, scores in data.items():
        responses[region] = {}
        for key, (forward, reverse) in scores.items():
            responses[region][key] = max(forward.max(), reverse.max())

    records = defaultdict(list)
    for (seqid, start, end), scores in responses.items():
        records['seqid'].append(seqid)
        records['roi-norm-start'].append(start)
        records['roi-norm-end'].append(end)
        for key, value in scores.items():
            records[key].append(value)
    return dict(records)


# Calculate motif response score for each promoter as a maximum of its forward and reverse complement scores
# Using 4 processes to avoid running out of memory
responses = Parallel(n_jobs=4, verbose=100, pre_dispatch='all', batch_size=1)(
    delayed(job)(fname) for fname in ld.scoring.raw_signal.iterdir()
)
merged = responses.pop()
for record in responses:
    for key, value in record.items():
        merged[key].extend(value)
responses = pd.DataFrame(merged)

# Get normalized ROI coordinates for all promoters
df = pd.read_pickle(cCRE.cCRE.overlaps.pkl)
df = df[df['roi-type'] == 'PLS'].copy()

sequences = pd.read_pickle(cCRE.sequences.saveto).drop(columns=['sequence'])
df = df.merge(sequences, on=['roi-type', 'seqid', 'roi-start', 'roi-end'], how='left')
df = df[['Gene ID', 'roi-type', 'seqid', 'roi-norm-start', 'roi-norm-end']]

# Join the calculate response data
df = df.merge(responses, on=['seqid', 'roi-norm-start', 'roi-norm-end'], how='left')
df = df.drop(columns=['seqid', 'roi-norm-start', 'roi-norm-end'])

# Calculate the final score -> max response across all matched gene promoters
df = df.groupby(['Gene ID', 'roi-type']).max().reset_index()
df.to_pickle(ld.scoring.response_per_motif, protocol=pickle.HIGHEST_PROTOCOL)
