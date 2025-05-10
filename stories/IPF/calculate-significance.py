import pickle

import pandas as pd
from scipy import stats

import ld
from stories.motifs import ld as motifs

motifs = pd.read_pickle(motifs.scoring.response_per_cluster)
assert (motifs['roi-type'] == 'PLS').all(), motifs['roi-type'].unique_values()
motifs = motifs.drop(columns=['roi-type']).set_index('Gene ID').copy()

# Mark motifs with Z-score above the threshold as True ('positive') and the rest as False
for m in motifs.columns:
    vals = motifs[m].values
    assert abs(vals.mean()) <= 1e-3
    assert abs(vals.std() - 1) <= 1e-3

    vals = (vals - vals.mean()) / vals.std()
    vals = vals >= ld.thresholds.zscore
    motifs[m] = vals

# Load the genes of interest
genes = pd.read_pickle(ld.SCRNA_FOLD_CHANGE)

df = []
for (tissue, subset), fc in genes.items():
    # Filter genes to only include those with matched scRNA data
    before = len(fc)
    data = pd.merge(motifs, fc, how='inner', left_index=True, right_index=True)
    print(f"{subset}: {before} -> {len(data)}({len(data) / before:.2%})")

    records = []
    for m in motifs.columns:
        hasmotif = data.loc[motifs[m], 'log2FoldChange']
        nomotif = data.loc[~motifs[m], 'log2FoldChange']

        if len(hasmotif) == 0 or len(nomotif) == 0:
            continue

        # Wilcoxon rank-sum (Mann Whitney U) test between the two groups
        pvalue = stats.mannwhitneyu(hasmotif, nomotif, alternative='two-sided').pvalue

        records.append({
            'motif': m, 'p-value': pvalue,
            'Median log2 fold change': hasmotif.median(),
            'Control log2 fold change': nomotif.median(),
            'Have motif': len(hasmotif), 'No motif': len(nomotif)
        })
    records = pd.DataFrame(records)
    records['tissue'] = tissue
    records['subset'] = subset
    df.append(records)

df = pd.concat(df)

ld.STAT_TESTS.parent.mkdir(parents=True, exist_ok=True)
df.to_pickle(ld.STAT_TESTS, protocol=pickle.HIGHEST_PROTOCOL)
