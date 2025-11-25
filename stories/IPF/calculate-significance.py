import pickle

import pandas as pd
from scipy import stats

import ld
from stories.motifs import ld as motifs

# Mark motifs with Z-score above the threshold as True ('positive') and below as False ('negative')
motifs = pd.read_pickle(motifs.scoring.response)
motif_groups = {}
for m in motifs.columns:
    motif_groups[m] = {
        "positive": set(motifs.index[motifs[m] >= ld.thresholds.zscore]),
        "negative": set(motifs.index[motifs[m] < ld.thresholds.zscore])
    }

# Load the genes of interest
genes = pd.read_pickle(ld.SCRNA_FOLD_CHANGE)

df = []
for (tissue, subset), fc in genes.items():
    records = []
    for motif, gids in motif_groups.items():
        hasmotif = fc[fc.index.isin(gids['positive'])]
        nomotif = fc[fc.index.isin(gids['negative'])]

        if len(hasmotif) == 0 or len(nomotif) == 0:
            continue

        # Wilcoxon rank-sum (Mann Whitney U) test between the two groups
        pvalue = stats.mannwhitneyu(hasmotif, nomotif, alternative='two-sided').pvalue

        records.append({
            'motif': motif, 'p-value': pvalue,
            'Δ(Median log2FoldChange)': hasmotif.median() - nomotif.median(),
            'Have motif': len(hasmotif), 'No motif': len(nomotif), 'Total genes': len(fc)
        })
    records = pd.DataFrame(records)
    records['tissue'] = tissue
    records['subset'] = subset
    df.append(records)

df = pd.concat(df)

# Calculate q-values (FDR) for each tissue
records = []
for tissue, subdf in df.groupby(['tissue']):
    subdf['q-value'] = stats.false_discovery_control(subdf['p-value'], method='bh')
    records.append(subdf)
df = pd.concat(records).sort_values(['tissue', 'subset', 'q-value'])

# Save all results
ld.STAT_TESTS.parent.mkdir(parents=True, exist_ok=True)
df.to_pickle(ld.STAT_TESTS, protocol=pickle.HIGHEST_PROTOCOL)

# Filter results by q-value threshold and save as CSV
df = df[df['q-value'] <= ld.thresholds.qvalue].copy()
df.to_csv(ld.STAT_TESTS.with_suffix('.csv'), index=False)
