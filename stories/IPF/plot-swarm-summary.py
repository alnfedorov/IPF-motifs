import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import ld

plt.rcParams['svg.fonttype'] = 'none'
ld.RESULTS.mkdir(parents=True, exist_ok=True)

PALETTE = {
    'Other': '#B3B3B3',
    'ISRE-like': '#D73027',
    'GAS-like': '#1A9850',
}
REMAP_CLUSTERS = {
    'cluster_062': 'RUNX-like',
    'cluster_005': 'RFX-like',
    'cluster_040': 'Yy1/ZFP42',
    'cluster_037': 'ZBTB12/ZNF524'
}
ORDER = [
    "cDC1", "pDC", "cDC2", "cMono", "ncMono", "SPP1+ Mac", "ISGhi AM", "AM_UD",
    "Intm AM", "FABP4hi AM", "IGF-1+ AM", "MT_AM", "HSP+ Mac", "IM", "CXCL10+ Mono-Mac",
    "MonoDC", "CD206hi FN1hi AM", "MonoMac"
]

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

Y = '-log10(p-value)'

df = pd.read_pickle(ld.STAT_TESTS)
df['-log10(p-value)'] = -df['p-value'].apply(lambda x: np.log10(x))

for tissue, subdf in df.groupby('tissue'):
    subdf = subdf.sort_values(by=['subset'])
    subdf['motif'] = subdf['motif'].replace(REMAP_CLUSTERS)
    print(subdf[subdf['q-value'] <= ld.thresholds.qvalue][[
        'tissue', 'subset', 'motif', 'p-value', 'q-value', 'Δ(Median log2FoldChange)', 'Have motif', 'No motif',
        'Total genes'
    ]])

    subdf['Category'] = 'Other'
    for motif in ['GAS-like', 'ISRE-like']:
        mask = subdf['motif'] == motif
        subdf.loc[mask, 'Category'] = motif

    # subdf['Median effect'] = 'None/Low'
    # subdf.loc[subdf['Median log2 fold change'] < -0.585, 'Median effect'] = 'Downregulated in IPF'
    # subdf.loc[subdf['Median log2 fold change'] > 0.585, 'Median effect'] = 'Upregulated in IPF'

    total_genes = {}
    for subset, subsubdf in subdf.groupby('subset'):
        assert subsubdf['Total genes'].nunique() == 1
        total_genes[subset] = subsubdf['Total genes'].max()
    order = sorted(
        subdf['subset'].unique(), key=lambda x: (ORDER.index(x), -1) if x in ORDER else (-1, x)
    )
    print(order)
    assert set(order) == set(subdf['subset']), set(subdf['subset']) - set(order)

    subdf['X'] = subdf['subset'].map({subset: i for i, subset in enumerate(order)})
    subdf['X'] += np.random.uniform(-0.2, 0.2, size=len(subdf))

    # subdf['Size'] = subdf['Have motif']
    # subdf['Size'] = subdf['Have motif'] / (subdf['Have motif'] + subdf['No motif'])
    # subdf['Size'] = (subdf['Size'] / subdf['Size'].max() * 5).clip(2, 5)

    width = 0.6 * len(order) + 0.5
    fig, ax = plt.subplots(figsize=(width, 6))
    sns.scatterplot(
        data=subdf, y=Y, x='X', ax=ax, hue='Category', s=90, edgecolor='white', linewidth=0.75, palette=PALETTE,
    )

    labels = [f"{subset}\n(n={total_genes[subset]})" for subset in order]
    ax.set_xticks(np.arange(len(order)), labels, rotation=90)
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel(None)
    ax.grid(axis='x')

    # Label ISRE/GAS-like motifs and all significant motifs
    for x, subset in enumerate(order):
        signif = subdf[(subdf['subset'] == subset) & (subdf['q-value'] <= ld.thresholds.qvalue)]
        for motif, y in signif[['motif', Y]].itertuples(index=False):
            fontsize = 6 if motif.startswith("cluster") else 12
            ax.text(
                x, y, motif, ha='left', va='center', color='black', fontweight='bold', fontsize=fontsize
            )

    # FDR line
    minpval = subdf[subdf['q-value'] <= ld.thresholds.qvalue][Y].min()
    ax.axhline(minpval, color='black', linestyle='--', linewidth=1)
    ax.text(
        1.0, minpval, f'FDR ≤ {ld.thresholds.qvalue}', ha='right', va='bottom', color='black',
        transform=ax.get_yaxis_transform(), fontsize=10, fontweight='bold'
    )

    # fig.show()
    for postfix in "svg", "png":
        fig.savefig(
            ld.RESULTS / f"swarm-summary-{tissue}.{postfix}",
            dpi=450, bbox_inches="tight", pad_inches=0, transparent=False
        )
    plt.close(fig)
