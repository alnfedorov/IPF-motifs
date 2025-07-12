import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

import ld

plt.rcParams['svg.fonttype'] = 'none'

PALETTE = {
    'Not significant': '#93a2cb',
    'Significant': '#da6434'
}

X = '-log10(p-value)'

df = pd.read_pickle(ld.STAT_TESTS)
df['-log10(p-value)'] = -df['p-value'].apply(lambda x: np.log10(x))

for tissue, subdf in df.groupby('tissue'):
    subdf = subdf.sort_values(by=['subset'])
    subdf['q-value'] = stats.false_discovery_control(subdf['p-value'], method='bh')

    subdf['Category'] = 'Not significant'

    mask = subdf['q-value'] <= ld.thresholds.qvalue
    subdf.loc[mask, 'Category'] = 'Significant'

    # subdf['Median effect'] = 'None/Low'
    # subdf.loc[subdf['Median log2 fold change'] < -0.585, 'Median effect'] = 'Downregulated in IPF'
    # subdf.loc[subdf['Median log2 fold change'] > 0.585, 'Median effect'] = 'Upregulated in IPF'

    total_genes = {}
    for subset, subsubdf in subdf.groupby('subset'):
        total = subsubdf['Have motif'] + subsubdf['No motif']
        assert total.nunique() == 1
        total_genes[subset] = total.max()
    subdf['subset'] = subdf['subset'].apply(lambda x: f"{x}\n[Genes={total_genes[x]}]")

    # order = [f"{subset}\n[Genes={total_genes[subset]}]" for subset in ORDER]
    order = sorted(subdf['subset'].unique())
    assert set(order) == set(subdf['subset']), set(subdf['subset']) - set(order)

    subdf['Y'] = subdf['subset'].map({subset: i for i, subset in enumerate(order)})
    subdf['Y'] += np.random.uniform(-0.2, 0.2, size=len(subdf))

    # subdf['Size'] = subdf['Have motif']
    # subdf['Size'] = subdf['Have motif'] / (subdf['Have motif'] + subdf['No motif'])
    # subdf['Size'] = (subdf['Size'] / subdf['Size'].max() * 5).clip(2, 5)

    height = 0.5 * len(order) + 0.5
    fig, ax = plt.subplots(figsize=(14.8, height))
    sns.scatterplot(
        data=subdf, y='Y', x=X, ax=ax, hue='Category', s=50, edgecolor='white', linewidth=0.75, palette=PALETTE,
    )

    ax.set_yticks(np.arange(len(order)), order)
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylabel(None)
    ax.grid(axis='y')

    # Label ISRE/GAS-like motifs and all significant motifs
    for y, subset in enumerate(order):
        notsignif = subdf[(subdf['subset'] == subset) & (subdf['Category'] == 'Not significant')]
        for motif in "ISRE-like", "GAS-like":
            x = notsignif[notsignif['motif'] == motif]
            if x.empty:
                continue
            x = x[X].values[0]
            if x > 0:
                ax.text(x, y, motif, ha='left', va='center', color='red')

        signif = subdf[(subdf['subset'] == subset) & (subdf['Category'] != 'Not significant')]
        for motif, corr, signif in signif[['motif', X, 'Category']].itertuples(index=False):
            fontsize = 6 if motif.startswith("cluster") else 12
            ax.text(
                corr, y, motif, ha='left', va='center', color='black', fontweight='bold', fontsize=fontsize
            )

    # FDR line
    minpval = subdf[subdf['Category'] != 'Not significant'][X].min()
    ax.axvline(minpval, color='black', linestyle='--', linewidth=1)
    ax.text(
        minpval, 0, f'FDR ≤ {ld.thresholds.qvalue}', ha='left', va='bottom', color='black',
        transform=ax.get_xaxis_transform(), fontsize=10, fontweight='bold'
    )

    fig.show()
    ld.RESULTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        ld.RESULTS / f"swarm-summary-{tissue}.svg", dpi=450, bbox_inches="tight", pad_inches=0, transparent=True
    )
    plt.close(fig)
