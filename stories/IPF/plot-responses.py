import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import GRCh38
import ld
from stories.motifs import ld as motifs

plt.rcParams['svg.fonttype'] = 'none'

# Define the target motifs and populations
MOTIFS = ['ISRE-like', 'GAS-like']
POPULATIONS = ["CXCL10+ AM", "IGF1AM"]

# Load the motif scores
motifs = pd.read_pickle(motifs.scoring.response_per_cluster)
assert (motifs['roi-type'] == 'PLS').all(), motifs['motif'].unique_values()
motifs = motifs.drop(columns=['roi-type']).set_index('Gene ID').copy()

# Load the genes of interest
genes = pd.read_pickle(ld.SCRNA_FOLD_CHANGE)

df = []
for (tissue, subset), fc in genes.items():
    before = len(fc)
    data = pd.merge(motifs, fc, how='inner', left_index=True, right_index=True)
    print(f"{subset}: {before} -> {len(data)}({len(data) / before:.2%})")

    data = data[[*MOTIFS, 'log2FoldChange']]
    data["Tissue"] = tissue
    data["Subset"] = subset
    df.append(data)
df = pd.concat(df)

# Make all plots
genes = {gene.ind.split('.')[0]: gene for gene in GRCh38.gencode.load().genes.values()}
for population in POPULATIONS:
    data = df[(df['Tissue'] == 'lung') & (df['Subset'] == population)].copy()
    for motif in MOTIFS:
        scale = 1.25
        fig, ax = plt.subplots(figsize=(6.4 * 1.25, 4.8 * 1.25))

        data['Has motif'] = data[motif] >= ld.thresholds.zscore

        alpha = data['Has motif'].apply(lambda x: 1 if x else 0.35)
        color = data['Has motif'].apply(lambda x: 'blue' if x else 'red')
        sns.scatterplot(
            data=data, x='log2FoldChange', y=motif, hue='Has motif', s=10, ax=ax,
            alpha=alpha, color=color, legend=False
        )
        sns.despine(fig=fig, ax=ax)

        N = len(data)
        ax.set_title(f"{population} [N={N}]")

        ax.axvline(0, color='black', lw=2)
        ax.axhline(ld.thresholds.zscore, color='black', lw=2)
        ax.set(xlim=(-5, 5), ylim=(-4, 4))

        # Label number of genes in each quadrant
        for masking, x, y in [
            (lambda x: (x['log2FoldChange'] < 0) & (x[motif] < ld.thresholds.zscore), 0, 0),
            (lambda x: (x['log2FoldChange'] < 0) & (x[motif] > ld.thresholds.zscore), 0, 1),
            (lambda x: (x['log2FoldChange'] > 0) & (x[motif] < ld.thresholds.zscore), 1, 0),
            (lambda x: (x['log2FoldChange'] > 0) & (x[motif] > ld.thresholds.zscore), 1, 1),
        ]:
            N = masking(data).sum()
            ha = 'left' if x == 0 else 'right'
            va = 'bottom' if y == 0 else 'top'
            ax.text(x, y, f"N={N} ({N / len(data):.1%})", transform=ax.transAxes, fontsize=12, ha=ha, va=va)

        # Label the genes in upper left and right quadrants
        for masking, x, y in [
            (lambda x: (x['log2FoldChange'] < 0) & (x[motif] > ld.thresholds.zscore), 0, 1),
            (lambda x: (x['log2FoldChange'] > 0) & (x[motif] > ld.thresholds.zscore), 1, 1),
        ]:
            for gid, row in data[masking].iterrows():
                gname = genes[gid].attrs.name
                ax.text(row['log2FoldChange'], row[motif], gname, fontsize=8, ha='center', va='center')

        fig.savefig(
            ld.RESULTS / f"{motif}-{population}-response.svg", dpi=300, bbox_inches="tight",
            pad_inches=0, transparent=True
        )
        fig.show()
        plt.close(fig)
