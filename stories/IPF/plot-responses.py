import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import GRCh38
import ld
from stories.motifs import ld as motifs

plt.rcParams['svg.fonttype'] = 'none'
SAVETO = ld.RESULTS / "responses"
SAVETO.mkdir(parents=True, exist_ok=True)

# Define target motifs
MOTIFS = ['ISRE-like', 'GAS-like']

# Load the motif scores for each gene
motifs = pd.read_pickle(motifs.scoring.response)

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
for (tissue, subset), data in df.groupby(['Tissue', 'Subset']):
    for motif in MOTIFS:
        scale = 1.25
        fig, (ax_kde, ax) = plt.subplots(
            2, 1,
            sharex=True,
            figsize=(6.4 * scale, 6.4 * scale),
            gridspec_kw={'height_ratios': [1, 4], 'hspace': 0.1},
        )

        data['Has motif'] = data[motif] >= ld.thresholds.zscore

        # Scatter plot
        alpha = data['Has motif'].apply(lambda x: 1 if x else 0.35)
        color = data['Has motif'].apply(lambda x: 'blue' if x else 'red')
        sns.scatterplot(
            data=data, x='log2FoldChange', y=motif, hue='Has motif', s=10, ax=ax,
            alpha=alpha, color=color, legend=False
        )
        sns.despine(fig=fig, ax=ax)

        # KDE plot
        sns.kdeplot(
            data=data, x='log2FoldChange', hue='Has motif', fill=True, common_norm=False,
            alpha=0.35, ax=ax_kde, legend=False, palette={False: 'blue', True: 'red'}
        )
        sns.despine(fig=fig, ax=ax_kde)

        # Add titles and limits
        N = len(data)
        ax_kde.set_title(f"{tissue}: {subset} [N={N}]")

        ax.axvline(0, color='black', lw=2)
        ax.axhline(ld.thresholds.zscore, color='black', lw=2)
        ax.set(xlim=(-5, 5), ylim=(-4, 4))

        # Annotate medians on the KDE plot
        for has_motif, group in data.groupby('Has motif'):
            color, ha = ('red', 'left') if has_motif else ('blue', 'right')
            fc = group['log2FoldChange'].median()
            ax_kde.text(
                fc, ax_kde.get_ylim()[1], f"{fc:.2f}", color=color, ha=ha, va='top', fontsize=8
            )

        # Label the genes in upper left and right quadrants
        for masking, x, y in [
            (lambda x: (x['log2FoldChange'] < 0) & (x[motif] > ld.thresholds.zscore), 0, 1),
            (lambda x: (x['log2FoldChange'] > 0) & (x[motif] > ld.thresholds.zscore), 1, 1),
        ]:
            mask = masking(data)
            selected = data[mask].copy()
            selected['rank'] = selected[motif].rank(ascending=False) + \
                               selected['log2FoldChange'].abs().rank(ascending=False)
            for gid, row in selected.nsmallest(50, 'rank').iterrows():
                gname = genes[gid].attrs.name
                ax.text(row['log2FoldChange'], row[motif], gname, fontsize=6, ha='center', va='center')

        fig.savefig(
            SAVETO / f"{tissue}-{subset}-{motif}-response.svg", dpi=300, bbox_inches="tight",
            pad_inches=0, transparent=True
        )
        fig.show()
        plt.close(fig)
