import pickle

import pandas as pd

import ld

foldchanges = {}
for tissue, data in ld.resources.subtypes.items():
    print(tissue)
    for celltype, path in data.items():
        df = pd.read_csv(path, index_col=0)
        significant = set(df[df['padj'] <= 0.1].index)

        # Select only genes where estimated fold change is 'reliable'
        df = df[['log2FoldChange', 'lfcSE', 'baseMean']].dropna(how='any')
        df['metric'] = df['log2FoldChange'].abs() / df['lfcSE']
        df = df[df['metric'] >= 1].copy()

        # Sanity check
        assert set(df.index) & significant == significant

        foldchanges[tissue, celltype] = df['log2FoldChange']
        print(f"\t{celltype}: {len(df)} genes with log2(fold change) > log2(SE)")

# Map genes to Ensembl IDs
mapping = pd.read_csv(ld.resources.synonyms, dtype=str, sep='\t')
assert mapping["Approved symbol"].nunique() == len(mapping)
mapping = mapping.set_index("Approved symbol")["Ensembl gene ID"].to_dict()

# Extra manually resolved genes
for gname, gid in [
    ("QARS", "ENSG00000172053"), ("SSFA2", "ENSG00000138434"), ("FAM207A", "ENSG00000160256"),
    ("H2AFZ", "ENSG00000164032"), ("H2AFY", "ENSG00000113648"), ("LARS", "ENSG00000133706"),
    ("HARS", "ENSG00000170445"), ("AAED1", "ENSG00000158122"),
    ("H2AFV", "ENSG00000105968"), ("ATP5MPL", "ENSG00000156411"), ("FAM172A", "ENSG00000113391"),
    ("TROVE2", "ENSG00000116747"), ("FAM102B", "ENSG00000162636"), ("WDR60", "ENSG00000126870"),
    ("HIST1H1E", "ENSG00000168298"), ("HIST1H4C", "ENSG00000197061"), ("GATD3A", "ENSG00000160221"),
]:
    assert gname not in mapping
    mapping[gname] = gid

print("\nGenes mapped to Ensembl IDs")
for celltype, df in foldchanges.items():
    mask = df.index.isin(mapping.keys())
    print(f"{celltype}: {mask.sum() / len(df):.2%} genes")

    df = df[mask].copy()
    df.index = df.index.map(mapping)
    foldchanges[celltype] = df

# Save the classification
ld.SCRNA_FOLD_CHANGE.parent.mkdir(exist_ok=True, parents=True)
with open(ld.SCRNA_FOLD_CHANGE, "wb") as stream:
    pickle.dump(foldchanges, stream)
