from collections import defaultdict

import pandas as pd
from biobit.core.loc import Orientation
from biobit.io.bed import Bed6
from intervaltree import IntervalTree
from joblib import Parallel, delayed

import GRCh38
import ld
import utils

ANNOTATION = GRCh38.gencode.load()

tss_index = defaultdict(IntervalTree)
transcripts_index = defaultdict(IntervalTree)
all_transcripts = []

for rna in ANNOTATION.rnas.values():
    is_well_defined = len(rna.attrs.tags & {'MANE Select', 'MANE Plus Clinical'}) > 0
    is_well_defined |= rna.attrs.TSL == 1 or rna.attrs.TSL == 2
    if not is_well_defined:
        continue

    tid = rna.ind.split(".")[0]
    gene = ANNOTATION.genes[rna.gene]
    gid = gene.ind.split(".")[0]

    if gene.attrs.type != "protein_coding":
        continue

    record = {
        "Transcript ID": tid, "Gene ID": gid, "Gene name": gene.attrs.name, "Gene Type": gene.attrs.type,
        "seqid": rna.loc.seqid, "start": rna.loc.start, "end": rna.loc.end, "strand": rna.loc.strand,
    }

    tss = rna.loc.start if rna.loc.strand == "+" else rna.loc.end
    tss_index[rna.loc.seqid].addi(tss - 1, tss + 1, data=record)
    transcripts_index[rna.loc.seqid].addi(rna.loc.start, rna.loc.end, data=record)

    all_transcripts.append((rna, gene))

# Load all cCREs from ENCODE and select only those that are relevant for the analysis
df = pd.read_csv(GRCh38.cCRE, sep='\t', names=["seqid", "start", "end", "ind", "name", "ccre"])
df['ccre'] = df['ccre'].apply(lambda x: set(x.split(',')))

ROIs = {}
for roi, saveto, color in [
    ("PLS", ld.cCRE.PLS, "255,0,0"),
    ("pELS", ld.cCRE.pELS, "0,255,0"),
    ("DNase-H3K4me3", ld.cCRE.DNase_H3K4me3, "0,0,255")
]:
    ccre = df[df['ccre'].apply(lambda x: roi in x)]
    ccre = ccre[['seqid', 'start', 'end', 'name']].copy()

    # Save as BED
    bed = [
        Bed6(seqid, (start, end), name, 0, Orientation.Dual)
        for seqid, start, end, name in ccre.itertuples(index=False)
    ]
    utils.bed.tbindex(bed, Bed6, saveto)

    # Save the ROIs partitioned by contig
    for seqid, partition in ccre.groupby('seqid'):
        ROIs[seqid, roi] = partition


# Map ROIs to transcripts
def job(roi: str, seqid: str, offset: int, coordinates: pd.DataFrame, index: IntervalTree):
    results = []
    for roi_seqid, start, end, name in coordinates.itertuples(index=False):
        assert roi_seqid == seqid, f"{roi_seqid} != {seqid}"

        for rna in index.overlap(start - offset, end + offset):
            assert rna.data["seqid"] == roi_seqid, f"{rna.data['seqid']} != {roi_seqid}"

            results.append({
                "Transcript ID": rna.data["Transcript ID"],
                "Gene ID": rna.data["Gene ID"], "Gene name": rna.data["Gene name"],
                "seqid": roi_seqid,
                "rna-start": rna.data["start"], "rna-end": rna.data["end"],
                "rna-strand": rna.data["strand"],
                "roi-type": roi, "roi-start": start, "roi-end": end, "roi-name": name
            })
    return results


indices = {"PLS": tss_index, "pELS": tss_index, "DNase-H3K4me3": transcripts_index}
results = Parallel(n_jobs=-1)(
    delayed(job)(roi, seqid, ld.cCRE.overlaps.max_distances[roi], coordinates, indices[roi][seqid])
    for (seqid, roi), coordinates in ROIs.items()
)
results = pd.concat([pd.DataFrame(x) for x in results]).drop_duplicates()
results["imputed"] = False

# Impute promoters for transcripts without matched ENCODE PLS
records = []
has_pls = set(results.loc[results['roi-type'] == "PLS", "Transcript ID"])
primary_contigs = set(f"chr{i}" for i in range(1, 23)) | {"chrX", "chrY"}
for rna, gene in all_transcripts:
    ind = rna.ind.split('.')[0]
    if ind in has_pls or rna.loc.seqid not in primary_contigs:
        continue

    if rna.loc.strand == "+":
        start = rna.loc.start - ld.sequences.lengths["PLS"] // 2
        end = rna.loc.start + ld.sequences.lengths["PLS"] - (rna.loc.start - start)
    else:
        assert rna.loc.strand == '-'
        start = rna.loc.end - ld.sequences.lengths["PLS"] // 2
        end = rna.loc.end + ld.sequences.lengths["PLS"] - (rna.loc.end - start)
    assert end - start == ld.sequences.lengths["PLS"] and start < end
    if start < 0:
        continue

    records.append({
        "Transcript ID": ind,
        "Gene ID": gene.ind.split('.')[0], "Gene name": gene.attrs.name,
        "seqid": rna.loc.seqid, "rna-start": rna.loc.start, "rna-end": rna.loc.end, "rna-strand": rna.loc.strand,
        "roi-type": "PLS", "roi-start": start, "roi-end": end, "roi-name": f"Imputed[{ind}]", "imputed": True
    })
print(f"Total parsed transcripts: {len(all_transcripts):,}")
print(f"\tImputed PLS: {len(records):,}")
print(f"\tENCODE PLS: {len(has_pls):,}")
imputed = pd.DataFrame(records)

# Combine the data
assert (results.columns == imputed.columns).all(), set(results.columns).symmetric_difference(imputed.columns)
results = pd.concat([results, imputed])

# Save as pkl
ld.cCRE.overlaps.saveto.mkdir(exist_ok=True, parents=True)
results.to_pickle(ld.cCRE.overlaps.pkl)

# Save as BED
bed = []
for tid, gname, roi, seqid, start, end, imputed in results[
    ["Transcript ID", "Gene name", "roi-type", "seqid", "roi-start", "roi-end", "imputed"]
].itertuples(index=False, name=None):
    imputed = "Impute" if imputed else "cCRE"
    bed.append(Bed6(
        seqid, (start, end), f"{roi}-{gname}-{tid}-{imputed}", 0, Orientation.Dual,
    ))
utils.bed.tbindex(bed, Bed6, ld.cCRE.overlaps.bed)
