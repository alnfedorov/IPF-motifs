import pickle

import pandas as pd
from biobit import io
from joblib import Parallel, delayed

import GRCh38
import ld
from stories.cCRE import ld as cCRE
from utils import motifs, fasta


def screeen(seqid: str, start: int, end: int, allmotifs: motifs.ZeroOrderMotifsCollection):
    reader = io.fasta.IndexedReader(GRCh38.fasta)
    forward = fasta.fetch(reader, seqid, (start, end), strand='+')
    revcomp = fasta.fetch(reader, seqid, (start, end), strand='-')

    # Get scores for all motifs
    scores = motifs.score(forward.upper(), revcomp.upper(), allmotifs)
    scores = {(motif.ind, motif.target): score for motif, score in zip(allmotifs.motifs, scores)}
    return (seqid, start, end), scores


# Load all motifs and calculate PWMs
with open(ld.jaspar.nonredundant) as f:
    jaspar = motifs.parse.jaspar(f)

pwms = tuple(pfm.to_ppm().to_pwm() for pfm in jaspar.motifs)
database = motifs.ZeroOrderMotifsCollection("ACGT", motifs=pwms)

# Load promoter sequences
sequences = pd.read_pickle(cCRE.sequences.saveto)
sequences = sequences[sequences['roi-type'] == 'PLS']
regions = sequences[['seqid', 'roi-norm-start', 'roi-norm-end']].drop_duplicates()

for seqid, group in regions.groupby('seqid'):
    print(f"Processing {len(group)} regions on {seqid}")

    result = Parallel(n_jobs=-1, verbose=100, pre_dispatch='all', batch_size=128)(
        delayed(screeen)(seqid, start, end, database)
        for seqid, start, end in group.itertuples(index=False, name=None)
    )
    repacked = dict(result)

    print(f"Saving {len(repacked)} regions on {seqid}")
    saveto = ld.scoring.raw_signal / f"{seqid}.pkl"
    saveto.parent.mkdir(parents=True, exist_ok=True)
    with open(saveto, 'wb') as f:
        pickle.dump(repacked, f, protocol=pickle.HIGHEST_PROTOCOL)

    del result, repacked
