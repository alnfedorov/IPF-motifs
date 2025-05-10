from pathlib import Path

from . import gencode, seqid

name = "GRCh38"
organism = "Homo sapiens"

ROOT = Path(__file__).parent

fasta = ROOT / "GRCh38.primary_assembly.genome.fa.bgz"
cCRE = ROOT / "GRCh38-cCREs.bed.gz"  # ENCODE cCRE version 3
