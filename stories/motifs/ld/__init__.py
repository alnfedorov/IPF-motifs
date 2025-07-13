from pathlib import Path

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"
RESOURCES = ROOT / "resources"


class jaspar:
    nonredundant = RESOURCES / "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    redundant = RESOURCES / "JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt"

    clusters = RESOURCES / "clusters.tab"
    parsed_clusters = RESULTS / "parsed-clusters.pkl"


class scoring:
    raw_signal = RESULTS / "signal.pkl"
    response_per_cluster = RESULTS / "clusters.pkl"
    response_per_motif = RESULTS / "motifs.pkl"
