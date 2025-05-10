from pathlib import Path

ROOT = Path(__file__).parent
RESULTS = ROOT / "results"

STAT_TESTS = RESULTS / "stat_tests.pkl"
SCRNA_FOLD_CHANGE = RESULTS / "fold-change.pkl"


class thresholds:
    qvalue = 0.1
    zscore = 2


class resources:
    folder = ROOT / "resources"
    synonyms = folder / "HGNC.tsv"

    subtypes = {
        "lung": {
            "CD9+_FN1_CD206hiAMs": folder / "lung" / "CD9+_FN1_CD206hiAMs_DEG_0.1.csv",
            "CXCL10+ AM": folder / "lung" / "CXCL10+ AM_DEG_0.1.csv",
            "Cycling AM 1": folder / "lung" / "Cycling AM 1_DEG_0.1.csv",
            "Cycling AM 2": folder / "lung" / "Cycling AM 2_DEG_0.1.csv",
            "Effector CD4+ T cells": folder / "lung" / "Effector CD4+ T cells_DEG_0.1.csv",
            "FABP4AM": folder / "lung" / "FABP4AM_DEG_0.1.csv",
            "IGF1AM": folder / "lung" / "IGF1AM_DEG_0.1.csv",
            "Naive CD4+ T cells": folder / "lung" / "Naive CD4+ T cells_DEG_0.1.csv",
            "SPP1+ AM": folder / "lung" / "SPP1+ AM_DEG_0.1.csv",
        }
    }
