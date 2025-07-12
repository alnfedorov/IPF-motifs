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
        },
        "integrated": {
            "AM_UD": folder / "integrated" / "AM_UD_DEG_wRiboGenes.csv",
            "CD206hi FN1hi AM": folder / "integrated" / "CD206hi FN1hi AM_DEG_wRiboGenes.csv",
            "cDC1": folder / "integrated" / "cDC1_DEG_wRiboGenes.csv",
            "cDC2": folder / "integrated" / "cDC2_DEG_wRiboGenes.csv",
            "cMono": folder / "integrated" / "cMono_DEG_wRiboGenes.csv",
            "CXCL10+ Mono-Mac": folder / "integrated" / "CXCL10+ Mono-Mac_DEG_wRiboGenes.csv",
            "Cycling Macs": folder / "integrated" / "Cycling Macs_DEG_wRiboGenes.csv",
            "FABP4hi AM": folder / "integrated" / "FABP4hi AM_DEG_wRiboGenes.csv",
            "HSP+ Mac": folder / "integrated" / "HSP+ Mac_DEG_wRiboGenes.csv",
            "IGF-1+ AM": folder / "integrated" / "IGF-1+ AM_DEG_wRiboGenes.csv",
            "IM": folder / "integrated" / "IM_DEG_wRiboGenes.csv",
            "Intm AM": folder / "integrated" / "Intm AM_DEG_wRiboGenes.csv",
            "ISGhi AM": folder / "integrated" / "ISGhi AM_DEG_wRiboGenes.csv",
            "MonoDC": folder / "integrated" / "MonoDC_DEG_wRiboGenes.csv",
            "MonoMac": folder / "integrated" / "MonoMac_DEG_wRiboGenes.csv",
            "MT_AM": folder / "integrated" / "MT_AM_DEG_wRiboGenes.csv",
            "ncMono": folder / "integrated" / "ncMono_DEG_wRiboGenes.csv",
            "pDC": folder / "integrated" / "pDC_DEG_wRiboGenes.csv",
            "SPP1+ Mac": folder / "integrated" / "SPP1+ Mac_DEG_wRiboGenes.csv",
        }
    }
