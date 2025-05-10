from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

gff3 = ROOT / "gencode.v46.primary_assembly.annotation.gff3.gz"
index = ROOT / "gencode.v46.primary_assembly.annotation.index.pkl"

GeneSource = Literal['ENSEMBL', 'HAVANA']

GeneType = Literal[
    'transcribed_processed_pseudogene', 'snoRNA', 'TEC', 'vault_RNA', 'Mt_rRNA', 'scRNA', 'TR_J_pseudogene', 'scaRNA',
    'pseudogene', 'processed_pseudogene', 'misc_RNA', 'IG_J_pseudogene', 'TR_V_gene',
    'transcribed_unprocessed_pseudogene', 'IG_V_gene', 'TR_V_pseudogene', 'artifact', 'lncRNA', 'rRNA',
    'protein_coding', 'snRNA', 'unitary_pseudogene', 'IG_pseudogene', 'TR_J_gene', 'ribozyme', 'IG_C_gene',
    'IG_V_pseudogene', 'miRNA', 'sRNA', 'unprocessed_pseudogene', 'translated_processed_pseudogene', 'Mt_tRNA',
    'TR_D_gene', 'IG_D_gene', 'rRNA_pseudogene', 'transcribed_unitary_pseudogene', 'TR_C_gene', 'IG_J_gene',
    'IG_C_pseudogene'
]
GeneLevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrGene:
    source: GeneSource
    level: GeneLevel
    name: str
    type: GeneType

    def __attrs_post_init__(self):
        if self.type not in GeneType.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if self.source not in GeneSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


RNASource = Literal['ENSEMBL', 'HAVANA']

RNAType = Literal[
    'lncRNA', 'protein_coding_LoF', 'scRNA', 'translated_processed_pseudogene', 'nonsense_mediated_decay',
    'unitary_pseudogene', 'IG_D_gene', 'Mt_tRNA', 'TR_V_gene', 'IG_V_gene', 'processed_pseudogene',
    'transcribed_unprocessed_pseudogene', 'TR_J_gene', 'non_stop_decay', 'transcribed_unitary_pseudogene', 'ribozyme',
    'protein_coding_CDS_not_defined', 'IG_J_pseudogene', 'artifact', 'IG_J_gene', 'transcribed_processed_pseudogene',
    'sRNA', 'snRNA', 'scaRNA', 'rRNA', 'IG_C_gene', 'IG_V_pseudogene', 'processed_transcript', 'vault_RNA',
    'rRNA_pseudogene', 'unprocessed_pseudogene', 'misc_RNA', 'TR_D_gene', 'protein_coding', 'TEC', 'TR_V_pseudogene',
    'TR_C_gene', 'snoRNA', 'Mt_rRNA', 'pseudogene', 'retained_intron', 'miRNA', 'IG_pseudogene', 'TR_J_pseudogene',
    'IG_C_pseudogene'
]
RNATag = Literal[
    'Ensembl canonical', 'GENCODE basic', 'GENCODE primary', 'MANE Select', 'MANE Plus Clinical', 'alternative_3_UTR',
    'non_canonical_other', 'mRNA_start_NF', 'mRNA_end_NF', 'RNA_Seq_supported_only', 'bicistronic', 'upstream_uORF',
    'CCDS', 'non_canonical_U12', 'appris_principal_4', 'exp_conf', '3_nested_supported_extension',
    'non_canonical_genome_sequence_error', '454_RNA_Seq_supported', 'NAGNAG_splice_site', 'appris_principal_2',
    'inferred_transcript_model', 'seleno', 'downstream_ATG', 'non_submitted_evidence', 'dotter_confirmed',
    'retained_intron_CDS', 'cds_end_NF', 'non_canonical_TEC', 'alternative_5_UTR', 'inferred_exon_combination',
    'non_ATG_start', 'stop_codon_readthrough', 'retained_intron_final', 'appris_alternative_2', 'TAGENE',
    'RNA_Seq_supported_partial', 'non_canonical_conserved', 'appris_principal_3', 'overlapping_uORF',
    'low_sequence_quality', 'NMD_exception', 'nested_454_RNA_Seq_supported', '5_standard_supported_extension',
    'CAGE_supported_TSS', 'NMD_likely_if_extended', 'upstream_ATG', '3_standard_supported_extension',
    'appris_principal_5', 'pseudo_consens', '5_nested_supported_extension', 'retained_intron_first',
    'non_canonical_polymorphism', 'not_organism_supported', 'sequence_error', 'appris_principal_1', 'RP_supported_TIS',
    'cds_start_NF', 'readthrough_transcript', 'appris_alternative_1', 'not_best_in_genome_evidence'
]
RNATSL = Literal[1, 2, 3, 4, 5, None]
RNALevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrRNA:
    source: RNASource
    level: RNALevel
    name: str
    type: RNAType
    tags: frozenset[RNATag]
    TSL: RNATSL
    CDS: frozenset[str]

    def __attrs_post_init__(self):
        if self.source not in RNASource.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.type not in RNAType.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if any(x not in RNATag.__args__ for x in self.tags):
            raise ValueError(f"Invalid tags: {self.tags}")
        if self.TSL not in RNATSL.__args__:
            raise ValueError(f"Invalid TSL: {self.TSL}")


CDSSource = Literal['ENSEMBL', 'HAVANA']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrCDS:
    source: CDSSource
    transcripts: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in CDSSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


def load() -> at.Annotome[AttrGene, AttrRNA, AttrCDS]:
    return at.read_pkl(index.as_posix())
