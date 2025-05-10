#!/bin/bash
set -euxo pipefail

VERSION=$1
SAVETO=$2

########################################################################################################################
# GRCh38 assembly
########################################################################################################################

wget -O "${SAVETO}/GRCh38.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/GRCh38.primary_assembly.genome.fa.gz"

wget -O "${SAVETO}/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz"

########################################################################################################################
# Postprocessing
########################################################################################################################

file="${SAVETO}/GRCh38.primary_assembly.genome.fa.gz"
saveto="${file/.fa.gz/.fa.bgz}"

gunzip "${file}" --stdout | bgzip --threads "$(nproc)" -l 6 -o "${saveto}" /dev/stdin
rm "${file}"
samtools faidx "${saveto}"
