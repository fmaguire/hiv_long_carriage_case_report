#!/bin/bash
set -euo pipefail

# following this protocol
# https://www.protocols.io/view/introduction-to-calculating-dn-ds-ratios-with-code-qhwdt7e.html

# align genomes with nextclade and get amino acid alignments
# nextclade v1.10.3 DB:2022-02-07T12:00:00Z
cat ../genomes/JD-* > input/hiv_genomes.fasta
./nextclade-Linux-x86_64 dataset get --name sars-cov-2  --output-dir database
./nextclade-Linux-x86_64 --input-fasta input/hiv_genomes.fasta --include-reference --output-basename nextclade_hiv --output-dir nextclade_hiv --input-dataset database

# create concatenate nt and aa alignments of coding regions of genome
python create_concatenated_alignment.py --gff3 database/genemap.gff --nextclade nextclade_hiv

# run pal2nal (v14) to prepare for codeml
pal2nal.v14/pal2nal.pl concatenated_coding_alignment.faa concatenated_coding_alignment.fna -output paml -nogap > hiv_concatenated.pal2nal

# run codeml
codeml 

# parse output to generate codeml_kaks_values.tsv 
python parse_codeml_output.py hiv_codeml.txt
