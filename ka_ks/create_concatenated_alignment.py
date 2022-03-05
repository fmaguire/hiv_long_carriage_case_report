#!/usr/bin/env python

import argparse
from operator import itemgetter
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Script to generate concatenated whole genome nt and amino acid alignments from nextclade output')
    parser.add_argument('--gff3', required=True, help='Gene map used by nextclade (included in nextclade db as genemap.gff')
    parser.add_argument('--nextclade', required=True, help="Nextclade results folder")
    parser.add_argument('--include_stop', help="Include stop codons in concat",  action='store_false')
    args = parser.parse_args()

    # parse the gene location indices and sort by start
    gene_locations = {}
    with open(args.gff3) as fh:
        for line in fh:
            if not line.startswith('#'):
                line = line.strip().split()
                gene_name = line[-1].replace('gene_name=', '')
                # exclude ORF7ab as there are too many gaps/translation errors
                if not gene_name.startswith('ORF7'):
                    gene_locations[gene_name] = (int(line[3]), int(line[4]))
    gene_order = sorted([(k,*v) for k,v in gene_locations.items()], key=itemgetter(1))
    gene_order = [x[0] for x in gene_order]

    # parse full alignment and store IDs in order
    nextclade = Path(args.nextclade)
    genome_nt_afa_fp = list(nextclade.glob('*.aligned.fasta'))[0]
    with open(genome_nt_afa_fp) as fh:
        genome_nt_afa = AlignIO.read(fh, 'fasta')

    # for each gene concatenate the AA seq and extract+concatenate the coding nt
    coding_nt_afa_list = []
    coding_aa_afa_list = []

    for gene in gene_order:
        # slice the nt alignment using gene locations
        start, stop = gene_locations[gene]
        if args.include_stop:
            coding_nt_afa_list.append(genome_nt_afa[:, start - 1: stop])
        else:
            coding_nt_afa_list.append(genome_nt_afa[:, start - 1: stop - 3])
        #coding_nt_afa_list.append(genome_nt_afa[:, start - 1: stop])

        # parse the aa alignments
        gene_aa_fp = list(nextclade.glob(f"*.gene.{gene}.fasta"))[0]
        if gene_aa_fp.exists():
            with open(gene_aa_fp) as fh:
                gene_aa_afa = AlignIO.read(fh, 'fasta')
                if args.include_stop:
                    coding_aa_afa_list.append(gene_aa_afa)
                else:
                    gene_aa_afa = gene_aa_afa[:,:-1]
                    coding_aa_afa_list.append(gene_aa_afa)
                #coding_aa_afa_list.append(gene_aa_afa)
        else:
            print(gene_aa_fp)
            assert False

    coding_nt_afa = coding_nt_afa_list[0]
    coding_aa_afa = coding_aa_afa_list[0]
    for gene_ix in range(1, len(coding_nt_afa_list)):
        coding_nt_afa += coding_nt_afa_list[gene_ix]
        coding_aa_afa += coding_aa_afa_list[gene_ix]


    # write out the concatenate alignments in a consistent order
    print("Writing concatenated_coding_alignment.{faa,fna}")
    with open('concatenated_coding_alignment.fna', 'w') as out_fh:
        SeqIO.write(coding_nt_afa, out_fh, 'fasta-2line')
        #for record in coding_nt_afa:
        #    record.seq = record.seq.ungap('-')
        #    SeqIO.write(record, out_fh, 'fasta-2line')

    with open('concatenated_coding_alignment.faa', 'w') as out_fh:
        SeqIO.write(coding_aa_afa, out_fh, 'fasta-2line')


    # check ref seq
    #print(len(coding_nt_afa))
    #for record_ix in range(len(coding_nt_afa)):
    #    concat_ref_nt = coding_nt_afa[record_ix, :]
    #    concat_ref_aa = coding_aa_afa[record_ix, :]
    #    for aa, codon in enumerate(range(0, len(concat_ref_nt), 3)):
    #        ref_codon = concat_ref_nt[codon: codon+3]
    #        try:
    #            tran_ref = str(ref_codon.translate().seq)
    #        except:
    #            tran_ref = 'X'
    #        ref_codon = str(ref_codon.seq)
    #        ref_aa = str(concat_ref_aa.seq[aa])
    #        if tran_ref != ref_aa:
    #            print(concat_ref_nt.id, ref_codon, tran_ref, ref_aa)

