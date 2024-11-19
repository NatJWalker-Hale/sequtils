#! /usr/bin/python3


import os
import re
import sys
import argparse
import gffutils


def parse_cds_from_genomic_ncbi(path: str):
    """
    Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple
    """
    with open(path, "r", encoding="utf-8") as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = re.search(r'\[protein_id=(.*?)\]', line).group(1)
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_protein_ncbi(path: str):
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path, "r", encoding="utf-8") as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:].split(" ")[0]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_fasta_lengths(path: str):
    """
    parse FASTA line by line and return length of each sequence
    """
    lengths = {}
    with open(path, "r", encoding="utf-8") as f:
        length = 0
        line = f.readline().strip()  # here we read lines into mem one at a time
        while line:
            if line.startswith(">"):
                if length:
                    lengths[seqid] += length
                seqid = line.split()[0].strip().lstrip(">")
                length = 0
                lengths[seqid] = 0
            else:
                length += len(line)
            line = f.readline().strip()
        lengths[seqid] += length # write lengths of final contig
    return lengths


def get_chromosomes(feature_db: gffutils.FeatureDB) -> dict[tuple]:
    """
    returns a dictionary of seqid: (chromosome names, gene count)
    """
    seqids = sorted(
        list(
        feature_db.seqids()
            )
                   )

    # here we are assuming input from e.g. NCBI has chr number in order of increasing seq accession
    return {(id): (f"Chr{i+1}" if "NC_" in id else id) for i, id in enumerate(seqids)}


def get_lens_and_gff(genome_fasta: str, chrom_corres: dict, feature_db: gffutils.FeatureDB,
                     pref: str) -> dict:
    """
    parse FASTA line by line and return length of each sequence
    """

    lengths = parse_fasta_lengths(genome_fasta)

    out_lengths = {}
    out_gff = []
    chrom_id = 0
    unplaced_id = 0

    for seqid, length in lengths.items():
        gene_id = 0
        if "Chr" in chrom_corres[seqid]:
            chrom_id_str = f"{chrom_id+1}"
            chrom_id += 1
        else:
            chrom_id_str = f"u{unplaced_id+1}"
            unplaced_id += 1

        for f in feature_db.region(seqid, featuretype="gene"):
            tmp = 0
            protein_coding = False
            for m in feature_db.children(f):
                if m.featuretype == "mRNA":
                    protein_coding = True
                    # print(f"{seqid} {f.id} {c.id} protein coding")

                    # width = sum(c.end - c.start for c in feature_db.children(m)
                    #              if c.featuretype == "CDS")

                    width = m.end - m.start
                    # print(f"{width}\t{m.end - m.start}")

                    if width > tmp:
                        tmp = width
                        main = m

            if protein_coding:

                p_id = next(c["protein_id"] for c in feature_db.children(main)
                            if c.featuretype == "CDS")[0]
                # p_id = "test"

                out_gff.append((chrom_corres[seqid],
                                f"{pref}{chrom_id_str}g{gene_id+1:05d}",
                                f"{main.start}",
                                f"{main.end}",
                                main.strand,
                                f"{gene_id+1}",
                                p_id))
                gene_id += 1

        out_lengths[chrom_corres[seqid]] = (length, gene_id)
    return out_lengths, out_gff


def write_outputs(lengths_dict: dict[tuple],
                  gff_list: list[tuple],
                  cds_dict: dict,
                  pep_dict: dict,
                  prefix: str):
    """
    write the .lens file for WGDI
    """
    with open(f"{prefix}.lens", "w", encoding="utf-8") as out_lens:
        for seqid, length_count in lengths_dict.items():
            out_lens.write(f"{seqid}\t{length_count[0]}\t{length_count[1]}\n")

    with open(f"{prefix}.gff", "w", encoding="utf-8") as out_gff:
        with open(f"{prefix}.cds.fa", "w", encoding="utf-8") as out_cds:
            with open(f"{prefix}.pep.fa", "w", encoding="utf-8") as out_pep:
                for entry in gff_list:
                    out_gff.write(f"{'\t'.join(entry)}\n")
                    out_cds.write(f">{entry[1]} {entry[-1]}\n{cds_dict[entry[-1]]}\n")
                    out_pep.write(f">{entry[1]} {entry[-1]}\n{pep_dict[entry[-1]]}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
                                     "script to create simplified GFFs for input to WGDI",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=
"""
writes GFF with the following spec (taken from WGDI):
Column    Information    Explanation
     1            Chr    Chromosome number
     2             Id    Gene name
     3          Start    The starting location of a gene
     4            End    The ending location of a gene
     5      Direction    Direction of a gene sequence
     6          Order    Order of each chromosome, starting from 1
     7       Original    Original id and not read

also writes a .lens file containing:
Column    Information    Explanation
     1            Chr    Chromosome number
     2         Length    Length of chromosome sequences
     3         Number    Number of chromosome genes
"""
                                    )

    parser.add_argument("GFF", help="input GFF")
    parser.add_argument("fasta", help="genome FASTA")
    parser.add_argument("CDS", help="CDS FASTA")
    parser.add_argument("pep", help="protein FASTA")
    parser.add_argument("prefix", help="string to append to output files and gene_ids")
    args = parser.parse_args(sys.argv[1:] or ["--help"])

    dbname = args.GFF.rstrip(".gff") + ".db"
    if os.path.isfile(dbname) and os.path.getsize(dbname) > 0:
        db = gffutils.FeatureDB(dbname)
    else:
        db = gffutils.create_db(args.GFF, dbname, id_spec={'gene': 'ID',
                                                           'mRNA': 'ID'},
                               )

    # print(db['gene-LOC131165771'])

    # for m in db.children('gene-LOC131165771'):
    #     if m.featuretype == "mRNA":
    #         for c in db.children(m):
    #             if c.featuretype == "CDS":
    #                 print(c["protein_id"])      

    chroms = get_chromosomes(db)

    lens, gff = get_lens_and_gff(args.fasta, chroms, db, args.prefix)

    cds_seqs = dict(parse_cds_from_genomic_ncbi(args.CDS))
    pep_seqs = dict(parse_protein_ncbi(args.pep))

    write_outputs(lens, gff, cds_seqs, pep_seqs, args.prefix)
