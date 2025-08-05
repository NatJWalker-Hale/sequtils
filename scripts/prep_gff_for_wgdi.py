#! /usr/bin/python3


import os
import re
import sys
import argparse
import gffutils
import sequence as sq


def get_chromosomes(feature_db: gffutils.FeatureDB) -> dict[tuple]:
    """
    returns a dictionary of seqid: chromosome names
    """
    seqids = sorted(
        list(
        feature_db.seqids()
            )
                   )

    # here we are assuming input from e.g. NCBI has chr number in order of increasing seq accession
    return {(id): (f"Chr{i+1}" if "NC_" in id or "CM" in id else id) for i, id in enumerate(seqids)}


def get_chromosomes_gwh(path: str):
    """
    parse FASTA headers to find sequence IDs with chromosome info
    """
    with open(path, "r", encoding="utf-8") as f:
        name = chrom = ""
        line = f.readline().strip()  # here we read lines into mem one at a time
        while line:
            if line.startswith(">"):
                if name:
                    yield name, chrom
                name = line[1:].split("\t")[0]
                chrom = re.search(r'OriSeqID=(.*?)\t', line).group(1)
            line = f.readline().strip()
        yield name, chrom


def get_chromosomes_ensembl(path: str):
    """
    parse FASTA headers to find sequence IDs with chromosome info
    """
    with open(path, "r", encoding="utf-8") as f:
        name = chrom = ""
        line = f.readline().strip()  # here we read lines into mem one at a time
        while line:
            if line.startswith(">"):
                if name:
                    yield name, chrom
                name = line[1:].split(" ")[0]
                chrom = f"Chr{name}"
            line = f.readline().strip()
        yield name, chrom


def get_chromosomes_fasta(path: str):
    """
    parse FASTA headers to find sequence IDs with chromosome info
    """
    with open(path, "r", encoding="utf-8") as f:
        name = chrom = ""
        line = f.readline().strip()  # here we read lines into mem one at a time
        i = 0
        while line:
            if line.startswith(">"):
                i += 1
                if name:
                    yield name, chrom
                name = line[1:].split("\t")[0]
                chrom = f"Chr{i}" if any(x in name for x in ["NC_", "CM", "chr", "Chr"]) else name
                # chrom = f"Chr{i}" if "NC_" in name or "CM" in name or "chr" in name else name
            line = f.readline().strip()
        yield name, chrom


def parse_cds_and_pep(cds_path: str, pep_path: str, ncbi=True, phytozome=False,
                      gwh=False, ensembl=False, maker=False) -> dict[str]:
    """
    parse input sequence files
    """
    if ncbi:
        cds_dict = dict(sq.parse_cds_from_genomic_ncbi(cds_path))
        pep_dict = dict(sq.parse_protein_ncbi(pep_path))
    elif phytozome:
        cds_dict = dict(sq.parse_cds_phytozome(cds_path))
        pep_dict = dict(sq.parse_protein_phytozome(pep_path))
    elif gwh:
        cds_dict = dict(sq.parse_fasta_gwh(cds_path))
        pep_dict = dict(sq.parse_fasta_gwh(pep_path))
    elif ensembl:
        cds_dict = dict(sq.parse_cds_ensembl(cds_path))
        pep_dict = dict(sq.parse_protein_ensembl(pep_path))
    elif maker:
        cds_dict = dict(sq.parse_fasta_first_header(cds_path))
        pep_dict = dict(sq.parse_fasta_first_header(pep_path))
    else:
        cds_dict = dict(sq.parse_fasta_first_header(cds_path))
        pep_dict = dict(sq.parse_fasta_first_header(pep_path))

    return cds_dict, pep_dict


def get_lens_and_gff(genome_fasta: str, chrom_corres: dict, feature_db: gffutils.FeatureDB,
                     pref: str, get_lengths=True, ncbi=True, phytozome=False, gwh=False,
                     ensembl=False, maker=False) -> dict:
    """
    UPDATE
    """

    # get the lengths from the genome fasta
    if get_lengths:
        lengths = sq.parse_fasta_lengths(genome_fasta)
    else:
        lengths = {v: 0 for v in chrom_corres.values()}

    out_lengths = {}
    out_gff = []
    chrom_id = 0
    unplaced_id = 0

    # iterate over the seqids in lengths
    for seqid, length in lengths.items():
        gene_id = 0
        if "Chr" in chrom_corres[seqid] or "chr" in chrom_corres[seqid]:
            # naming scheme for genes in chromosomes
            chrom_id_str = f"{chrom_id+1}"
            chrom_id += 1
        else:
            # naming scheme for unplaced scaffolds
            chrom_id_str = f"u{unplaced_id+1}"
            unplaced_id += 1

        # iterate over genes in the current sequence
        for f in feature_db.region(seqid, featuretype="gene"):
            tmp = 0
            p_id = ""
            protein_coding = False
            # iterate over children of gene entries
            for m in feature_db.children(f):
                # filter e.g. non-coding in NCBI
                if m.featuretype == "mRNA":
                    protein_coding = True

                    # use phytozomes internal designation of primary transcript for consistency
                    # with primary only FASTAs
                    if phytozome:
                        # print(m["longest"])
                        if m["longest"][0] == '1':
                            main = m
                            # break out
                            break
                    # print(f"{seqid} {f.id} {c.id} protein coding")

                    # otherwise, calculate length of alt transcripts as sum of CDS element lengths
                    width = sum(c.end - c.start for c in feature_db.children(m)
                                 if c.featuretype == "CDS")

                    # width = m.end - m.start
                    # print(f"{width}\t{m.end - m.start}")

                    # pick the longest transcript
                    if width > tmp:
                        tmp = width
                        main = m

            # skip e.g. ncRNA
            if protein_coding:

                if ncbi:
                    # extract protein ID from CDS elements, these are identical in NCBI formats
                    p_id = next(c["protein_id"] for c in feature_db.children(main)
                                if c.featuretype == "CDS")[0]
                    # p_id = "test"
                if phytozome:
                    # ID that matches phytozome FASTAs is element Name of Attributes
                    p_id = main["Name"][0]
                    # print(p_id)
                if gwh or maker:
                    p_id = main.id
                if ensembl:
                    p_id = main["transcript_id"][0] + ".1"

                out_gff.append((chrom_corres[seqid],  # replace seqid with new chrom designation
                                f"{pref}{chrom_id_str}g{gene_id+1:05d}", # 1-indexed
                                f"{main.start}",  # use mRNA coordinates
                                f"{main.end}",
                                main.strand,
                                f"{gene_id+1}", # 1-indexed
                                p_id))
                gene_id += 1

        out_lengths[chrom_corres[seqid]] = (length, gene_id)
    return out_lengths, out_gff


def write_outputs(lengths_dict: dict[tuple],
                  gff_list: list[tuple],
                  cds_dict: dict,
                  pep_dict: dict,
                  prefix: str,
                  write_lens = True):
    """
    write the .lens file for WGDI
    """
    if write_lens:
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

finally, modifies the provided CDS and peptide FASTAs to include the 
new and the old names
"""
                                    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ncbi", help="mode for NCBI formatted files", action="store_true")
    group.add_argument("--phytozome", help="mode for Phytozome formatted files",
                       action="store_true")
    group.add_argument("--gwh", help="mode for Genome Warehouse formatted files",
                       action="store_true")
    group.add_argument("--ensembl", help="mode for ensembl formatted files", action="store_true")
    group.add_argument("--maker", help="mode for maker formatted files", action="store_true")
    parser.add_argument("--no_lens", help="do not parse genome or write .lens file",
                        action="store_false")
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
    cds_seqs, pep_seqs = parse_cds_and_pep(args.CDS, args.pep, args.ncbi, args.phytozome, args.gwh,
                                           args.ensembl, args.maker)

    if args.ncbi or args.phytozome:
        chroms = get_chromosomes(db)
    elif args.ensembl:
        chroms = dict(get_chromosomes_ensembl(args.fasta))
    elif args.maker:
        chroms = dict(get_chromosomes_fasta(args.fasta))
    else:
        chroms = dict(get_chromosomes_gwh(args.fasta))

    print(chroms)

    lens, gff = get_lens_and_gff(args.fasta, chroms, db, args.prefix, args.no_lens, args.ncbi,
                                 args.phytozome, args.gwh, args.ensembl, args.maker)

    write_outputs(lens, gff, cds_seqs, pep_seqs, args.prefix, args.no_lens)
