#!/usr/bin/env python3


"""
utilties for parsing, writing and processing common sequence file formats
"""


import re
from collections import Counter
import numpy as np


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
                if re.search(r'\[pseudo=true\]', line):
                    continue
                name = re.search(r'\[protein_id=(.*?)\]', line).group(1)
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_cds_phytozome(path: str):
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
                name = line[1:].split(" ")[0]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_protein_ncbi(path: str):
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
                name = line[1:].split(" ")[0]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_protein_phytozome(path: str):
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
                name = re.search(r'transcript=(.*?) ', line).group(1)
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_fasta_gwh(path: str):
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
                name = re.search(r'OriID=(.*?)\t', line).group(1)
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_cds_ensembl(path: str):
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
                name = line[1:].split(" ")[0]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_protein_ensembl(path: str):
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
                name = re.search(r'transcript:(.*?) ', line).group(1)
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_fasta_first_header(path: str):
    """
    courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846

    given a path tries to parse a fasta file. Returns an iterator which yields a (name, sequence) 
    tuple
    """
    with open(path, "r", encoding="utf-8") as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:].split(" ")[0]  # we take the first name before a space to split
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


def check_aligned(seq_dict: dict[str: str]) -> bool:
    """
    checks if sequences in a sequence dictionary are the same length
    """
    if len(set(len(s) for s in seq_dict.values())) > 1:
        return False
    return True


def parse_phylip(path: str):
    """
    parse relaxed PHYLIP from file. Will fail with interleaved
    """
    with open(path, "r", encoding="utf-8") as inf:
        next(inf, (None, None))
        header = seq = ""
        for line in inf:
            line = line.strip()
            if bool(re.search(r"\s", line)):
                line = line.split()
                if header:
                    yield header, seq
                    header = line[0]
                    seq = line[1]
                else:
                    header = line[0]
                    seq = line[1]
            else:
                seq += line.strip()
        yield header, seq


def parse_phylip_str(phy_str: str):
    """
    parse relaxed PHYLIP from string. Will fail with interleaved
    """
    lines = phy_str.strip().split("\n")
    header = seq = ""
    for line in lines[1:]:
        line = line.strip()
        if bool(re.search(r"\s", line)):  # sequence header
            line = line.split()
            if header:
                yield header, seq
                header = line[0]
                seq = line[1]
            else:
                header = line[0]
                seq = line[1]
        else:
            seq += line.strip()
    yield header, seq


def get_phylip_str(seq_dict: dict[str: str]) -> str:
    """
    writes a PHYLIP-formatted string from an aligned sequence dictionary {header: sequence}. First
    allows a particular sequence to be placed at the start
    """
    if not check_aligned(seq_dict):
        raise ValueError("sequences are not aligned, write to FASTA instead")
    nseq = len(seq_dict)
    seql = set(len(v) for v in seq_dict.values()).pop()
    out = ""
    out += f" {nseq} {seql}\n"
    for header, seq in seq_dict.items():
        out += f"{header}    {seq}\n"
    return out


def parse_fasta(path: str):  # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
    """
    given a path tries to parse a fasta file. Returns an iterator which yields a (name, sequence) 
    tuple
    """
    with open(path, "r", encoding="utf-8") as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_fasta_str(fa_str: str):
    """
    Given a str representing the content of a FASTA-formatted file, return an iterator yielding 
    a (name, sequence tuple)
    """
    name = sequence = ""
    for line in fa_str.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name:
                yield name, sequence.upper()
            name = line[1:]
            sequence = ""
            continue
        sequence += line
    if name and sequence:
        yield name, sequence.upper()


def parse_partition_file(path: str):
    """
    given a path tries to parse a partition file. Returns an iterator yielding a
    (name, (start, end)) tuple
    """
    with open(path, "r", encoding="utf-8") as inf:
        for line in inf:
            line = line.strip().split(" ")
            name = line[1]
            sites = line[-1]
            start, end = [int(x) for x in sites.split("\\")[0].split("-")]
            yield name, (start, end)


def col_dict_to_seq_dict(col_dict: dict[int: dict[str: str]]) -> dict[str: str]:
    """
    function to return column dictionary to sequence dictionary
    """
    out = {}
    for col in col_dict.values():
        for header, char in col.items():
            try:
                out[header] += char
            except KeyError:
                out[header] = char
    return out


def get_fasta_str(seq_dict: dict[str: str]) -> str:
    """
    writes sequence dictionary to multiline string
    """
    ret = ""
    for header, seq in seq_dict.items():
        ret += f">{header}\n{seq}\n"
    return ret


def write_fasta(seq_dict: dict, out_file: str):
    """
    writes a fasta file from sequence dictionary
    """
    with open(out_file, "w", encoding="utf-8") as outf:
        outf.write(get_fasta_str(seq_dict=seq_dict))


def write_clustalw_conservation(seq_dict: dict) -> str:
    """
    write clustalw-style conservation string
    """
    strong = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
    weaker = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NEQHRK", "FVLIM", "HFY"]
    cols = get_columns(seq_dict)
    cons_str = ""
    for col in cols.values():
        unique_chars = set(col.values())
        if len(unique_chars) == 1:
            cons_str += "*"
        elif any(unique_chars <= set(group) for group in strong):
            cons_str += ":"
        elif any(unique_chars <= set(group) for group in weaker):
            cons_str += "."
        else:
            cons_str += " "
    return cons_str


def write_clustalw(seq_dict: dict) -> str:
    """
    write clustalw formatted alignment from seq_dict
    """
    if len(set(len(s) for s in seq_dict.values())) > 1:
        raise ValueError("sequences are not aligned!")
    conservation = write_clustalw_conservation(seq_dict)
    chars = len(list(seq_dict.values())[0])
    tot = 0
    out_str = ""
    out_str += "CLUSTAL W (X.XX)\n\n\n"
    longest_header = max(len(h) for h in seq_dict) + 6
    while tot < chars:
        for header, seq in seq_dict.items():
            pad = longest_header - len(header)
            if tot + 60 < chars:
                out_str += f"{header}" + " "*pad + seq[tot:tot+60] + "\n"
            else:
                out_str += f"{header}" + " "*pad + seq[tot:chars] + "\n"
        if tot + 60 < chars:
            out_str += " "*longest_header + conservation[tot:tot+60] + "\n\n"
        else:
            out_str += " "*longest_header + conservation[tot:chars] + "\n\n"
        tot += 60
    return out_str


def get_columns(seq_dict: dict) -> dict[int: dict[str: str]]:
    """
    takes a dictionary of an alignment (key: name, value: sequence) and returns a dictionary of 
    columns (key: position, value: dict{name: state})
    """
    col_dict = {}
    pos = 0
    for header, seq in seq_dict.items():
        for char in seq:
            try:
                col_dict[pos][header] = char
            except KeyError:
                col_dict[pos] = {}
                col_dict[pos][header] = char
            pos += 1
        pos = 0
    return col_dict


def get_site_specific_frequencies(seq_dict: dict, smooth: bool=True) -> dict:
    """
    takes a dictionary from get_columns and returns a dictionary of site-specific frequencies
    """
    aa = "ARNDCQEGHILKMFPSTWYV"
    col_dict = get_columns(seq_dict)
    freq_dict = {}  # key is pos, value is np array of freqs
    for pos, column in col_dict.items():
        counts = Counter(char for char in column.values() if char in set(aa))
        tot = sum(counts.values())
        if smooth:  # smoothing
            freqs = np.fromiter((counts[char] + 0.5 for char in aa), dtype=float) / (tot +
                                                                                        10)
            # Jeffrey's prior
        else:  # no smoothing
            freqs = np.fromiter((counts[char] for char in aa), dtype=float) / tot
        freq_dict[pos] = freqs
    return freq_dict


def write_site_specific_frequencies(site_freqs: dict, outpath: str="site_frequencies.tsv"):
    """
    write dictionary of site-specific frequencies to a tab-separated file of pos\tcomma-sep freqs
    """
    with open(outpath, "w", encoding="utf-8") as sff:
        for pos, freq in site_freqs.items():
            sff.write(f"{pos}\t")
            sff.write(",".join([f"{f:.4f}" for f in freq]) + "\n")


def calc_gene_frequencies(seq_dict: dict) -> np.array:
    """
    calculates the frequencies of amino acid character states in sequence dictionary, ignoring gaps
    """
    aa = "ARNDCQEGHILKMFPSTWYV"
    counts = Counter(char for seq in seq_dict.values() for char in seq if
                     char in set(aa))
    tot = sum(counts.values())
    frequencies = np.fromiter((counts[char] for char in aa), dtype=float) / tot
    return frequencies


def parse_site_frequencies_file(inf: str) -> dict:
    """
    parse a TSV of pos\tfreqs, where freqs is comma-separated
    """
    frequencies = {}
    with open(inf, "r", encoding="utf-8") as sff:
        for line in sff:
            line = line.strip().split("\t")
            frequencies[int(line[0])] = np.fromiter(line[1].split(","),
                                                    dtype=float)
    return frequencies


def parse_paml_rates(path: str) -> dict:
    """takes path to a paml 'rates' file containing estimated posterior mean 
    rate per site and returns a dictionary of {site: rate}"""
    out = {}
    with open(path, "r", encoding='utf-8') as inf:
        going = False
        reading = False
        for line in inf:
            if line.strip().startswith("Site"):
                going = True
                continue
            if line == "\n" and going and not reading:
                reading = True
                continue
            if line == "\n" and going and reading:
                break

            if going and reading:
                line = line.strip().split()
                out[int(line[0]) - 1] = float(line[-2])  # for 0-indexing
    return out


def insert_gaps_by_seq(ref:dict, query:dict):
    """clones gaps from a reference sequence into a query sequence, for example
    to insert gaps in simulated sequences to match an empirical alignment"""
    out = {}
    for header, seq in query.items():
        out[header] = ""
        try:
            for n, _ in enumerate(seq):
                if ref[header][n] == "-":
                    out[header] += "-"
                else:
                    out[header] += query[header][n]
        except KeyError as e:
            raise KeyError(header + " not in reference data, skipping\n") from e
    return out


if __name__ == "__main__":
    seqs = dict(parse_fasta(
        "DODAa_combined_no_og_strict_for_synth.cds.fa.nostop.name.noF.best.fas.trans"))
    sf = get_site_specific_frequencies(seqs)
    print(sf)
