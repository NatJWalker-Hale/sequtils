#! /usr/bin/python3

import sys
from Bio import Seq

STOP = ["TAA", "taa", "TGA", "tga", "TAG", "tag"]  # for case insensitivity
START = ["ATG", "atg"]  # for case insensitivity


def find_orf(seq, any_start=True, overlaps=True):
    i = 0
    currorf = ""
    orflist = []
    going = False
    while i < len(seq):
        codon = seq[i:i+3]
        codon.replace("-", "N")
        # first, establish start of the ORF
        if not going:
            if any_start:
                currorf = codon
                going = True
            elif codon in START:
                currorf = codon
                going = True
        # second, extend the ORF
        else:
            if codon not in STOP:
                if i+3 >= len(seq):
                    if len(codon) == 3:
                        currorf += codon
                        orflist.append(currorf)
                    else:
                        orflist.append(currorf)
                else:
                    currorf += codon
                    if overlaps:
                        orflist.extend(find_orf(seq[i:], any_start=False,
                                                overlaps=False))
            elif codon in STOP:
                currorf += codon
                orflist.append(currorf)
                going = False
                currorf = ""
        i += 3
    return orflist


def find_orf_all(seq):
    """Given a string representing a coding sequence, find all ORFs for
    all frames on the + and - strand"""
    all_orfs = []
    all_orfs.extend(find_orf(seq))
    all_orfs.extend(find_orf(seq[1:]))
    all_orfs.extend(find_orf(seq[2:]))
    revc = Seq.reverse_complement(seq)
    all_orfs.extend(find_orf(revc))
    all_orfs.extend(find_orf(revc[1:]))
    all_orfs.extend(find_orf(revc[2:]))
    return list(set(all_orfs))


def parse_fasta(path):  # courtesy of Jonathan Chang
    # https://gist.github.com/jonchang/6471846
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path) as handle:
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


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python "+sys.argv[0]+" cds")
        sys.exit()

    fa_dict = dict([x for x in parse_fasta(sys.argv[1])])
    orf_dict = {}
    for key, value in fa_dict.items():
        orfs = find_orf_all(value)
        orflist = []
        for orf in orfs:
            orflist.append((orf, Seq.translate(orf)))
        orf_dict[key] = sorted(orflist, reverse=True, key=lambda x: len(x[0]))

    for k, v in orf_dict.items():
        count = 0
        for i in v:
            print(">" + str(k) + "_ORF" + str(count) + "_cds")
            print(i[0])
            print(">" + str(k) + "_ORF" + str(count) + "_pep")
            print(i[1])
            count += 1
