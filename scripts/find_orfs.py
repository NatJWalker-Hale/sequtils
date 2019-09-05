import sys
from Bio import Seq

STOP = ["TAA","taa","TGA","tga","TAG","tag"] # for case insensitivity 
START = ["ATG","atg"] # for case insensitivity 

def find_orf(seq,ignore_overlap=False):
    """Given a string representing a single frame of a coding sequence from + or -
    strand, return all ORFs, including overlaps"""
    i = 0
    orfnum = 0
    currorf = ""
    orflist = []
    start = False
    while i < len(seq):
        codon = str(seq[i:i+3])
        if codon in START:
            if start:
                currorf += codon
                if not ignore_overlap:
                    orflist.extend(find_orf(seq[i:],ignore_overlap=False))
            else:
                start = True
                currorf = codon
        elif start and codon not in STOP:
            currorf += codon
        elif start:
            currorf += codon
            orflist.append(currorf)
            start = False
            currorf = ""
        i += 3
    return orflist

def find_orf_all(seq):
    """Given a string representing a coding sequence, find all ORFs for 
    all frames on the + and - strand"""
    all_orfs = []
    all_orfs.extend(find_orf(seq[0:]))
    all_orfs.extend(find_orf(seq[1:]))
    all_orfs.extend(find_orf(seq[2:]))
    revc = Seq.reverse_complement(seq)
    all_orfs.extend(find_orf(revc[0:]))
    all_orfs.extend(find_orf(revc[1:]))
    all_orfs.extend(find_orf(revc[2:]))
    return list(set(all_orfs))

def parse_fasta(path): # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
    with open(path) as handle:
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
            orflist.append((orf,Seq.translate(orf)))
        orf_dict[key] = sorted(orflist, reverse=True, key=lambda x: len(x[0]))
    print(orf_dict)
    





