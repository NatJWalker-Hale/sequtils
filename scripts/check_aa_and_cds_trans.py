import sys
import Bio
from Bio import Seq
import warnings
warnings.filterwarnings("error")


STOP = ["TAA", "taa", "TGA", "tga", "TAG", "tag"]
START = ["ATG", "atg"]


def find_orf(seq, any_start=True, overlaps=True):
    i = 0
    currorf = ""
    orflist = []
    going = False
    while i < len(seq):
        codon = seq[i:i+3]
        # if "N" in codon or "n" in codon or "-" in codon:
        #     codon = "NNN"  # make all semi-ambiguous codons ambiguous
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


def parse_fasta(path):
    # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
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
    if len(sys.argv) != 3:
        print("usage: python "+sys.argv[0]+" cds aa")
        sys.exit()

    """This script extends my ORF finder script to take a cds and
    corresponding aa translation, calculate all ORFs, and find the
    CDS that matches the aa for that sequence. It is designed for
    the case where one might have a mix of cds and transcripts and
    doesn't want to manually figure out the ORFs."""

    cds_dict = dict([x for x in parse_fasta(sys.argv[1])])
    aa_dict = dict([x for x in parse_fasta(sys.argv[2])])
    # verify that files match
    # cds_seq = [x for x in cds_dict.keys()]
    # aa_seq = [x for x in aa_dict.keys()]
    # if set(cds_seq) != set(aa_seq):
    #     print("These files do not contain the same set of sequences!")
    #     sys.exit()

    # getting ORFs
    # orf_dict = {}
    # for key, value in cds_dict.items():
    #     orfs = find_orf_all(value)
    #     orflist = []
    #     for orf in orfs:
    #         orflist.append((orf, Seq.translate(orf.upper(), gap="N")))
    #     orf_dict[key] = sorted(orflist, reverse=True, key=lambda x: len(x[0]))
    # LEAVING THIS AS EVIDENCE OF HOW STUPID I ONCE WAS

    # check translations against aa
    for key, value in cds_dict.items():
        # first check against aa
        if key not in aa_dict.keys():
            sys.stderr.write("%s not in AAs, printing original CDS\n" % key)
            print(">"+key)
            print(cds_dict[key])
            continue
        match = False
        modmatch = False
        try:
            trans = Seq.translate(value.upper(), gap="N")
        except Bio.BiopythonWarning as w:
            if w.args[0].startswith("Partial codon"):
                sys.stderr.write("Sequence %s is not a multiple of 3. Trimming.\n" % key)
                # check which side has mismatch
                start = Seq.translate(value[:15])
                if start == aa_dict[key][:5]:  # l matches, r must be problem
                    trim = len(value) % 3
                    value = value[:-trim]
                    trans = Seq.translate(value.upper(), gap="N")
                else:  # assume l is problem
                    trim = len(value) % 3
                    value = value[trim:]
                    trans = Seq.translate(value.upper(), gap="N")
            else:
                sys.stderr.write("Warning: %s. Skipping.\n" % w.args[0])
                continue
        if trans.rstrip("*") == aa_dict[key].rstrip("*"):
            match = True
            modmatch = True
            print(">"+key)
            print(value)
        if not match:
            orfs = find_orf_all(value)
            orflist = []
            for orf in orfs:
                orflist.append((orf, Seq.translate(orf.upper(), gap="N")))
            for f in orflist:
                if f[1].rstrip("*") == aa_dict[key].rstrip("*"):
                    print(">"+key)
                    print(f[0])
                    # if f[0][-3:] in STOP:
                    #     print(f[0][:-6])
                    # else:
                    #     print(f[0][:-3])
                    # sys.stderr.write("Warning: %s matched up to final amino acid! Printing modified inferred CDS\n" % key)
                    # NOT A CLUE
                    modmatch = True
        if not modmatch:  # this will also catch cases where Biopython translates ambiguous codons but peptide file is X
            print(">"+key)
            print(cds_dict[key])
            sys.stderr.write("Warning: no match found for "+key+"! Printing original CDS\n")
