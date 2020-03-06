import sys 
import os
import argparse

STOP = ["TAA","taa","TGA","tga","TAG","tag"] # for case insensitivity 

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

def remove_stop(seq):
    i = 0
    while i < len(seq):
        codon = seq[i:i+3]
        if codon in STOP:
            if i+3 < len(seq):
                seq = seq[:i]+"NNN"+seq[i+3:]
            else:
                seq = seq[:-3]
        i += 3
    return seq

if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--sequence",help="Codon sequence file.")
    args = parser.parse_args()

    seq_dict = dict([x for x in parse_fasta(args.sequence)])
    for key, value in seq_dict.items():
        if len(value) % 3 != 0:
            print("Sequence length not evenly divisible by 3. Check file.")
            sys.exit()
        print(">"+key)
        print(remove_stop(value))


