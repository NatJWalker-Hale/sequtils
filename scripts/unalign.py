import sys

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
        print("usage: python "+sys.argv[0]+" aln_fasta")
        sys.exit()

    for x in parse_fasta(sys.argv[1]):
        print(">"+x[0])
        print(x[1].replace("-",""))